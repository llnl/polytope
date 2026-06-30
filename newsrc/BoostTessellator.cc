//------------------------------------------------------------------------
// BoostTessellator
//------------------------------------------------------------------------
#include "BoostTessellator.hh"
#include "EdgeUtils.hh"

#include <iostream>

// Handy Boost stuff
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "polytope_internal.hh" // Pulls in POLY_ASSERT
#include "RegisterBoostPolygonTypes.hh"
#include "Shapes.hh"
#include "QuantPLC.hh"
#include "Clipping2D.hh"

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

namespace polytope {

namespace bp = boost::polygon;

namespace {
using IntType = HashKey<2>::IntType;
using PolygonWithHoles = bp::polygon_with_holes_data<IntType>;
using PolygonSet = bp::polygon_set_data<IntType>;
void printpoint(const Quantizer<2>& Q,
                const Point2<IntType>& point) {
  auto qp = Q.dequantize(point);
  std::cout << qp << std::endl;
}
void printPolygon(const Quantizer<2>& Q,
                  const PolygonWithHoles& polygon) {
  auto hole_points = bp::innerPoints(polygon);
  auto points = bp::outerPoints(polygon);
  std::cout << "v = [";
  for (const auto& point : points) {
    printpoint(Q, point);
  }
  std::cout << "vertices.append(v)" << std::endl;
  for (const auto& hole : hole_points) {
    std::cout << "hole " << std::endl;
    for (const auto& p : hole) {
      printpoint(Q, p);
    }
  }
}
void printCell(const Quantizer<2>& Q,
               const Cell<2, IntType>::CellType& pcell) {
  PolygonWithHoles polygon;
  bp::set_points(polygon, pcell.begin(), pcell.end());
  printPolygon(Q, polygon);
}
}

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//
//------------------------------------------------------------------------------
void
BoostTessellator::
tessellateQuantized(const QuantPLC<2>& qplc,
                    QuantizedTessellation& result) const {
  POLY_CONTRACT_VAR(qplc);
  // Type aliases
  using IntType = typename QuantTessellation<2>::IntType;
  using IntPoint = typename QuantTessellation<2>::IntPoint;
  using VD = boost::polygon::voronoi_diagram<RealType>;
  const Quantizer<2>& Q = result.m_Q;
  // Get the generators
  std::vector<IntPoint> generators = result.getIntPoints();
  const size_t numGenerators = generators.size();

  // Invoke the Boost.Voronoi diagram constructor
  VD voronoi;
  construct_voronoi(generators.begin(), generators.end(), &voronoi);

  // Build the tessellation data structures
  // In 2D: nodes are Voronoi vertices, faces are edges, cells are Voronoi cells
  result.m_nodes.reserve(voronoi.num_vertices()+4);
  result.m_faces.reserve(voronoi.num_edges());
  result.m_cells.resize(numGenerators);

  // Map IntPoint coordinates to our node indices (for deduplication)
  std::map<IntPoint, int> node2id;

  // Map canonical edges to face indices for oriented edge tracking
  edge::EdgeToFaceMap edgeToFace;

  // Add nodes for the box extent and keep track of their indices
  std::vector<IntPoint> box = shapes::createBoxPoints(Q.minBound, Q.maxBound);
  std::map<shapes::BoxSide, unsigned> cornerIndices; // Ordered lower left and CCW
  shapes::BoxSides sides;
  for (unsigned i = 0; i < 4; i++) {
    const auto n = result.m_nodes.size();
    cornerIndices[sides.corner(i)] = n;
    node2id[box[i]] = n;
    result.m_nodes.push_back(box[i]);
  }

  // Process each Voronoi cell
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin();
       cellItr != voronoi.cells().end();
       ++cellItr) {

    if (!cellItr->contains_point()) continue;

    const int cellIndex = cellItr->source_index();
    if (cellIndex >= numGenerators) continue;

    // Walk edges CCW around this cell
    const typename VD::edge_type* firstEdge = cellItr->incident_edge();
    const typename VD::edge_type* edge = firstEdge;

    // List of local edges
    std::vector<edge::Edge> localEdges;
    std::vector<std::pair<int, int>> clippedNodeSides;
    do {
      const VD::edge_type* nextEdge = edge->next();
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // An edge is considered infinite if Boost provides a null pointer
      // gen0 should always be the current cell's generator
      auto gindx1 = edge->cell()->source_index();
      auto gindx2 = edge->twin()->cell()->source_index();
      Clip2D<IntType> clipper;
      clipper.gen0 = result.m_points[gindx1];
      clipper.gen1 = result.m_points[gindx2];
      if (v0) {
        clipper.rp0 = Point2<double>(v0->x(), v0->y());
      } else {
        clipper.inf0 = true;
        clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1);
      }
      if (v1) {
        clipper.rp1 = Point2<double>(v1->x(), v1->y());
      } else {
        clipper.inf1 = true;
        clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1);
      }
      if (v1 && v0) {
        clipper.normalRay = pointDirection<IntType>(clipper.rp0, clipper.rp1);
      }
      if (clipper.doClipping(Q)) {
        edge = nextEdge;
        continue;
      }
      int startSide = -1;
      int endSide = -1;
      if (clipper.inf0) {
        startSide = static_cast<int>(clipper.firstSide);
      }
      if (clipper.inf1) {
        endSide = static_cast<int>(clipper.secondSide);
      }
      edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.m_nodes);

      if (curEdge.first == curEdge.second) {
        edge = nextEdge;
        continue;
      }
      clippedNodeSides.push_back(std::make_pair(startSide, endSide));
      localEdges.push_back(curEdge);

      edge = nextEdge;
    } while (edge != firstEdge);
    // Walk edges and clipped nodes to connect them
    std::vector<edge::Edge> finalEdges = shapes::closeClippedEdges(localEdges, clippedNodeSides, cornerIndices);
    // Create faces and cells from local edges
    removeCollinear(finalEdges, result.m_nodes);
    for (const auto& cedge : finalEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
      result.m_cells[cellIndex].push_back(signedFaceIndex);
    }
  }
}

} //end polytope namespace
