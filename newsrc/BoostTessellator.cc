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
  result.m_nodes.reserve(voronoi.num_vertices());
  result.m_faces.reserve(voronoi.num_edges());
  result.m_cells.resize(numGenerators);

  // Map IntPoint coordinates to our node indices (for deduplication)
  std::map<IntPoint, int> node2id;

  // Map canonical edges to face indices for oriented edge tracking
  edge::EdgeToFaceMap edgeToFace;

  // Add nodes for the box extent and keep track of their indices
  std::vector<IntPoint> box = shapes::createBoxPoints(Q.minBound, Q.maxBound);
  std::map<shapes::BoxCorner, unsigned> cornerIndices; // Ordered lower left and CCW
  for (unsigned i = 0; i < 4; i++) {
    const auto n = result.m_nodes.size();
    cornerIndices[static_cast<shapes::BoxCorner>(i)] = n;
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

    // Flags for handling infinite edges
    int firstClippedNode = -1;
    // Which box side each infinite edge intersects
    shapes::BoxSide firstBoxSide;
    do {
      bool flippedEdge = false;
      const VD::edge_type* nextEdge = edge->next();
      if (!edge->vertex0() && edge->vertex1()) {
        edge = edge->twin();
        flippedEdge = true;
      }
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // An edge is considered infinite if Boost provides a null pointer
      auto gindx1 = edge->cell()->source_index();//cellIndex;
      auto gindx2 = edge->twin()->cell()->source_index();

      Clip2D<IntType> clipper;
      clipper.gen0 = result.m_points[gindx1];
      clipper.gen1 = result.m_points[gindx2];
      clipper.normalRay = outwardRay(clipper.gen1, clipper.gen0);
      if (v0) {
        clipper.rp0 = Point2<double>(v0->x(), v0->y());
      } else {
        clipper.inf0 = true;
      }
      if (v1) {
        clipper.rp1 = Point2<double>(v1->x(), v1->y());
      } else {
        clipper.inf1 = true;
      }
      if (clipper.doClipping(Q)) {
        edge = nextEdge;
        continue;
      }
      if (flippedEdge) {
        std::swap(clipper.p0, clipper.p1);
        std::swap(clipper.inf0, clipper.inf1);
      }
      bool isInfinite = (clipper.inf0 || clipper.inf1);
      edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.m_nodes);

      if (curEdge.first == curEdge.second) {
        edge = nextEdge;
        continue;
      }
      // If both edges are infinite, make the start point n0
      if (clipper.inf0 && clipper.inf1) {
        firstClippedNode = curEdge.first;
        firstBoxSide = clipper.firstSide;
        clipper.inf0 = false;
      }
      // Care must be taken when dealing with infinite edges
      if (isInfinite && firstClippedNode < 0) {
        // First infinite edge, remember the intersecting edge
        firstBoxSide = clipper.curSide;
        firstClippedNode = (clipper.inf0) ? curEdge.first : curEdge.second;
        localEdges.push_back(curEdge);
      } else if (isInfinite) {
        // Second infinite edge, walk the exterior and accumulate corner points
        shapes::walkBoxEdges(firstBoxSide, clipper.curSide, cornerIndices,
                             curEdge, firstClippedNode, clipper.inf0, localEdges);
      } else {
        localEdges.push_back(curEdge);
      }

      edge = nextEdge;
    } while (edge != firstEdge);
    // Create faces and cells from local edges
    //removeCollinear(localEdges, result.m_nodes);
    for (const auto& cedge : localEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
      result.m_cells[cellIndex].push_back(signedFaceIndex);
    }
  }
}

} //end polytope namespace
