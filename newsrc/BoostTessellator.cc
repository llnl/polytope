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

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

namespace polytope {

namespace bp = boost::polygon;

namespace {
using IntType = HashKey<2>::IntType;
void printpoint(const Quantizer<2>& Q,
                const Point2<IntType>& point) {
  auto qp = Q.dequantize(point);
  std::cout << qp << std::endl;
}
}

void
BoostTessellator::
tessellateQuantized(const QuantPLC<2>& qplc,
                    QuantizedTessellation& result) const {
  POLY_CONTRACT_VAR(qplc);
  tessellateQuantized(result);
}

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//
// There are two ways to handle infinite edges provided by Boost.
// 1. If doClipping is disabled, infinite edges are extended to the bounding box
//    and edges are added to close the shape. If infinite edges span corners,
//    two edges are added around the corner.
// 2. If doClipping is enabled, extra generators (called guard generators) are
//    added at the maximum extents. Any infinite edges will be relative to these
//    generators. These generators are ignored altogether.
// During the first pass of tessellation, doClipping should be disabled. If
// tessellating is necessary during the clipping routine, doClipping should be
// enabled
//------------------------------------------------------------------------------
void
BoostTessellator::
tessellateQuantized(QuantizedTessellation& result, bool doClipping) const {
  // Type aliases
  using IntType = typename QuantTessellation<2>::IntType;
  using IntPoint = typename QuantTessellation<2>::IntPoint;
  using VD = boost::polygon::voronoi_diagram<RealType>;
  const Quantizer<2>& Q = result.m_Q;
  // Get the generators
  std::vector<IntPoint> generators = result.getIntPoints();
  const size_t numGenerators = generators.size();
  if (doClipping) {
   result.guardGenerators(generators);
  }

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
    int nodeIndx = 0;

    // List of local edges
    std::vector<edge::Edge> localEdges;

    // Flags for handling infinite edges
    int firstClippedNode = -1;
    // Which box side each infinite edge intersects
    shapes::BoxSide firstBoxSide;
    bool didskip = false;
    do {
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // TODO: Generalize infinite edge clipping
      // An edge is considered infinite if Boost provides a null pointer
      bool v0inf = (!v0) ? true : false;
      bool v1inf = (!v1) ? true : false;
      Point2<double> vp0, vp1;
      IntPoint p0, p1;
      bool bounds0 = true;
      bool bounds1 = true;
      shapes::BoxSide curSide;
      // There are two types of infinite points: Boost null pointer infinite points
      // and points provided by Boost that exceed our bounding box.
      // Each are handled differently
      if (v0) {
        vp0 = Point2<double>(v0->x(), v0->y());
        p0 = round<2, IntType>(vp0);
        bounds0 = Q.inQBounds(vp0);
      }
      if (v1) {
        vp1 = Point2<double>(v1->x(), v1->y());
        p1 = round<2, IntType>(vp1);
        bounds1 = Q.inQBounds(vp1);
      }
      // If both points are now infinite, we must extract the vertices
      auto gindx1 = cellIndex;
      auto gindx2 = edge->twin()->cell()->source_index();
      bool bothInf = false;
      if (v0inf && v1inf) {
        p1 = midPoint(result.m_points[gindx2], result.m_points[gindx1]);
        bothInf = true;
      } else if ((!bounds0 && !v1) || (!bounds1 && !v0) || (!bounds1 && !bounds0)) {
        didskip = true;
        edge = edge->next();
        continue;
      }
      if (v0inf) {
        // If vertex 0 is a Boost infinite vertex
        IntPoint outwardRay = normalRay(result.m_points[gindx2],
                                        result.m_points[gindx1]);
        clipInfiniteRay(p1, outwardRay, Q.minBound, Q.maxBound, p0, curSide);
        if (bothInf) {
          firstBoxSide = curSide;
        }
      } else if (!bounds0) {
        // If vertex 0 is outside our bounding box
        auto diff = vp0 - vp1;
        clipInfiniteRay(vp1, diff, Q.rminBound, Q.rmaxBound, vp0, curSide);
        p0 = round<2, IntType>(vp0);
        v0inf = true;
      }
      if (v1inf) {
        IntPoint outwardRay = normalRay(result.m_points[gindx1],
                                        result.m_points[gindx2]);
        clipInfiniteRay(p0, outwardRay, Q.minBound, Q.maxBound, p1, curSide);
      } else if (!bounds1) {
        auto diff = vp1 - vp0;
        clipInfiniteRay(vp0, diff, Q.rminBound, Q.rmaxBound, vp1, curSide);
        p1 = round<2, IntType>(vp1);
        v1inf = true;
      }
      bool isInfinite = (v0inf || v1inf);
      // Deduplicate by final IntPoint value, not by Boost vertex pointer
      auto it0 = node2id.find(p0);
      int n0;
      if (it0 == node2id.end()) {
        n0 = result.m_nodes.size();
        node2id[p0] = n0;
        p0.index = nodeIndx++;
        result.m_nodes.push_back(p0);
      } else {
        n0 = it0->second;
      }

      auto it1 = node2id.find(p1);
      int n1;
      if (it1 == node2id.end()) {
        n1 = result.m_nodes.size();
        node2id[p1] = n1;
        p1.index = nodeIndx++;
        result.m_nodes.push_back(p1);
      } else {
        n1 = it1->second;
      }

      if (n0 == n1) {
        edge = edge->next();
        continue;
      }
      edge::Edge curEdge = std::make_pair(n0, n1);
      // If both edges are infinite, make the start point n0
      if (bothInf) {
        firstClippedNode = n0;
        v0inf = false;
      }
      // Care must be taken when dealing with infinite edges
      if (isInfinite && firstClippedNode < 0) {
        // First infinite edge, remember the intersecting edge
        firstBoxSide = curSide;
        firstClippedNode = (v0inf) ? n0 : n1;
        localEdges.push_back(curEdge);
      } else if (isInfinite) {
        // Second infinite edge, walk the exterior and accumulate corner points
        shapes::walkBoxEdges(firstBoxSide, curSide, cornerIndices,
                             curEdge, firstClippedNode, v0inf, localEdges);        
      } else {
        localEdges.push_back(curEdge);
      }

      edge = edge->next();
    } while (edge != firstEdge);
    // Create faces and cells from local edges
    removeCollinear(localEdges, result.m_nodes);
    for (const auto& cedge : localEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
      result.m_cells[cellIndex].push_back(signedFaceIndex);
    }
  }
}

} //end polytope namespace
