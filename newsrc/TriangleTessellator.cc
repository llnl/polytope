//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include "TriangleTessellator.hh"
#include "EdgeUtils.hh"

#include <iostream>

#include "polytope_internal.hh" // Pulls in POLY_ASSERT
#include "Shapes.hh"
#include "QuantPLC.hh"
#include "GeomUtils.hh"
#include "Intersections.hh"
#include "Cell.hh"
#include "Clipping2D.hh"

#define TRILIBRARY
#define ANSI_DECLARATORS
#define CDT_ONLY

extern "C" {
#include "triangle.h"
}

namespace polytope {

namespace {
using IntType = HashKey<2>::IntType;
void printpoint(const Quantizer<2>& Q,
                const Point2<IntType>& point) {
  auto qp = Q.dequantize(point);
  std::cout << qp << std::endl;
}
}

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
void
TriangleTessellator::
tessellateQuantized(const QuantPLC<2>& qplc,
                    QuantizedTessellation& result) const {
  // Type aliases
  using IntType = typename QuantTessellation<2>::IntType;
  using RealPoint = Point2<double>;
  using IntPoint = typename QuantTessellation<2>::IntPoint;
  const Quantizer<2>& Q = result.m_Q;
  // Get the generators
  std::vector<double> generators = flattenCoords(result.getRealPoints());
  const auto numGenerators = generators.size()/2;

  // Build tessellation data structures (common for both cases)
  result.m_cells.resize(numGenerators);

  // Map IntPoint coordinates to node indices for deduplication
  std::map<IntPoint, int> node2id;

  // Map canonical edges to face indices for orientation tracking
  edge::EdgeToFaceMap edgeToFace;

  // Add nodes for the box extent
  std::vector<IntPoint> box = shapes::createBoxPoints(Q.minBound, Q.maxBound);
  std::map<shapes::BoxSide, unsigned> cornerIndices;
  shapes::BoxSides sides;
  for (unsigned i = 0; i < 4; i++) {
    const auto n = result.m_nodes.size();
    cornerIndices[sides.corner(i)] = n;
    node2id[box[i]] = n;
    result.m_nodes.push_back(box[i]);
  }
  // Special case: Triangle cannot handle 2 generators, so construct Voronoi manually
  if (numGenerators == 2) {
    Clip2D<IntType> clipper;
    // Compute perpendicular bisector between the two generators
    clipper.gen0 = result.m_points[1];
    clipper.gen1 = result.m_points[0];
    clipper.inf0 = true;
    clipper.inf1 = true;
    if (clipper.doClipping(Q)) {
      return;
    }

    edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.m_nodes);
    std::vector<edge::Edge> localEdges;
    shapes::walkBoxEdges(clipper.firstSide, clipper.curSide, cornerIndices,
                         curEdge, curEdge.first, false, localEdges);
    for (const auto& cedge : localEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
      result.m_cells[1].push_back(signedFaceIndex);
    }
    curEdge = edge::updateNodeMap(clipper.p1, clipper.p0, node2id, result.m_nodes);
    localEdges.clear();
    shapes::walkBoxEdges(clipper.curSide, clipper.firstSide, cornerIndices,
                         curEdge, curEdge.second, false, localEdges);
    for (const auto& cedge : localEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
      result.m_cells[0].push_back(signedFaceIndex);
    }
    return;
  }

  // Prepare Triangle input structure
  triangulateio in, out;
  in.segmentlist = nullptr;
  in.segmentmarkerlist = nullptr;
  in.holelist = nullptr;
  in.numberofholes = 0;
  in.numberofsegments = 0;
  in.numberofpoints = numGenerators;
  in.pointlist = new RealType[2*in.numberofpoints];
  std::copy(generators.begin(), generators.end(), in.pointlist);
  in.numberofpointattributes = 0;
  in.pointmarkerlist = nullptr;
  in.trianglelist = nullptr;
  in.triangleattributelist = nullptr;
  in.neighborlist = nullptr;
  in.regionlist = nullptr;
  in.numberofregions = 0;
  out.pointlist = 0;
  out.pointattributelist = 0;
  out.pointmarkerlist = 0;
  out.trianglelist = 0;
  out.triangleattributelist = 0;
  out.neighborlist = 0;
  out.segmentlist = 0;
  out.segmentmarkerlist = 0;
  out.edgelist = 0;
  out.edgemarkerlist = 0;
  // Normal case: use Triangle for 3+ generators
  triangulate((char*)"Qzn", &in, &out, 0);
  unsigned ntri = out.numberoftriangles;

  std::vector<std::set<unsigned>> gen2tri(numGenerators);
  std::vector<Point2<double>> centers;
  centers.reserve(ntri);

  // Extract the circumcenters of triangles (these become Voronoi vertices)
  for (auto i = 0; i < ntri; ++i) {
    int ia = out.trianglelist[3*i];
    int ib = out.trianglelist[3*i+1];
    int ic = out.trianglelist[3*i+2];
    auto a = result.m_points[ia].template type_cast<double>();
    auto b = result.m_points[ib].template type_cast<double>();
    auto c = result.m_points[ic].template type_cast<double>();
    Point2<double> rcen = circumcenter(a, b, c);
    centers.push_back(rcen);
    // Use original coordinates and quantize them but keep them as doubles
    // Point2<double> a(generators[2*ia], generators[2*ia+1]);
    // Point2<double> b(generators[2*ib], generators[2*ib+1]);
    // Point2<double> c(generators[2*ic], generators[2*ic+1]);
    // Point2<double> rcen = circumcenter(a, b, c);
    // Point2<double> rcenq = rcen.convertXi<double, double>(Q.m_xlo_o, Q.m_dx_o);
    // centers.push_back(rcenq);
    gen2tri[ia].insert(i);
    gen2tri[ib].insert(i);
    gen2tri[ic].insert(i);
  }

  // Process each generator to build its Voronoi cell
  // TODO: Retain edges for other generators to eliminate redundant calculations
  for (int cellIndex = 0; cellIndex < numGenerators; ++cellIndex) {
    // Walk edges around this generator point
    auto genit = gen2tri[cellIndex].begin();
    int curTri = *genit;
    bool ccwDir = true;
    std::vector<edge::Edge> localEdges;
    // List of sides associated with clipped nodes
    std::vector<std::pair<int, int>> clippedNodeSides;
    // Walk the edges, if there is an infinite edge in the
    // CW direction of this cell, start there
    for (auto it : gen2tri[cellIndex]) {
      int v0 = out.trianglelist[3*it];
      int v1 = out.trianglelist[3*it+1];
      // Find which vertex is the generator
      int localIndex = (v0 == cellIndex) ? 0 : (v1 == cellIndex) ? 1 : 2;
      int prevSide = (localIndex + 2)%3;
      bool curBound = Q.inQBounds(centers[it]);
      int prevTri = out.neighborlist[3*it+prevSide];
      bool prevBound = true;
      if (prevTri != -1) {
        prevBound = Q.inQBounds(centers[prevTri]);
      }
      if (curBound && (prevTri == -1 || !prevBound)) {
        curTri = it;
      }
    }
    int startTri = curTri;
    do {
      POLY_ASSERT2(curTri >= 0, "Cannot have negative curTri");
      int v0 = out.trianglelist[3*curTri];
      int v1 = out.trianglelist[3*curTri+1];
      int v2 = out.trianglelist[3*curTri+2];
      int tri[3] = {v0, v1, v2};

      // Find which vertex is the generator
      int localIndex = (v0 == cellIndex) ? 0 : (v1 == cellIndex) ? 1 : 2;

      int ccwSide = (localIndex + 1)%3;
      int cwSide = (localIndex + 2)%3;
      int localSide = (ccwDir) ? ccwSide : cwSide;
      int nextTri = out.neighborlist[3*curTri+localSide];
      Clip2D<IntType> clipper;
      int otherGen = tri[3 - localIndex - localSide];
      int indx1 = (ccwDir) ? cellIndex : otherGen;
      int indx2 = (ccwDir) ? otherGen : cellIndex;
      clipper.gen0 = result.m_points[indx1];
      clipper.gen1 = result.m_points[indx2];
      clipper.rp0 = centers[curTri];
      if (nextTri == -1) {
        clipper.inf1 = true;
        auto thirdPoint = result.m_points[tri[localSide]];
        clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1, thirdPoint);
      } else {
        clipper.rp1 = centers[nextTri];
        clipper.normalRay = pointDirection<IntType>(clipper.rp0, clipper.rp1);
      }
      if (clipper.doClipping(Q)) {
        if (nextTri == -1) {
          if (!ccwDir) break;
          ccwDir = false;
          nextTri = startTri;
        }
        curTri = nextTri;
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
      if (!ccwDir) {
        std::swap(clipper.p0, clipper.p1);
        std::swap(clipper.inf0, clipper.inf1);
        std::swap(startSide, endSide);
      }
      edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1,
                                               node2id, result.m_nodes);
      if (curEdge.first == curEdge.second) {
        if (nextTri == -1) {
          if (!ccwDir) break;
          ccwDir = false;
          nextTri = startTri;
        } else {
          curTri = nextTri;
          if (curTri == startTri) break;
        }
        continue;
      }
      if (ccwDir) {
        localEdges.push_back(curEdge);
        clippedNodeSides.push_back(std::make_pair(startSide, endSide));
      } else {
        localEdges.insert(localEdges.begin(), curEdge);
        clippedNodeSides.insert(clippedNodeSides.begin(),
                                std::make_pair(startSide, endSide));
      }
      if (nextTri == -1) {
        if (!ccwDir) break;
        curTri = startTri;
        ccwDir = false;
      } else {
        curTri = nextTri;
        if (curTri == startTri) break;
      }
    } while (true);
    std::vector<edge::Edge> finalEdges = shapes::closeClippedEdges(localEdges, clippedNodeSides, cornerIndices);
    // Remove collinear points from the edge loop
    removeCollinear(finalEdges, result.m_nodes);
    // Create faces and add to cell
    for (const auto& cedge : finalEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
      result.m_cells[cellIndex].push_back(signedFaceIndex);
    }
  }
  // Clean up Triangle memory
  delete[] in.pointlist;
  // Note: Triangle allocates out.* arrays, but they're cleaned up by Triangle internally
}

} //end polytope namespace
