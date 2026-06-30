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

void initTriangleData(triangulateio& in) {
  in.pointlist = nullptr;
  in.pointattributelist = nullptr;
  in.pointmarkerlist = nullptr;
  in.numberofpoints = 0;
  in.numberofpointattributes = 0;
  in.trianglelist = nullptr;
  in.triangleattributelist = nullptr;
  in.trianglearealist = nullptr;
  in.neighborlist = nullptr;
  in.numberoftriangles = 0;
  in.numberofcorners = 0;
  in.numberoftriangleattributes = 0;
  in.segmentlist = nullptr;
  in.segmentmarkerlist = nullptr;
  in.numberofsegments = 0;
  in.holelist = nullptr;
  in.numberofholes = 0;
  in.regionlist = nullptr;
  in.numberofregions = 0;
  in.edgelist = nullptr;
  in.edgemarkerlist = nullptr;
  in.numberofedges = 0;
  in.segmentlist = nullptr;
  in.segmentmarkerlist = nullptr;
  in.holelist = nullptr;
  in.numberofholes = 0;
  in.numberofsegments = 0;
  in.numberofedges = 0;
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
  const auto N = generators.size()/2;

  // Build tessellation data structures (common for both cases)
  result.m_cells.resize(N);

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
  // Prepare Triangle input structure
  triangulateio in, out;
  initTriangleData(in);
  initTriangleData(out);
  in.numberofpoints = N;
  in.pointlist = new RealType[2*in.numberofpoints];
  std::copy(generators.begin(), generators.end(), in.pointlist);
  unsigned ntri = 0;
  if (N > 2) {
    // Normal case: use Triangle for 3+ generators
    triangulate((char*)"Qzn", &in, &out, 0);
    ntri = out.numberoftriangles;
  }
  //-------------------------------------------------------------------
  // Special collinear or 2 generators cases
  //-------------------------------------------------------------------
  if (ntri == 0) {
    // Points are already ordered by hash so walk them in order and solve
    std::vector<std::vector<edge::Edge>> localEdges(N);
    // List of sides associated with clipped nodes
    std::vector<std::vector<std::pair<int, int>>> clippedNodeSides(N);
    for (int cellIndex = 0; cellIndex < N-1; ++cellIndex) {
      int nextPoint = cellIndex + 1;
      Clip2D<IntType> clipper;
      clipper.gen0 = result.m_points[cellIndex];
      clipper.gen1 = result.m_points[nextPoint];
      clipper.inf0 = true;
      clipper.inf1 = true;
      clipper.normalRay = outwardRay<IntType>(clipper.gen0, clipper.gen1);
      if (clipper.doClipping(Q)){
        continue;
      }
      int startSide = static_cast<int>(clipper.firstSide);
      int endSide = static_cast<int>(clipper.secondSide);
      edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.m_nodes);
      localEdges[cellIndex].push_back(curEdge);
      clippedNodeSides[cellIndex].push_back(std::make_pair(startSide, endSide));
      curEdge = edge::updateNodeMap(clipper.p1, clipper.p0, node2id, result.m_nodes);
      localEdges[nextPoint].push_back(curEdge);
      clippedNodeSides[nextPoint].push_back(std::make_pair(endSide, startSide));
    }
    for (int cellIndex = 0; cellIndex < N; ++cellIndex) {
      std::vector<edge::Edge> finalEdges = shapes::closeClippedEdges(localEdges[cellIndex], clippedNodeSides[cellIndex], cornerIndices);
      removeCollinear(finalEdges, result.m_nodes);
      for (const auto& cedge : finalEdges) {
        int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.m_faces, edgeToFace);
        result.m_cells[cellIndex].push_back(signedFaceIndex);
      }
    }
    return;
  }
  std::vector<std::set<unsigned>> gen2tri(N);
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
    gen2tri[ia].insert(i);
    gen2tri[ib].insert(i);
    gen2tri[ic].insert(i);
  }

  // Process each generator to build its Voronoi cell
  // TODO: Retain edges for other generators to eliminate redundant calculations
  for (int cellIndex = 0; cellIndex < N; ++cellIndex) {
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
      if (curTri == -1) {
        if (!ccwDir) break;
        curTri = startTri;
        ccwDir = false;
      }
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
      // gen0 should always be the current cell's generator
      clipper.gen0 = result.m_points[cellIndex];
      clipper.gen1 = result.m_points[otherGen];
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
      edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.m_nodes);
      if (curEdge.first == curEdge.second) {
        curTri = nextTri;
        continue;
      }
      localEdges.push_back(curEdge);
      clippedNodeSides.push_back(std::make_pair(startSide, endSide));
      curTri = nextTri;
    } while (curTri != startTri);
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
