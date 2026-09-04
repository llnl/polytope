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
#include "Communicator.hh"

#define TRILIBRARY
#define ANSI_DECLARATORS
#define CDT_ONLY

#define REAL double
#define VOID void

extern "C" {
#include "triangle.h"
}

namespace polytope {

namespace {
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
tessellateQuantizedImpl(QuantizedTessellation& result) {
  // Type aliases
  const auto& Q = Quantizer<2>::instance();
  // Get quantized generators cast as doubles and flattened
  std::vector<double> generators = flattenCoords(result.getRealQPoints());
  const auto N = generators.size()/2;

  // Build tessellation data structures (common for both cases)
  result.cells.resize(N);

  // Map QuantizedPoint coordinates to node indices for deduplication
  std::map<QuantizedPoint<2>, int> node2id;

  // Map canonical edges to face indices for orientation tracking
  edge::EdgeToFaceMap edgeToFace;

  // Map generator pairs to edges
  edge::GenPairToEdgeDataMap genPairToEdge;

  // Add nodes for the box extent and keep track of their indices
  auto cornerIndices = addBoxPoints(Q, node2id, result.nodes);

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
  if (ntri == 0u) {
    // Points are already ordered by hash so walk them in order and solve
    std::vector<std::vector<edge::Edge>> localEdges(N);
    // List of sides associated with clipped nodes
    std::vector<std::vector<std::pair<int, int>>> clippedNodeSides(N);
    for (auto cellIndex = 0u; cellIndex < N-1; ++cellIndex) {
      int nextPoint = cellIndex + 1;
      Clip2D<QuantizedCoordinate<2>> clipper;
      clipper.gen0 = result.points[cellIndex];
      clipper.gen1 = result.points[nextPoint];
      clipper.inf0 = true;
      clipper.inf1 = true;
      clipper.normalRay = outwardRay<QuantizedCoordinate<2>>(clipper.gen0, clipper.gen1);
      if (clipper.doClipping()){
        continue;
      }
      int startSide = static_cast<int>(clipper.firstSide);
      int endSide = static_cast<int>(clipper.secondSide);
      edge::Edge curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.nodes);
      localEdges[cellIndex].push_back(curEdge);
      clippedNodeSides[cellIndex].push_back(std::make_pair(startSide, endSide));
      curEdge = edge::updateNodeMap(clipper.p1, clipper.p0, node2id, result.nodes);
      localEdges[nextPoint].push_back(curEdge);
      clippedNodeSides[nextPoint].push_back(std::make_pair(endSide, startSide));
    }
    for (auto cellIndex = 0u; cellIndex < N; ++cellIndex) {
      std::vector<edge::Edge> finalEdges = closeClippedEdges(localEdges[cellIndex], clippedNodeSides[cellIndex], cornerIndices);
      removeCollinear(finalEdges, result.nodes);
      for (const auto& cedge : finalEdges) {
        int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.faces, edgeToFace);
        result.cells[cellIndex].push_back(signedFaceIndex);
      }
    }
    return;
  }
  std::vector<std::set<unsigned>> gen2tri(N);
  std::vector<Point2<double>> centers;
  centers.reserve(ntri);

  // Extract the circumcenters of triangles (these become Voronoi vertices)
  for (auto i = 0u; i < ntri; ++i) {
    int ia = out.trianglelist[3*i];
    int ib = out.trianglelist[3*i+1];
    int ic = out.trianglelist[3*i+2];
    auto a = result.points[ia].template type_cast<double>();
    auto b = result.points[ib].template type_cast<double>();
    auto c = result.points[ic].template type_cast<double>();
    Point2<double> rcen = circumcenter(a, b, c);
    centers.push_back(rcen);
    gen2tri[ia].insert(i);
    gen2tri[ib].insert(i);
    gen2tri[ic].insert(i);
  }

  // Process each generator to build its Voronoi cell
  for (auto cellIndex = 0u; cellIndex < N; ++cellIndex) {
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
      int localIndex = (v0 == int(cellIndex)) ? 0 : (v1 == int(cellIndex)) ? 1 : 2;
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
      int localIndex = (v0 == int(cellIndex)) ? 0 : (v1 == int(cellIndex)) ? 1 : 2;

      int ccwSide = (localIndex + 1)%3;
      int cwSide = (localIndex + 2)%3;
      int localSide = (ccwDir) ? ccwSide : cwSide;
      int nextTri = out.neighborlist[3*curTri+localSide];
      int otherGen = tri[3 - localIndex - localSide];
      edge::GenPair gp = edge::orderPair(cellIndex, otherGen);
      edge::Edge curEdge;
      int startSide = -1;
      int endSide = -1;
      // Check if this pair of generators has already been computed
      auto cacheIt = genPairToEdge.find(gp);
      if (cacheIt != genPairToEdge.end()) {
        const edge::EdgeData& ed = cacheIt->second;
        curEdge = std::make_pair(ed.curEdge.second, ed.curEdge.first);
        startSide = ed.endSide;
        endSide = ed.startSide;
      } else {
        Clip2D<QuantizedCoordinate<2>> clipper;
        // gen0 should always be the current cell's generator
        clipper.gen0 = result.points[cellIndex];
        clipper.gen1 = result.points[otherGen];
        clipper.rp0 = centers[curTri];
        if (nextTri == -1) {
          clipper.inf1 = true;
          auto thirdPoint = result.points[tri[localSide]];
          clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1, thirdPoint);
        } else {
          clipper.rp1 = centers[nextTri];
          clipper.normalRay = pointDirection<QuantizedCoordinate<2>>(clipper.rp0, clipper.rp1);
        }
        if (clipper.doClipping()) {
          curTri = nextTri;
          continue;
        }
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
        curEdge = edge::updateNodeMap(clipper.p0, clipper.p1, node2id, result.nodes);
        if (curEdge.first == curEdge.second) {
          curTri = nextTri;
          continue;
        }
        genPairToEdge[gp] = {curEdge, startSide, endSide};
      }
      localEdges.push_back(curEdge);
      clippedNodeSides.push_back(std::make_pair(startSide, endSide));

      curTri = nextTri;
    } while (curTri != startTri);
    if (localEdges.size() > 2) {
      // Order the edges to have a consistent sequence of clipped node sides
      edge::orderClippedNodes(clippedNodeSides, localEdges);
    }
    std::vector<edge::Edge> finalEdges = closeClippedEdges(localEdges, clippedNodeSides, cornerIndices);
    // Remove collinear points from the edge loop
    removeCollinear(finalEdges, result.nodes);
    // Create faces and add to cell
    for (const auto& cedge : finalEdges) {
      int signedFaceIndex = edge::addOrientedEdge(cedge.first, cedge.second, result.faces, edgeToFace);
      result.cells[cellIndex].push_back(signedFaceIndex);
    }
    // Check for nearly duplicate nodes
    POLY_ASSERT2(!edge::hasNearDuplicates(result.points[cellIndex], node2id),
                 "Found nearly duplicate nodes.");
  }
  // Clean up Triangle memory
  delete[] in.pointlist;
  // Note: Triangle allocates out.* arrays, but they're cleaned up by Triangle internally
}

} //end polytope namespace
