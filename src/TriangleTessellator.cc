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
#include "VoronoiConstructor.hh"
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
  if (result.points.empty()) {
    return;
  }
  // Type aliases
  const auto& Q = Quantizer<2>::instance();
  // Get quantized generators cast as doubles and flattened
  std::vector<double> generators = flattenCoords(result.getRealQPoints());
  const auto N = generators.size()/2;

  // Build tessellation data structures (common for both cases)
  result.cells.resize(N);

  VoronoiConstructor constructor(result);

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
    // A single generator has no Voronoi primitive: its cell is the box.
    if (N == 1u) {
      BoxSides sides;
      for (unsigned i = 0; i < 4; ++i) {
        const auto n0 = constructor.cornerIndices.at(sides.corner(i));
        const auto n1 = constructor.cornerIndices.at(sides.corner((i + 1) % 4));
        result.cells[0].push_back(
            addOrientedEdge(n0, n1, result.faces, constructor.edgeToFace));
      }
      delete[] in.pointlist;
      return;
    }

    // Points are already ordered by hash.  Each adjacent pair produces an
    // unbounded Voronoi bisector; VoronoiConstructor clips, orients, closes,
    // and turns those primitives into cells.
    std::vector<std::vector<VoronoiPrimitive>> cellPrimitives(N);
    for (auto cellIndex = 0u; cellIndex < N-1; ++cellIndex) {
      cellPrimitives[cellIndex].emplace_back(cellIndex, cellIndex + 1);
      cellPrimitives[cellIndex + 1].emplace_back(cellIndex, cellIndex + 1);
    }
    for (auto cellIndex = 0u; cellIndex < N; ++cellIndex) {
      constructor.constructEdges(cellPrimitives[cellIndex], cellIndex);
    }
    delete[] in.pointlist;
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
    std::vector<VoronoiPrimitive> vps;
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
      VoronoiPrimitive vp(cellIndex, otherGen);
      vp.setV0(centers[curTri]);
      if (nextTri != -1) {
        vp.setV1(centers[nextTri]);
      } else {
        vp.setThirdPoint(tri[localSide]);
      }
      vps.push_back(vp);
      curTri = nextTri;
    } while (curTri != startTri);
    // Build voronoi from set of clipped edges
    constructor.constructEdges(vps, cellIndex);
  }
  // Clean up Triangle memory
  delete[] in.pointlist;
  // Note: Triangle allocates out.* arrays, but they're cleaned up by Triangle internally
}

} //end polytope namespace
