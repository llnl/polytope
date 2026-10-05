//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include "TriangleTessellator.hh"
#include "EdgeUtils.hh"

#include <iostream>
#include <utility>

#include "polytope_internal.hh" // Pulls in POLY_ASSERT
#include "GeomUtils.hh"
#include "Intersections.hh"

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
VoronoiAssembler<2>
TriangleTessellator::
tessellateQuantizedImpl(QuantizedTessellation& result) const {
  VoronoiAssembler<2> assembler(result);
  if (result.points.empty()) {
    return assembler;
  }
  // Get quantized generators cast as doubles and flattened
  std::vector<double> generators = flattenCoords(result.getRealQPoints());
  const auto N = generators.size()/2;

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
    assembler.degenerateAssembly();
    delete[] in.pointlist;
    return assembler;
  }
  std::vector<std::array<int, 3>> triangles;
  std::vector<std::array<int, 3>> neighbors;
  triangles.reserve(ntri);
  neighbors.reserve(ntri);

  for (auto i = 0u; i < ntri; ++i) {
    // triangles[i][v] is the generator index at Triangle local vertex
    // 'v'. The local vertex ordering must agree with neighbors below.
    int ia = out.trianglelist[3*i];
    int ib = out.trianglelist[3*i+1];
    int ic = out.trianglelist[3*i+2];
    triangles.push_back({ia, ib, ic});

    // neighbors[i][side] is the triangle across the side opposite
    // triangles[i][side]; Triangle uses -1 when that side is on the hull.
    neighbors.push_back({out.neighborlist[3*i],
                         out.neighborlist[3*i+1],
                         out.neighborlist[3*i+2]});
  }
  assembler.assembleDelaunay(triangles, neighbors);
  // Clean up Triangle memory
  delete[] in.pointlist;
  // Note: Triangle allocates out.* arrays, but they're cleaned up by Triangle internally
  return assembler;
}

} //end polytope namespace
