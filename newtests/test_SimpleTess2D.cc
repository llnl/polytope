// Comprehensive unit tests for simple 2D tessellations
//

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "polytope.hh"
#include "Intersections.hh"
#include "Point.hh"
#include "polytope_test_utilities.hh"
#include "Boundary2D.hh"
#include "QuantTessellation.hh"
#include "BoostTessellator.hh"
#include "TriangleTessellator.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace polytope;
using namespace std;

namespace {

using CoordType = typename HashKey<2>::IntType;;
using IntPoint = Point2<CoordType>;

//------------------------------------------------------------------------------
// Simple square
//------------------------------------------------------------------------------
void testSquare(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Square ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  vector<double> points = {-0.05, -0.05, 0.05, -0.05, 0.05, 0.05, -0.05, 0.05};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Simple square with diagonal in top right
//------------------------------------------------------------------------------
void testCutSquare(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Cut Square ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  vector<double> points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Cut square with a triangle hole
//------------------------------------------------------------------------------
void testCutSquareHole(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Cut Square Triangle Hole ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  vector<double> points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
  // Add a triangle hole
  vector<double> newPoints = {0.3, -0.4, 0.2, 0.4, 0.2, -0.4};
  auto Nf = newPoints.size()/2;
  auto N = boundary.mPLCpoints.size()/2;
  copy(newPoints.begin(), newPoints.end(), back_inserter(boundary.mPLCpoints));
  boundary.mPLC.holes = vector<vector<vector<int>>>(1);
  boundary.mPLC.holes[0].resize(Nf);
  for (int i = 0; i < Nf; ++i) {
    unsigned fbegin = N + i;
    unsigned fend = N + (i+1)%Nf;
    boundary.mPLC.holes[0][i].resize(2);
    boundary.mPLC.holes[0][i][0] = fbegin;
    boundary.mPLC.holes[0][i][1] = fend;
  }
  //boundary.finalize();
  Quantizer<2> Q(boundary.mQ);
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Diamond tessellation
//------------------------------------------------------------------------------
void testDiamond(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Diamond ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  double lov = -0.4;
  double hiv = 0.4;
  vector<double> points = {lov, 0, hiv, 0, 0, lov, 0, hiv};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Obtuse triangle generators
//------------------------------------------------------------------------------
void testObtuse(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Obtuse Triangle ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.mDiff = 1.;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);

  vector<double> points = {0.67, -0.14,
                           0.91, 0.3,
                           0.49, -0.4};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Obtuse triangle generators
//------------------------------------------------------------------------------
void testClippedObtuse(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Clipped Obtuse Triangle ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(5);
  Quantizer<2> Q(boundary.mQ);
  vector<Point2<int>> ipoints = {{388767934, 243196867},
                                 {432427229, 322641370},
                                 {356518563, 195394647}};
  vector<double> points;
  for (auto& p : ipoints) {
    auto pp = Q.dequantize(p);
    points.push_back(pp[0]);
    points.push_back(pp[1]);
  }
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  QuantTessellation<2> quantMesh(Q, points);
  quantMesh.cullExternalPoints(QPLC);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Nearly external edge
//------------------------------------------------------------------------------
void testExternalEdge(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": External Edge ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  vector<double> points = {0.49, -0.1,
                           0.495, 0.,
                           0.485, 0.,
                           0.49, 0.1};
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  QuantTessellation<2> quantMesh(Q, points);
  quantMesh.cullExternalPoints(QPLC);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Only two generators
//------------------------------------------------------------------------------
void testTwoPoints(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Two Generators ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.mDiff = 1.4;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  vector<double> points = {1.0, 0.3,
                           0.5, -0.4};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    boost.tessellateQuantized(QPLC, quantMesh);
  } else {
    TriangleTessellator tri(Q);
    tri.tessellateQuantized(QPLC, quantMesh);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Clipped point
//------------------------------------------------------------------------------
void testClippedPoint(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Clipped Point ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.mDiff = 1.0;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  vector<double> points = {0.9, 0.3,
                           0.5, -0.4,
                           1.2, 0.9}; // This generator should be removed
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    quantMesh.cullExternalPoints(QPLC);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    quantMesh.cullExternalPoints(QPLC);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

//------------------------------------------------------------------------------
// Collinear
//------------------------------------------------------------------------------
void testCollinear(const int tnum, bool boostTess) {
  cout << "\n=== Test " << tnum << ": Collinear ===" << endl;
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.mDiff = 1.0;
  boundary.setDefaultBoundary(0);
  Quantizer<2> Q(boundary.mQ);
  vector<double> points = {0.1, 0.3, 0.7, 0.3, 0.4, 0.3};
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
  if (boostTess) {
    BoostTessellator boost(Q);
    quantMesh.cullExternalPoints(QPLC);
    boost.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
  } else {
    TriangleTessellator tri(Q);
    quantMesh.cullExternalPoints(QPLC);
    tri.tessellateQuantized(QPLC, quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, 0.0);
}

} // anonymous namespace

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  try {
    bool boost = true;
    for (int i = 0; i < 2; ++i) {
      int test = 1;
      testSquare(test++, boost);
      testCutSquare(test++, boost);
      testCutSquareHole(test++, boost);
      testDiamond(test++, boost);
      testObtuse(test++, boost);
      testClippedObtuse(test++, boost);
      testExternalEdge(test++, boost);
      testTwoPoints(test++, boost);
      testClippedPoint(test++, boost);
      testCollinear(test++, boost);
      boost = false;
    }

    cout << "\n=== ALL TESTS PASSED ===" << endl;
  } catch (const exception& e) {
    cout << "\n=== TEST FAILED WITH EXCEPTION ===" << endl;
    cout << e.what() << endl;
#ifdef POLYTOPE_ENABLE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
