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

void tests(const int tnum, bool boostTess) {
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D<double> boundary;
  boundary.mDiff = 1.;
  std::vector<double> points;
  std::string testname;
  switch (tnum) {
  case 1: // Square
    testname = "Square";
    boundary.setDefaultBoundary(0);
    points = {-0.05, -0.05, 0.05, -0.05, 0.05, 0.05, -0.05, 0.05};
    break;
  case 2: // Cut square
    testname = "Cut square";
    boundary.setDefaultBoundary(0);
    points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
    break;
  case 3: // Cut square with hole
    {
      testname = "Cut square with hole";
      boundary.setDefaultBoundary(0);
      points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
      // Add a triangle hole
      vector<double> newPoints = {0.6, -0.8, 0.4, 0.8, 0.4, -0.8};
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
    }
    break;
  case 4: // Diamond
    testname = "Diamond";
    boundary.setDefaultBoundary(0);
    points = {-0.4, 0, 0.4, 0, 0, -0.4, 0, 0.4};
    break;
  case 5: // Obtuse
    testname = "Obtuse triangle";
    boundary.setDefaultBoundary(0);
    points = {0.67, -0.14, 0.91, 0.3, 0.49, -0.4};
    break;
  case 6: // Obtuse with star
    testname = "Obtuse triangle clipped by star";
    boundary.setDefaultBoundary(5);
    points = {0.67, -0.14, 0.91, 0.3, 0.49, -0.4};
    break;
  case 7: // Mod pad obtuse with star
    testname = "Mod pad obtuse triangle clipped by star";
    boundary.m_pad = 0.5;
    boundary.setDefaultBoundary(5);
    points = {0.67, -0.14, 0.91, 0.3, 0.49, -0.4};
    break;
  case 8: // Nearly external edge test
    testname = "Nearly external edge";
    boundary.setDefaultBoundary(0);
    points = {0.99, 0.4, 0.995, 0.5, 0.985, 0., 0.99, 0.6};
    break;
  case 9: // Nearly external edge test
    testname = "Nearly external edge 2";
    boundary.setDefaultBoundary(0);
    points = {0.99, 0.4, 0.995, 0.5, 0.985, 0.5, 0.99, 0.6};
    break;
  case 10: // Two generators
    testname = "Two generators";
    boundary.setDefaultBoundary(0);
    points = {1.0, 0.3, 0.5, -0.4};
    break;
  case 11: // Clipped point
    testname = "Clipped point";
    boundary.setDefaultBoundary(0);
    points = {0.9, 0.3, 0.5, -0.4, 1.01, 0.9}; // The last generator should be removed
    break;
  case 12: // Collinear points
    testname = "Collinear points";
    boundary.setDefaultBoundary(0);
    points = {0.1, 0.3, 0.7, 0.3, 0.4, 0.3};
    break;
  case 13: // Clipped collinear points
    {
      testname = "Clipped collinear points";
      boundary.setDefaultBoundary(5);
      vector<int> starIndices = {4, 5, 6};
      for (auto& pi : starIndices) {
        int i = boundary.mPLC.holes[0][pi][0];
        points.push_back(boundary.mPLCpoints[2*i]);
        points.push_back(-0.7);
      }
      break;
    }
  }
  cout << "\n=== " << outname << " Test " << tnum << ":  " << testname << " ===" << endl;
  auto Q = boundary.mQ;
  QuantTessellation<2> quantMesh(Q, points);
  QuantPLC<2> QPLC(boundary.mPLC, Q, boundary.mPLCpoints);
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

} // anonymous namespace

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif
  const int numtest = 13;
  try {
    bool boost = true;
    for (int i = 0; i < 2; ++i) {
      for (int test = 1; test <= numtest; ++test) {
        tests(test, boost);
      }
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
