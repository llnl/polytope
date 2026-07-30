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
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif
#include "Communicator.hh"

using namespace polytope;
using namespace std;

namespace {

using CoordType = typename HashKey<2>::IntType;;
using IntPoint = Point2<CoordType>;

void tests(const int tnum, bool boostTess) {
  std::string outname = (boostTess) ? "boost" : "triangle";
  Boundary2D boundary;
  boundary.mDiff = 1.;
  std::vector<double> points;
  std::vector<IntPoint> qpoints;
  std::string testname;
  // Number of nodes expected from test if not -1
  int numNodes = -1;
  int circleNodes = 90; // Number of nodes used to make circles
  // Note: quantizer is initialized in the setDefaultBoundary calls
  switch (tnum) {
  case 1: // Square
    testname = "Square";
    boundary.setDefaultBoundary(0);
    points = {-0.05, -0.05, 0.05, -0.05, 0.05, 0.05, -0.05, 0.05};
    numNodes = 9;
    break;
  case 2: // Cut square
    {
      testname = "Cut square";
      boundary.setDefaultBoundary(0);
      points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
      numNodes = 11;
      break;
    }
  case 3: // Cut square with hole
    {
      testname = "Cut square with hole";
      boundary.setDefaultBoundary(11);
      points = {0.05, 0.025, 0.025, 0.05, 0.05, -0.05, -0.05, -0.05, -0.05, 0.05};
      numNodes = 17;
    }
    break;
  case 4: // Diamond
    testname = "Diamond";
    boundary.setDefaultBoundary(0);
    points = {-0.4, 0, 0.4, 0, 0, -0.4, 0, 0.4};
    numNodes = 5;
    break;
  case 5: // Obtuse
    testname = "Obtuse triangle";
    boundary.setDefaultBoundary(0);
    points = {0.67, -0.14, 0.91, 0.3, 0.49, -0.4};
    numNodes = 8;
    break;
  case 6: // Obtuse with star
    boundary.mDiff = 0.5;
    testname = "Obtuse triangle clipped by star";
    boundary.setDefaultBoundary(5);
    points = {0.67, -0.14, 0.91, 0.3, 0.49, -0.4};
    numNodes = circleNodes + 16;
    break;
  case 7: // Mod pad obtuse with star
    testname = "Mod pad obtuse triangle clipped by star";
    boundary.m_pad = 0.5;
    boundary.setDefaultBoundary(10);
    points = {0.67, -0.14, 0.91, 0.3, 0.49, -0.4};
    numNodes = 20;
    break;
  case 8: // Nearly external edge test
    testname = "Nearly external edge";
    boundary.setDefaultBoundary(0);
    points = {0.99, 0.4, 0.995, 0.5, 0.985, 0., 0.99, 0.6};
    numNodes = 10;
    break;
  case 9: // Nearly external edge test
    testname = "Nearly external edge 2";
    boundary.setDefaultBoundary(0);
    points = {0.99, 0.4, 0.995, 0.5, 0.985, 0.5, 0.99, 0.6};
    numNodes = 10;
    break;
  case 10: // Two generators
    testname = "Two generators";
    boundary.setDefaultBoundary(0);
    points = {1.0, 0.3, 0.5, -0.4};
    numNodes = 6;
    break;
  case 11: // Clipped point
    testname = "Clipped point";
    boundary.setDefaultBoundary(0);
    points = {0.9, 0.3, 0.5, -0.4, 1.01, 0.9}; // The last generator should be removed
    numNodes = 6;
    break;
  case 12: // Collinear points
    testname = "Collinear points";
    boundary.setDefaultBoundary(0);
    points = {0.1, 0.3, 0.7, 0.3, 0.4, 0.3};
    numNodes = 8;
    break;
  case 13: // Clipped collinear points
    {
      testname = "Clipped collinear points";
      boundary.setDefaultBoundary(10);
      vector<int> starIndices = {4, 5, 6};
      for (auto& pi : starIndices) {
        int i = boundary.mPLC.holes[0][pi][0];
        points.push_back(boundary.mPLCpoints[2*i]);
        points.push_back(-0.7);
      }
      break;
    }
  case 14: // Difficult set of quantized points discovered during rotation test
    {
      testname = "Difficult quantized points";
      boundary.setDefaultBoundary(0);
      std::vector<int> rqp = {547445970, 380834812,
                              545899097, 374153061,
                              544323133, 368524286,
                              542713618, 709756041,
                              547445970, 692907011,
                              630394586, 932588752};
      qpoints = extractCoords<2, CoordType>(rqp);
      break;
    }
  case 15: // Difficult set of quantized points discovered during rotation test
    {
      testname = "Difficult quantized points 2";
      boundary.setDefaultBoundary(0);
      std::vector<int> rqp = {958459259, 368147224,
                              947586761, 544198943,
                              772074713, 824295508,
                              950857238, 541497367,
                              960886415, 365110638};
      qpoints = extractCoords<2, CoordType>(rqp);
      break;
    }
  }
  cout << "\n=== " << outname << " Test " << tnum << ":  " << testname << " ===" << endl;
  QuantTessellation<2> quantMesh;
  if (qpoints.size() > 0) {
    quantMesh.init(qpoints);
  } else {
    quantMesh.init(points);
  }
  QuantPLC<2> QPLC(boundary.mPLC, boundary.mPLCpoints);
  quantMesh.cullExternalPoints(QPLC);
  if (boostTess) {
    BoostTessellator boost;
    boost.tessellateQuantized(quantMesh);
    quantMesh.clipTessellation(QPLC, boost);
#ifdef POLYTOPE_ENABLE_TRIANGLE
  } else {
    TriangleTessellator tri;
    tri.tessellateQuantized(quantMesh);
    quantMesh.clipTessellation(QPLC, tri);
#endif
  }
  Tessellation<2, double> mesh;
  quantMesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  outputMesh(mesh, outname, tnum, double(tnum));
  compareArea(boundary, mesh);
  testWatertight(mesh, boundary.mPLC.holes.size());
  if (numNodes > 0) {
    // Ideally we would match nodes exactly but determine collinearity exactly is not really
    // possible with quantized coordinates
    POLY_CHECK2(mesh.nodes.size() >= numNodes, "We must have at least " << numNodes
                << " but we only have " << mesh.nodes.size()/2);
  }
}

} // anonymous namespace

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  const int numtest = 15;
  try {
#ifdef POLYTOPE_ENABLE_TRIANGLE
    for (bool boost : {true, false}) {
      for (int test = 1; test <= numtest; ++test) {
        tests(test, boost);
      }
    }
#else
    for (int test = 1; test <= numtest; ++test) {
      tests(test, true);
    }
#endif

    cout << "\n=== ALL TESTS PASSED ===" << endl;
  } catch (const exception& e) {
    cout << "\n=== TEST FAILED WITH EXCEPTION ===" << endl;
    cout << e.what() << endl;
    comm.finalize();
    return 1;
  }

  comm.finalize();
  return 0;
}
