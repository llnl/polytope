// test_RandomPoints
//
// Stress test for meshing complicated PLC boundaries with/without holes.
// Iterate over each of the default boundaries defined in Boundary2D.hh
// and tessellate using N randomly-distributed generators for N=10,100,1000.
// Can test both Triangle and Boost 2D tessellators.

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"
#include "BoostTessellator.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

#include "Communicator.hh" 

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Boundary2D& boundary,
                  Tessellator<2, double>& tessellator,
                  int boundaryID, int numSweeps) {
  // output name
  ostringstream os;
  os << "RandomPoints_" << tessellator.name();
  string testName = os.str();

  Generators<2> generators(boundary);
  unsigned nPoints = 10;
  Tessellation<2,double> mesh;
  int seed = 1049600;
  for (unsigned n = 0; n < numSweeps; ++n) {
    POLY_CHECK(mesh.empty());
    nPoints = nPoints * 10;
    int plotIndex = 3*boundaryID + n;

    cout << nPoints << " points..." << endl;
    generators.randomPoints( nPoints, seed );
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, plotIndex);
    compareArea(boundary, mesh);
    mesh.clear();
    plotIndex++;
  }
}

// -----------------------------------------------------------------------
// testAllBoundaries
// -----------------------------------------------------------------------
void testAllBoundaries(Tessellator<2, double>& tessellator, int numSweeps) {
  for (int bid = 0; bid < 11; ++bid) {
    cout << "Testing boundary type " << bid << endl;
    Boundary2D boundary;
    boundary.mDiff = 0.5;
    boundary.setDefaultBoundary(bid);
    testBoundary(boundary, tessellator, bid, numSweeps);
  }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  int numSweeps = 2;
  if (argc >= 2) {
    numSweeps = std::stoi(argv[1]);
  }

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator tessellator;
    testAllBoundaries(tessellator, numSweeps);
  }
#endif

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator tessellator;
    testAllBoundaries(tessellator, numSweeps);
  }
#endif
  comm.finalize();
  return 0;
}
