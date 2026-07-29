// test_RandomPoints
//
// Stress test for meshing complicated PLC boundaries with/without holes.
// Iterate over each of the default boundaries defined in Boundary2D.hh
// and tessellate using N randomly-distributed generators for N=10,100,1000.
// Can test both Triangle and Voro++ 2D tessellators. Voro++ has been
// commented out since it currently lacks PLC capabilities.

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"

#include "TetgenTessellator.hh"

#include "Communicator.hh" 

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Tessellator<3, double>& tessellator,
                  int numSweeps) {
  // output name
  ostringstream os;
  os << "TetgenTessellator_" << tessellator.name();
  string testName = os.str();

  Point3<double> pmin(0., 0., 0.);
  Point3<double> pmax(1., 1., 1.);
  // Quantizer is initialized in Generators constructor
  Generators<3> generators(pmin, pmax);
  unsigned nPoints = 10;
  Tessellation<3,double> mesh;
  for (unsigned n = 0; n < numSweeps; ++n) {
    POLY_CHECK(mesh.empty());
    cout << nPoints << " points..." << endl;
    generators.randomPoints(nPoints);
    tessellator.tessellate(generators.mPoints, mesh);
    outputMesh(mesh, testName, n);
    nPoints = nPoints * 10;
    mesh.clear();
  }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

  int numSweeps = 2;
  {
    cout << "\nTetgen Tessellator:\n" << endl;
    TetgenTessellator tessellator;
    testBoundary(tessellator, numSweeps);
  }
  comm.finalize();
  return 0;
}
