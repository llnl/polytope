// test_twoGenerators
//
// If coordMax in TriangleTessellator is set too low, then a wrap-around effect
// will become visible when trying to tessellate a unit-square boundary

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {

  // output name
  string testName = "TwoGenerators_" + tessellator.name();

  // the domain boundary
  Boundary2D<double> boundary;
  boundary.setUnitSquare();

  // generators
  vector<double> points;
  points.push_back(-0.25);  points.push_back(-0.125);
  points.push_back( 0.25);  points.push_back( 0.125);
  
  // Make that mesh!
  Tessellation<2,double> mesh;
  tessellator.tessellate( points, boundary.mPLCpoints, boundary.mPLC, mesh );
  outputMesh(mesh, testName, points);

  // Some post-conditions
  POLY_CHECK(mesh.nodes.size()/2 == 6);
  POLY_CHECK(mesh.cells.size()   == 2);
  POLY_CHECK(mesh.faces.size()   == 7);
  for (unsigned i = 0; i != 2; ++i) POLY_CHECK(mesh.cells[i].size() == 4);
  const double area = computeTessellationArea( mesh );
  const double fracerr = std::abs(boundary.mArea - area)/boundary.mArea;
  const double tol = 1.0e-6;
  POLY_CHECK2( fracerr < tol, ""
               << "              Area  = " << area << endl
               << "              Error = " << boundary.mArea - area << endl
               << "   Fractional error = " << fracerr << endl );
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    test(tessellator);
  }
#endif

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    test(tessellator);
  }
#endif      


  
  
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
