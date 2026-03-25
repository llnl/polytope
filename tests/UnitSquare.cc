// test_UnitSquare
//
// Mesh a unit square with (nx-by-nx) Cartesian generators for nx in [2,100].
// Perform checks on the resulting tessellation to see if it is indeed Cartesian.
// Triangle and Boost tessellators are tested here.
// -----------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"

#include "timingUtilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// checkCartesianMesh
// -----------------------------------------------------------------------
void checkCartesianMesh(Tessellation<2,double>& mesh, unsigned nx, unsigned ny)
{
   POLY_CHECK(mesh.nodes.size()/2 == (nx + 1)*(ny + 1) );
   POLY_CHECK(mesh.cells.size()   == nx*ny );
   POLY_CHECK(mesh.faces.size()   == nx*(ny + 1) + ny*(nx + 1) );
   for (unsigned i = 0; i != nx*ny; ++i) POLY_CHECK(mesh.cells[i].size() == 4);
   
   std::vector<std::set<unsigned> > nodeCells = mesh.computeNodeCells();
   for (unsigned i = 0; i != (nx+1)*(ny+1); ++i)
   {
      POLY_CHECK( (nodeCells[i].size() == 4) ||
                  (nodeCells[i].size() == 2) ||
                  (nodeCells[i].size() == 1) );
   }
}

// -----------------------------------------------------------------------
// generateMesh
// -----------------------------------------------------------------------
void generateMesh(Tessellator<2,double>& tessellator)
{
  // Set the boundary
  Boundary2D<double> boundary;
  boundary.setUnitSquare();
  Generators<2,double> generators( boundary );
  
  const unsigned Nmin   = 2;
  const unsigned Nmax   = 101;

  
  for (unsigned nx = Nmin; nx != Nmax; ++nx) {
    cout << "Testing nx=" << nx << endl;

    // Create generators
    std::vector<unsigned> nxny(2,nx);
    generators.cartesianPoints(nxny);
    Tessellation<2,double> mesh;
     
    Timing::Time start = Timing::currentTime();
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    double time = Timing::difference(start, Timing::currentTime());
    
    cout << "Mesh construction took " << time << " seconds." << endl;
    
    // CHECKS:
    cout << "   num mesh nodes : " << mesh.nodes.size()/2 << endl;
    cout << "   num mesh cells : " << mesh.cells.size()   << endl;
    cout << "   num mesh faces : " << mesh.faces.size()   << endl;
    checkCartesianMesh(mesh,nx,nx);

    // ofstream os;
    // os.open("oneMillionRandomGenerators.txt");
    // os << "# x \t\t y \t\t z" << endl;
    // for (unsigned j = 0; j != generators.mPoints.size()/2; ++j) {
    //   os << generators.mPoints[2*j] << " \t " 
    //      << generators.mPoints[2*j+1] << " \t " << 0.0 << endl;
    // }
    // os.close();
  }
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
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
    generateMesh(tessellator);
  }
#endif   

#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    generateMesh(tessellator);
  }
#endif

   // cout << "\nVoro 2D Tessellator:\n" << endl;
   // VoroPP_2d<double> voro;
   // generateMesh(voro);
   // cout << "Voro 2D: PASS" << endl;

  cout << "PASS" << endl;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
