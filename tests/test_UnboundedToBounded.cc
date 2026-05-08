// test_TriangleUnboundedToBounded

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// runTest
// -----------------------------------------------------------------------
void runTest(Tessellator<2,double>& tessellator) {
  
  // Parameters
  const int bType = 3;  // M-shape with holes
  const bool addZeroAreaHole = false;
  //srand(58494385);
    
  ostringstream os;
  os << "UnboundedToBounded_" << tessellator.name();
  const string testName = os.str();
  
  // Boundary data
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(bType);
  
  if (addZeroAreaHole && bType==0) {
    boundary.mPLCpoints.push_back(-0.4);  boundary.mPLCpoints.push_back( 0.0);
    boundary.mPLCpoints.push_back(-0.1);  boundary.mPLCpoints.push_back( 0.1);
    boundary.mPLCpoints.push_back( 0.2);  boundary.mPLCpoints.push_back(-0.2);
    boundary.mPLCpoints.push_back(-0.1);  boundary.mPLCpoints.push_back( 0.1);
    
    boundary.mPLC.holes = vector<vector<vector<int> > >(1);
    boundary.mPLC.holes[0].resize(4);
    for (unsigned i = 0; i != 4; ++i) {
      boundary.mPLC.holes[0][i].resize(2);
      boundary.mPLC.holes[0][i][0] = 4 + i;
      boundary.mPLC.holes[0][i][1] = 4 + (i+1)%4;
    }

    boundary.boostMyBoundary();
  }
  
  // Generator data
  Generators<2,double> generators(boundary);
  generators.randomPoints(10);

  Tessellation<2,double> mesh;
  
  // Do the unbounded tessellation
  tessellator.tessellate(generators.mPoints, mesh);
  outputMesh(mesh, testName, generators.mPoints, 1);
  
  // Now do the bounded version
  mesh.clear();
  tessellator.tessellate(generators.mPoints, 
			 boundary.mPLCpoints, 
			 boundary.mPLC, 
			 mesh);
  outputMesh(mesh, testName, generators.mPoints, 2);
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
    runTest(tessellator);
  }
#endif
   
#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    runTest(tessellator);
  }
#endif      
  
  

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
   return 0;
}
