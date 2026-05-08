// test_DistributedRandomPoints
//
// For a given default boundary in Boundary2D.hh, we generate N randomly-
// distributed points and arbitrarily assign them to processors.

#include <numeric>
#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "checkDistributedTessellation.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


//------------------------------------------------------------------------------
// Compute the square of the distance.
//------------------------------------------------------------------------------
double distance2(const double x1, const double y1,
                 const double x2, const double y2) {
  return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
}


//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
void runTest(const int bType,
             Tessellator<2,double>& tessellator) {

  // Seed the random number generator the same on all processes.
  srand(10483991);
  
  bool assignRandomly = false;
  
  // Figure out our parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Initialize the bounding PLC
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(bType);

  if (rank == 0) cerr << "Meshing boundary " << bType << endl;
  
  // Generate random points on all processors
  unsigned N = 2000;
  Generators<2,double> generators( boundary );
  generators.randomPoints( N );
  
  if (rank == 0) cerr << "Points generated" << endl;
  
  vector<double> myGenerators;
  
  // Assign random points to each processor
  if (assignRandomly){
    for (unsigned i = 0; i < N; ++i) {
      if (i%numProcs == rank) {
        myGenerators.push_back(generators.mPoints[2*i]  );
        myGenerators.push_back(generators.mPoints[2*i+1]);
      }
    }
  }
  // Assign points to processors in quasi-Voronoi fashion
  else{
    vector<double> xproc, yproc;
    double p[2];
    xproc.reserve(numProcs);
    yproc.reserve(numProcs);
    for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
      boundary.getPointInside(p);
      xproc.push_back(p[0]);
      yproc.push_back(p[1]);
    }
    for (unsigned i = 0; i < N; ++i){
      unsigned owner = 0;
      double minDist2 = distance2(generators.mPoints[2*i], 
                                  generators.mPoints[2*i+1], 
                                  xproc[0], yproc[0]);
      for (unsigned iproc = 1; iproc < numProcs; ++iproc) {
        const double d2 = distance2(generators.mPoints[2*i], 
                                    generators.mPoints[2*i+1], 
                                    xproc[iproc], yproc[iproc]);
        if( d2 < minDist2 ){
          owner = iproc;
          minDist2 = d2;
        }
      }
      if (rank == owner) {
        myGenerators.push_back(generators.mPoints[2*i  ]);
        myGenerators.push_back(generators.mPoints[2*i+1]);
      }
    }
  }
  
  if( rank == 0 ) cerr << "Points assigned to processors" << endl;

  POLY_CHECK(!myGenerators.empty());

  Tessellation<2, double> mesh;
  tessellator.tessellate(myGenerators, boundary.mPLCpoints, boundary.mPLC, mesh);

  if( rank == 0 ) cerr << "Generated tessellation" << endl;

  // Do some sanity checks on the stuff in the shared info.
  unsigned numNeighborDomains = mesh.neighborDomains.size();
  unsigned ncells = mesh.cells.size();
  unsigned nnodes = mesh.nodes.size()/2;
  unsigned nfaces = mesh.faces.size();
  POLY_CHECK(mesh.sharedNodes.size() == numNeighborDomains);
  POLY_CHECK(mesh.sharedFaces.size() == numNeighborDomains);
  POLY_CHECK(mesh.neighborDomains.size() == 0 or *max_element(mesh.neighborDomains.begin(), mesh.neighborDomains.end()) < numProcs);
  for (unsigned k = 0; k != numNeighborDomains; ++k) {
     POLY_CHECK(mesh.sharedNodes[k].size() > 0);
     POLY_CHECK(*max_element(mesh.sharedNodes[k].begin(), mesh.sharedNodes[k].end()) < nnodes);
     POLY_CHECK(mesh.sharedFaces[k].size() == 0 or *max_element(mesh.sharedFaces[k].begin(), mesh.sharedFaces[k].end()) < nfaces);
  }

  if( rank == 0 ) cerr << "Checks passed" << endl;

  // Figure out which of our nodes and faces we actually own.
  vector<unsigned> ownNodes(nnodes, 1), ownFaces(nfaces, 1);
  for (unsigned k = 0; k != mesh.sharedNodes.size(); ++k) {
     if (mesh.neighborDomains[k] < rank) {
        for (unsigned j = 0; j != mesh.sharedNodes[k].size(); ++j) {
           POLY_CHECK(mesh.sharedNodes[k][j] < ownNodes.size());
           ownNodes[mesh.sharedNodes[k][j]] = 0;
        }
        for (unsigned j = 0; j != mesh.sharedFaces[k].size(); ++j) {
           POLY_CHECK(mesh.sharedFaces[k][j] < ownFaces.size());
           ownFaces[mesh.sharedFaces[k][j]] = 0;
        }
     }
  }
  unsigned nnodesOwned = (nnodes == 0U ?
                          0U :
                          accumulate(ownNodes.begin(), ownNodes.end(), 0U));
  unsigned nfacesOwned = (nfaces == 0U ?
                          0U :
                          accumulate(ownFaces.begin(), ownFaces.end(), 0U));

  // Gather some global statistics.
  unsigned ncellsGlobal, nnodesGlobal, nfacesGlobal;
  MPI_Allreduce(&ncells, &ncellsGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nnodesOwned, &nnodesGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nfacesOwned, &nfacesGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

  // Spew the mesh statistics.
  if (rank == 0) {
     cout << "   num mesh cells : " << ncells << " " << ncellsGlobal << endl;
     cout << "   num mesh nodes : " << nnodes << " " << nnodesGlobal << endl;
     cout << "   num mesh faces : " << nfaces << " " << nfacesGlobal << endl;
  }
   
  std::string testName = "DistributedRandomPoints_" + tessellator.name();
  outputMesh(mesh, testName, myGenerators, bType);


  // A helper method checks the correctness of the parallel data structures
  const string parCheck = checkDistributedTessellation(mesh);
  POLY_CHECK2(parCheck == "ok", parCheck);
}


//------------------------------------------------------------------------------
// test all the boundaries
//------------------------------------------------------------------------------
void testAllBoundaries(Tessellator<2,double>& tessellator) {
   for (int btype = 0; btype != 10; ++btype) runTest(btype, tessellator);
}


//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
  // Initialize MPI.
  MPI_Init(&argc, &argv);
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    if (rank == 0) cout << "\nTriangle Tessellator:\n" << endl;
    SerialDistributedTessellator<2, double> tessellator
       (new TriangleTessellator<double>(), true, true);
    testAllBoundaries(tessellator);
  }
#endif


#ifdef POLYTOPE_ENABLE_BOOST
  {
    if (rank == 0) cout << "\nBoost Tessellator:\n" << endl;
    DistributedTessellator<2, double> tessellator
       (new BoostTessellator<double>(), true, true);
    testAllBoundaries(tessellator);
  }
#endif


   
   MPI_Finalize();
   return 0;
}
