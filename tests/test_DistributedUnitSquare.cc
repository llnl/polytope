// Try tessellating a simple lattice of generators in a box in parallel.
// We use randomly chosen seed locations to divide up the generators
// between processors.

#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <stdlib.h>
#include <limits>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "checkDistributedTessellation.hh"
#include "Point.hh"

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
// The test itself
//------------------------------------------------------------------------------
void runTest(Tessellator<2,double>& tessellator) {

  // Parameters
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  const unsigned Nmin = 50;
  const unsigned Nmax = 100;

  // Figure out our parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (rank == 0) cout << "\nTesting Distributed " << tessellator.name() << endl;

  // Seed the random number generator the same on all processes.
  srand(10489591);

  // Try tessellating increasing numbers of generators.
  for (unsigned nx = Nmin; nx <= Nmax; ++nx) {
    if (rank == 0) cout << "Testing nx=" << nx << endl;

    // Create the seed positions for each domain.  Note we rely on this sequence
    // being the same for all processors and therefore don't need to communicate
    // this information.
    vector<double> xproc, yproc;
    xproc.reserve(numProcs);
    yproc.reserve(numProcs);
    for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
      xproc.push_back(x1 + random01()*(x2 - x1));
      yproc.push_back(y1 + random01()*(y2 - y1));
    }

    // Create the local generators.  Note this is not efficient in a couple of ways!  
    // All processes are walking all generators and checking which ones belong to them,
    // and the processor search process is N^2 in the number of processors.  But crimine,
    // this is just supposed to be a little unit test!
    vector<double> generators;
    const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
    unsigned ix, iy;
    double xi, yi;
    for (iy = 0; iy != nx; ++iy) {
      yi = std::max(y1, std::min(y2, y1 + (iy + 0.5)*dy));
      for (ix = 0; ix != nx; ++ix) {
        xi = std::max(x1, std::min(x2, x1 + (ix + 0.5)*dx));
        unsigned owner = 0;
        double minDist2 = distance2(xi, yi, xproc[0], yproc[0]);
        for (unsigned iproc = 1; iproc < numProcs; ++iproc) {
          const double d2 = distance2(xi, yi, xproc[iproc], yproc[iproc]);
          if (d2 < minDist2) {
            owner = iproc;
            minDist2 = d2;
          }
        }
        if (rank == owner) {
          generators.push_back(xi);
          generators.push_back(yi);
        }
      }
    }
    
    
    unsigned numGen = generators.size()/2;
    POLY_CHECK2( numGen > 1, "Processor " << rank << " only has " << numGen << " generators. This is not enough to tessellate!");

    // Create the tessellation.
    double xmin[2] = { x1, y1 };
    double xmax[2] = { x2, y2 };
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators, xmin, xmax, mesh);

    // vector<double> PLCpoints;
    // PLCpoints.push_back(xmin[0]);  PLCpoints.push_back(xmin[1]);
    // PLCpoints.push_back(xmax[0]);  PLCpoints.push_back(xmin[1]);
    // PLCpoints.push_back(xmax[0]);  PLCpoints.push_back(xmax[1]);
    // PLCpoints.push_back(xmin[0]);  PLCpoints.push_back(xmax[1]);
    // polytope::PLC<2> boundary;
    // boundary.facets.resize(4, vector<int>(2));
    // for (int i = 0; i != 4; ++i) { boundary.facets[i][0] = i;  boundary.facets[i][1] = (i+1)%4; }
    // distTest.tessellate(generators, PLCpoints, boundary, mesh);


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
    // cout << "Node positions: " << endl;
    // for (unsigned i = 0; i != mesh.nodes.size()/2; ++i) {
    //    cout << "   " << rank << "Node " << i << " @ (" << mesh.nodes[2*i] << " " << mesh.nodes[2*i + 1] << ")" << endl;
    // }
    // cout << "Face node sets: " << endl;
    // for (unsigned i = 0; i != nx*nx; ++i) {
    //   cout << "   FACES for mesh cell " << i << " :";
    //   for (unsigned j = 0; j != mesh.cells[i].size(); ++j) cout << " " << mesh.cells[i][j];
    //   cout << endl;
    // }
    // for (unsigned i = 0; i != mesh.faces.size(); ++i) {
    //   double xf = 0.0, yf = 0.0;
    //   cout << "   NODES for mesh face " << i << " :";
    //   for (unsigned j = 0; j != mesh.faces[i].size(); ++j) {
    //     unsigned k = mesh.faces[i][j];
    //     cout << " " << k;
    //     xf += mesh.nodes[2*k];
    //     yf += mesh.nodes[2*k + 1];
    //   }
    //   xf /= mesh.faces[i].size();
    //   yf /= mesh.faces[i].size();
    //   cout << " @ (" << xf << " " << yf << ")"  << endl;
    // }
    

    // for( unsigned i = 0; i != mesh.cells.size(); ++i ){
    //    cout << endl << rank << "Cell " << i << ": ";
    //    for( std::vector<int>::const_iterator faceItr = mesh.cells[i].begin();
    //         faceItr != mesh.cells[i].end(); ++faceItr){
    //       const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
    //       POLY_CHECK( mesh.faces[iface].size() == 2 );
    //       const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] : mesh.faces[iface][0];
    //       cout << " (" << mesh.nodes[2*inode] << "," << mesh.nodes[2*inode+1] << ") ";
    //    }
    // }

    
    // for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc ){
    //    unsigned otherProc = mesh.neighborDomains[iproc];
    //    for (std::vector<unsigned>::const_iterator nodeItr = mesh.sharedNodes[iproc].begin();
    //         nodeItr != mesh.sharedNodes[iproc].end(); ++nodeItr){
    //       cout << rank << "<-->" << otherProc << ":    ( " 
    //            << mesh.nodes[2*(*nodeItr)  ] << " , "
    //            << mesh.nodes[2*(*nodeItr)+1] << " ) " << endl;
    //    }
    // }
    
    // for (unsigned iproc = 0; iproc < numProcs; ++iproc){
    //    if( rank == iproc ){
    //       cout << endl << "DOMAIN " << iproc << endl;
    //       for (unsigned i = 0; i < mesh.nodes.size()/2; ++i ){
    //          cout << i << ":   ( " << mesh.nodes[2*i  ] 
    //               << " , " << mesh.nodes[2*i+1] << " ) " << endl;
    //       }
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    // }
    
    
    // Blago!
#ifdef HAVE_SILO
    {
      if (nx % 5 == 0) {
      vector<double> r2(mesh.cells.size(), rank), rownNodes(nnodes), rownFaces(nfaces);
      for (unsigned i = 0; i != nnodes; ++i) rownNodes[i] = ownNodes[i];
      for (unsigned i = 0; i != nfaces; ++i) rownFaces[i] = ownFaces[i];
      map<string, double*> nodeFields, edgeFields, faceFields, cellFields;
      cellFields["domain"] = &r2[0];
      nodeFields["ownNodes"] = &rownNodes[0];
      faceFields["ownFaces"] = &rownFaces[0];
      ostringstream os;
      os << "DistributedUnitSquare_" << tessellator.name() << "_" << numProcs << "domains";
      polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                             faceFields, cellFields, os.str(),
                                             nx, 0.0);
      }
    }
#endif
    // Blago!

    // Check the global sizes.
    POLY_CHECK2(nnodesGlobal == (nx + 1)*(nx + 1), nnodesGlobal << " != " << (nx + 1)*(nx + 1));
    POLY_CHECK2(ncellsGlobal == nx*nx, ncellsGlobal << " != " << nx*nx);
    for (unsigned i = 0; i != ncells; ++i) POLY_CHECK2(mesh.cells[i].size() == 4, mesh.cells[i].size() << " != " << 4);
    POLY_CHECK2(nfacesGlobal == 2*nx*(nx + 1), nfacesGlobal << " != " << 2*nx*(nx + 1));

    // We can delegate checking the correctness of the parallel data structures to a helper method.
    const string parCheck = checkDistributedTessellation(mesh);
    POLY_CHECK2(parCheck == "ok", parCheck);


    MPI_Barrier(MPI_COMM_WORLD);
  }
}


//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Initialize MPI.
  MPI_Init(&argc, &argv);

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    DistributedTessellator<2, double> tessellator
       (new TriangleTessellator<double>(), true, true);
    runTest(tessellator);
  }
#endif

#ifdef POLYTOPE_ENABLE_BOOST
  {
    DistributedTessellator<2, double> tessellator
       (new BoostTessellator<double>(), true, true);
    runTest(tessellator);
  }
#endif

  
  MPI_Finalize();
  return 0;
}
