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
#include "Point.hh"

#ifdef POLYTOPE_ENABLE_MPI
// extern "C" {
#include "mpi.h"
// }
#endif

#define POLY_CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

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
int main(int argc, char** argv) {

  // Initialize MPI.
  MPI_Init(&argc, &argv);

  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  const double degeneracy = 1.0e-10;
  const double lx = (x2 - x1), ly = (y2 - y1);

  // Figure out our parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Seed the random number generator the same on all processes.
  srand(10489592);

  // Try tessellating increasing numbers of generators.
  for (unsigned nx = 10; nx != 50; ++nx) {
    if (rank == 0) cout << "Testing nx=" << nx << endl;

    // Create the seed positions for each domain.  Note we rely on this sequence
    // being the same for all processors and therefore don't need to communicate
    // this information.
    vector<double> xproc, yproc;
    xproc.reserve(numProcs);
    yproc.reserve(numProcs);
    for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
      xproc.push_back(x1 + 0.5*(iproc + 0.5)*(x2 - x1));
      yproc.push_back(y1 + 0.5*(y2 - y1));
      // xproc.push_back(x1 + random01()*(x2 - x1));
      // yproc.push_back(y1 + random01()*(y2 - y1));
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

    // Create the tessellation.
    double xmin[2] = { x1, y1 };
    double xmax[2] = { x2, y2 };
    polytope::Tessellation<2, double> mesh;
    
    polytope::DistributedTessellator<2, double> distVoro(new polytope::VoroPP_2d<double>(),
                                                         true, true);
    distVoro.tessellate(generators, xmin, xmax, mesh);

    // Figure out which of our nodes and faces we actually own.
    unsigned ncells = mesh.cells.size();
    unsigned nnodes = mesh.nodes.size()/2;
    unsigned nfaces = mesh.faces.size();
    vector<unsigned> ownNodes(nnodes, 1), ownFaces(nfaces, 1);
    for (unsigned k = 0; k != mesh.sharedNodes.size(); ++k) {
      cerr << "Talking to " << mesh.neighborDomains[k] << endl;
      if (mesh.neighborDomains[k] < rank) {
        for (unsigned j = 0; j != mesh.sharedNodes[k].size(); ++k) ownNodes[mesh.sharedNodes[k][j]] = 0;
        for (unsigned j = 0; j != mesh.sharedFaces[k].size(); ++k) ownFaces[mesh.sharedFaces[k][j]] = 0;
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
//     cout << "Node positions: " << endl;
//     for (unsigned i = 0; i != mesh.nodes.size()/2; ++i) {
//       cout << "   Node " << i << " @ (" << mesh.nodes[2*i] << " " << mesh.nodes[2*i + 1] << ")" << endl;
//     }
//     cout << "Face node sets: " << endl;
//     for (unsigned i = 0; i != nx*nx; ++i) {
//       cout << "   FACES for mesh cell " << i << " :";
//       for (unsigned j = 0; j != mesh.cells[i].size(); ++j) cout << " " << mesh.cells[i][j];
//       cout << endl;
//     }
//     for (unsigned i = 0; i != mesh.faces.size(); ++i) {
//       double xf = 0.0, yf = 0.0;
//       cout << "   NODES for mesh face " << i << " :";
//       for (unsigned j = 0; j != mesh.faces[i].size(); ++j) {
//         unsigned k = mesh.faces[i][j];
//         cout << " " << k;
//         xf += mesh.nodes[2*k];
//         yf += mesh.nodes[2*k + 1];
//       }
//       xf /= mesh.faces[i].size();
//       yf /= mesh.faces[i].size();
//       cout << " @ (" << xf << " " << yf << ")"  << endl;
//     }

#ifdef HAVE_SILO
    // Blago!
    {
      vector<double> r2(mesh.cells.size(), rank), rownNodes(nnodes), rownFaces(nfaces);
      for (unsigned i = 0; i != nnodes; ++i) rownNodes[i] = ownNodes[i];
      for (unsigned i = 0; i != nfaces; ++i) rownFaces[i] = ownFaces[i];
      map<string, double*> nodeFields, edgeFields, faceFields, cellFields;
      cellFields["domain"] = &r2[0];
      nodeFields["ownNodes"] = &rownNodes[0];
      faceFields["ownFaces"] = &rownFaces[0];
      ostringstream os;
      os << "test_DistributedTessellator_" << nx << "x" << nx << "_lattice_" << numProcs << "domains";
      polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, os.str());
    }
    // Blago!
#endif

    // Check the global sizes.
    POLY_CHECK(nnodesGlobal == (nx + 1)*(nx + 1));
    POLY_CHECK(ncellsGlobal == nx*nx);
    for (unsigned i = 0; i != ncells; ++i) POLY_CHECK(mesh.cells[i].size() == 4);
    POLY_CHECK(nfacesGlobal == 2*nx*(nx + 1));

    // Check that everyone agrees who talks to who.
    for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
      unsigned numOthers = mesh.neighborDomains.size();
      vector<unsigned> otherNeighbors(mesh.neighborDomains);
      MPI_Bcast(&numOthers, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      if (numOthers > 0) {
        otherNeighbors.resize(numOthers);
        MPI_Bcast(&(otherNeighbors.front()), numOthers, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        POLY_CHECK(rank == sendProc or
               count(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), sendProc) == 
               count(otherNeighbors.begin(), otherNeighbors.end(), rank));
      }
    }

    // Create hashes of our local node and face positions.
    typedef polytope::Point2<uint64_t> Point;
    vector<Point> localNodeHashes, localFaceHashes;
    const double lmax = max(lx, ly);
    for (unsigned i = 0; i != nnodes; ++i) localNodeHashes.push_back(Point((mesh.nodes[2*i]     - x1)/lmax,
                                                                           (mesh.nodes[2*i + 1] - y1)/lmax,
                                                                           degeneracy, i));
    for (unsigned i = 0; i != nfaces; ++i) {
      POLY_CHECK(mesh.faces[i].size() == 2);
      const unsigned i1 = mesh.faces[i][0];
      const unsigned i2 = mesh.faces[i][1];
      const double xf = (0.5*(mesh.nodes[2*i1] + mesh.nodes[2*i2]) - x1)/lmax;
      const double yf = (0.5*(mesh.nodes[2*i1 + 1] + mesh.nodes[2*i2 + 1]) - y1)/lmax;
      localFaceHashes.push_back(Point(xf, yf, degeneracy, i));
    }

    // Each processor sends it's local hashes to each of it's process neighbors that have a greater rank.
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mesh.neighborDomains.size());
    list<vector<char> > localBuffers;
    list<unsigned> localBufSizes;
    for (unsigned i = 0; i != mesh.neighborDomains.size(); ++i) {
      const unsigned otherProc = mesh.neighborDomains[i];
      if (otherProc > rank) {
        vector<Point> sendNodeHashes, sendFaceHashes;
        for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[i].begin();
             itr != mesh.sharedNodes[i].end();
             ++itr) sendNodeHashes.push_back(localNodeHashes[*itr]);
        for (vector<unsigned>::const_iterator itr = mesh.sharedFaces[i].begin();
             itr != mesh.sharedFaces[i].end();
             ++itr) sendFaceHashes.push_back(localFaceHashes[*itr]);
        localBuffers.push_back(vector<char>());
        serialize(sendNodeHashes, localBuffers.back());
        serialize(sendFaceHashes, localBuffers.back());
        localBufSizes.push_back(localBuffers.back().size());
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBufSizes.back(), 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &sendRequests.back());
        if (localBufSizes.back() > 0) {
          sendRequests.push_back(MPI_Request());
          MPI_Isend(&localBuffers.back().front(), localBufSizes.back(), MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
        }
      }
    }

    // Now go over each of our neighbors and look for who we're receiving from.
    for (unsigned i = 0; i != mesh.neighborDomains.size(); ++i) {
      const unsigned otherProc = mesh.neighborDomains[i];
      if (mesh.neighborDomains[i] < rank) {

        // Get the other processors hashed positions.
        unsigned bufSize;
        MPI_Status recvStatus1, recvStatus2;
        MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &recvStatus1);
        if (bufSize > 0) {
          vector<char> buffer(bufSize, '\0');
          MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &recvStatus2);
          vector<Point> otherNodeHashes, otherFaceHashes;
          vector<char>::const_iterator itr = buffer.begin();
          deserialize(otherNodeHashes, itr, buffer.end());
          deserialize(otherFaceHashes, itr, buffer.end());
          POLY_ASSERT(itr == buffer.end());

          // Check that the other processes node and face positions line up with ours.
          const unsigned nn = mesh.sharedNodes[i].size(), nf = mesh.sharedFaces[i].size();
          POLY_CHECK(otherNodeHashes.size() == nn);
          POLY_CHECK(otherFaceHashes.size() == nf);
          for (unsigned j = 0; j != nn; ++j) {
            const unsigned k = mesh.sharedNodes[i][j];
            cerr << "Checking " << localNodeHashes[k] << " " << otherNodeHashes[j] << endl;
            POLY_CHECK(localNodeHashes[k] == otherNodeHashes[j]);
          }
          for (unsigned j = 0; j != nf; ++j) {
            const unsigned k = mesh.sharedFaces[i][j];
            POLY_CHECK(localFaceHashes[k] == otherFaceHashes[j]);
          }
        }
      }
    }

    // Make sure all our sends are completed.
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
  }

  
  MPI_Finalize();
  return 0;
}
