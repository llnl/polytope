// test_DistributedVisible2D
//
// Can test both Triangle and Boost 2D tessellators.
// This tests the generator exchange operations in the
// DistributedTessellator. Requires 6 ranks to run properly.

#include <cmath>
#include <exception>
#include <iostream>
#include <vector>

#include "polytope.hh"

#include "Communicator.hh"

#include "BoostTessellator.hh"
#include "DistributedTessellator.hh"
#include "PLC.hh"
#include "Tessellation.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace polytope;

namespace {

void test(int testNum, Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  int nranks = Communicator::getNRanks();

  const int Ngen = 16;
  const int activeRanks = 6;
  int seed = 1049609;
  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(0);
  Generators<2> generators(boundary);
  generators.randomPoints(Ngen, seed);
  // Distribute points to specific ranks
  const int usedRanks = std::min(nranks, activeRanks);
  const int nperRank = Ngen/usedRanks;
  std::vector<unsigned> procIndex(usedRanks);
  for (int i = 0; i < usedRanks; ++i) {
    procIndex[i] = i*nperRank;
  }
  auto finalRanks = generators.distributePointsAmongRanks(procIndex);

  // Reference values for a 6 processor test run
  std::vector<unsigned> finalProcN(activeRanks, 0);
  for (auto& fr : finalRanks) {
    finalProcN[fr.second]++;
  }
  // List of expected neighbor ranks for each rank
  std::vector<std::vector<int>> refNeighbors(activeRanks);
  refNeighbors[0] = {4, 1, 5};
  refNeighbors[1] = {4, 0, 5};
  refNeighbors[2] = {4};
  refNeighbors[3] = {5};
  refNeighbors[4] = {2, 0, 1};
  refNeighbors[5] = {0, 1, 3};
  // Expected number of generator points after exchanging
  std::vector<unsigned> refNPoints(activeRanks, 0);
  for (int p = 0; p < activeRanks; ++p) {
    refNPoints[p] = finalProcN[p];
    for (auto& np : refNeighbors[p]) {
      refNPoints[p] += finalProcN[np];
    }
  }

  Communicator::Barrier();
  DistributedTessellator<2> distributed(tessellator);
  Tessellation<2, double> localMesh;
  distributed.tessellate(generators.mPoints, localMesh);

  const auto localCells = static_cast<int>(localMesh.cells.size());
  int totalCells = 0;
  MPI_Allreduce(&localCells, &totalCells, 1, MPI_INT, MPI_SUM, Communicator::communicator());

  std::string outname = "DistributedVisible_" + tessellator.name();
  outputMesh(localMesh, outname, testNum, 0.);

  // Now get the current mesh including all it's neighbor generators
  QuantTessellation<2> qmesh(generators.mPoints);
  distributed.tessellateQuantized(qmesh);
  if (rank < activeRanks && nranks == 6) {
    Tessellation<2, double> procMesh;
    qmesh.fillTessellation(procMesh);
    SiloWriter<2, Tessellation<2, double>>::write(procMesh, "ProcMesh" + tessellator.name(), 1, rank);
    POLY_CHECK2(int(qmesh.points.size()) == refNPoints[rank],
                "Number of generator points for rank " << rank << "incorrect, "
                << "expected " << refNPoints[rank] << " actual " << qmesh.points.size());
  }

  Communicator::Barrier();
  if (rank == root) {
    cout << "\nTest passed." << endl;
  }
}
} // anonymous namespace

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);
  const int root = Communicator::getRoot();

#ifdef POLYTOPE_ENABLE_TRIANGLE
   {
     if (Communicator::getRank() == root) {
       cout << "\nTriangle Tessellator:\n" << endl;
     }
     TriangleTessellator tessellator;
     test(0, tessellator);
   }
#endif

   {
     if (Communicator::getRank() == root) {
       cout << "\nBoost Tessellator:\n" << endl;
     }
     BoostTessellator tessellator;
     test(0, tessellator);
   }
   if (Communicator::getRank() == root) {
     std::cout << "=== DistributedVoronoi2D passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
