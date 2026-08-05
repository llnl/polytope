// test_DistributedVisible2D
//
// Can test both Triangle and Boost 2D tessellators.
// This tests the generator exchange operations in the
// DistributedTessellator

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

void test(Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();

  const int Ngen = 100;
  int seed = 19843;
  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(0);
  Generators<2> generators(boundary);
  generators.randomPoints(Ngen, seed);
  // Distribute points to specific ranks
  auto procIndex = generators.assignRandomPointToRank(seed);
  auto finalRanks = generators.distributePointsAmongRanks(procIndex);

  Communicator::Barrier();
  DistributedTessellator<2> distributed(tessellator);
  Tessellation<2, double> localMesh;
  distributed.tessellate(generators.mPoints, boundary.mPLCpoints,
                         boundary.mPLC, localMesh);

  const auto localCells = static_cast<int>(localMesh.cells.size());
  int totalCells = 0;
  MPI_Allreduce(&localCells, &totalCells, 1, MPI_INT, MPI_SUM, Communicator::communicator());

  double localArea = computeTessellationArea(localMesh);
  double distributedArea = 0.0;
  MPI_Allreduce(&localArea, &distributedArea, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());

  std::string outname = "DistributedVisible_" + tessellator.name();
  outputMesh(localMesh, outname, 0, 0.);

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
     test(tessellator);
   }
#endif

   {
     if (Communicator::getRank() == root) {
       cout << "\nBoost Tessellator:\n" << endl;
     }
     BoostTessellator tessellator;
     test(tessellator);
   }
   if (Communicator::getRank() == root) {
     std::cout << "=== DistributedVoronoi2D passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
