// test_DistributedOrderedPoints
//
// Stress test for meshing complicated PLC boundaries with/without holes.
// Iterate over each of the default boundaries defined in Boundary2D.hh
// Can test both Triangle and Boost 2D tessellators.
// This test mimics test_DistributedRandomPoints except points for each
// rank are clumped together.

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

void test(const int btype, Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int nprocs = Communicator::getNProcs();
  int root = Communicator::getRoot();

  const int Ngen = 4*nprocs;
  int seed = 1049600;
  Boundary2D boundary;
  boundary.setDefaultBoundary(btype);
  Generators<2> generators(boundary);
  generators.randomPoints(Ngen, seed);
  std::vector<double> allPoints = generators.mPoints;
  // Distribute points to specific ranks
  auto procIndex = generators.assignRandomPointToRank(seed);
  auto finalRanks = generators.distributePointsAmongRanks(procIndex);

  double serialArea = 0.0;
  if (rank == root) {
    Tessellation<2, double> serialMesh;    
    tessellator.tessellate(allPoints, boundary.mPLCpoints, boundary.mPLC, serialMesh);
    serialArea = computeTessellationArea(serialMesh);
    std::string serialoutname = "SerialOrdered_" + tessellator.name();
    std::vector<double> finalRanksD;
    for (auto& p : serialMesh.points) {
      auto rankItr = finalRanks.find(p);
      POLY_CHECK2(rankItr != finalRanks.end(),
                  "Could not find final rank for serial generator " << p);
      finalRanksD.push_back(static_cast<double>(rankItr->second));
    }
    outputMesh(serialMesh, serialoutname, finalRanksD, "rank", btype, static_cast<double>(btype), 1);
  }

  Communicator::Barrier();
  DistributedTessellator<2> distributed(tessellator);
  Tessellation<2, double> localMesh;
  distributed.tessellate(generators.mPoints, boundary.mPLCpoints,
                         boundary.mPLC, localMesh);
  std::string outname = "DistributedOrdered_" + tessellator.name();
  outputMesh(localMesh, outname, btype, static_cast<double>(btype));

  const auto localCells = static_cast<int>(localMesh.cells.size());
  int totalCells = 0;
  MPI_Allreduce(&localCells, &totalCells, 1, MPI_INT, MPI_SUM, Communicator::communicator());

  POLY_CHECK2(totalCells == static_cast<int>(allPoints.size()/2),
              "Total number of cells " << totalCells
              << " is incorrect; expected " << allPoints.size()/2);

  compareArea(boundary, serialArea, "Serial area failure");
  compareArea(boundary, localMesh, "Distributed area failure");
  if (rank == root) {
    cout << "\nTest " << btype << " passed." << endl;
  }
}
} // anonymous namespace

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);
  const int root = Communicator::getRoot();
  const int Nbstart = 0;
  const int Nbtype = 11;

#ifdef POLYTOPE_ENABLE_TRIANGLE
   {
     if (Communicator::getRank() == root) {
       cout << "\nTriangle Tessellator:\n" << endl;
     }
     TriangleTessellator tessellator;
     for (int btype = Nbstart; btype < Nbtype; ++btype) {
       test(btype, tessellator);
     }
   }
#endif

   {
     if (Communicator::getRank() == root) {
       cout << "\nBoost Tessellator:\n" << endl;
     }
     BoostTessellator tessellator;
     for (int btype = Nbstart; btype < Nbtype; ++btype) {
       test(btype, tessellator);
     }
   }
   if (Communicator::getRank() == root) {
     std::cout << "=== DistributedVoronoi2D passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
