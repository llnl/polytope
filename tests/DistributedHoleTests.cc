// DistributedHoleTests
//
// Stress test of changing where points are and
// how they are distributed across ranks.
// Running with 8 ranks demonstrates issues with distributed clipping

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
#include "Partitioner.hh"

using namespace polytope;

namespace {

void test(Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  int nranks = Communicator::getNProcs();

  const int Ngen = 32;
  // Test case is square with a star hole
  Boundary2D boundary;
  boundary.setDefaultBoundary(10);
  // Number of different seeds for randomizing generator points
  const int NGseed = 10;
  // Number of different seeds for randomizing rank assignment
  const int NDseed = 10;
  int btype = 0;
  // Loop for randomizing the generator layout
  for (int Gseed = 0; Gseed < NGseed; ++Gseed) {
    Generators<2> generators(boundary);
    generators.randomPoints(Ngen, Gseed);
    std::vector<double> allPoints = generators.mPoints;
    // Loop for randomizing how generators are assigned to ranks
    for (int Dseed = 0; Dseed < NDseed; ++Dseed) {
      // Distribute points to specific ranks
      QuasiVoronoiPartitioner<2> qvpart(Dseed);
      // Solve the Voronoi in serial
      double serialArea = 0.0;
      if (rank == root) {
        Tessellation<2, double> serialMesh;
        tessellator.tessellate(allPoints, boundary.mPLCpoints, boundary.mPLC, serialMesh);
        serialArea = computeTessellationArea(serialMesh);
        qvpart.assignCellRanks(serialMesh);
        std::string serialoutname = "SerialHole";
        outputMesh(serialMesh, serialoutname, serialMesh.cellRank, "rank", btype, double(btype), 1);
      }
      Communicator::Barrier();
      DistributedTessellator<2> distributed(tessellator);
      Tessellation<2, double> localMesh;
      auto localPoints = qvpart.computeLocalPartition(allPoints);
      distributed.tessellate(localPoints, boundary.mPLCpoints, boundary.mPLC, localMesh);
      std::string outname = "DistributedHole";
      outputMesh(localMesh, outname, btype, static_cast<double>(btype));
      // Uncomment this section to output the mesh for each rank
      // QuantTessellation<2> qmesh(localPoints);
      // QuantPLC<2> qplc(boundary.mPLC, boundary.mPLCpoints);
      // distributed.tessellateQuantized(qmesh);
      // if (qmesh.points.size() > 0) {
      //   qmesh.clipTessellation(qplc, tessellator);
      //   Tessellation<2, double> procMesh;
      //   qmesh.fillTessellation(procMesh);
      //   SiloWriter<2, Tessellation<2, double>>::write(procMesh, "ProcMesh", 1, rank);
      // }
      compareArea(boundary, serialArea, "Serial area failure");
      compareArea(boundary, localMesh, "Distributed area failure");

      if (rank == root) {
        std::cout << btype << " test: gen seed " << Gseed << " partition seed " << Dseed << std::endl;
      }
      btype++;
    }
  }
}
} // anonymous namespace

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);
  const int root = Communicator::getRoot();

   {
     if (Communicator::getRank() == root) {
       cout << "\nBoost Tessellator:\n" << endl;
     }
     BoostTessellator tessellator;
     test(tessellator);
   }
   if (Communicator::getRank() == root) {
     std::cout << "=== DistributedHoleTests passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
