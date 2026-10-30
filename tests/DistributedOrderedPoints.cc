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
#include "QuasiVoronoiPartitioner.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace polytope;

namespace {

void test(const int btype, Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int nranks = Communicator::getNRanks();
  int root = Communicator::getRoot();

  const int Ngen = 4*nranks;
  int seed = 1049600;
  Boundary2D boundary;
  boundary.setDefaultBoundary(btype);
  Generators<2> generators(boundary);
  generators.randomPoints(Ngen, seed);
  std::vector<double>& allPoints = generators.mPoints;
  QuasiVoronoiPartitioner<2> qvpart(seed);

  Tessellation<2, double> serialMesh;
  if (rank == root) {
    serialMesh = tessellator.tessellate(allPoints, boundary.mPLCpoints, boundary.mPLC);
    std::string serialoutname = "SerialOrdered_" + tessellator.name();
    qvpart.assignCellRanks(serialMesh);
    SiloWriter<2, Tessellation<2, double>> serial_writer(serialMesh);
    // Output the serial mesh with a rank variable showing the ranks for the
    // distributed test
    serial_writer.write(serialoutname, btype, 0., 1);
    compareArea(boundary, serialMesh, tessellator.quantizer(),
                "Serial area failure", false);
  }
  Communicator::Barrier();
  DistributedTessellator<2> distributed(tessellator);
  auto localPoints = qvpart.computeLocalPartition(allPoints);
  Tessellation<2, double> localMesh =
    distributed.tessellate(localPoints, boundary.mPLCpoints, boundary.mPLC);
  std::string outname = "DistributedOrdered_" + tessellator.name();
  outputMesh(localMesh, outname, btype, static_cast<double>(btype));

  auto totalCells = localMesh.totalPointCount(true);
  POLY_CHECK2(totalCells == allPoints.size()/2,
              "Total number of cells " << totalCells
              << " is incorrect; expected " << allPoints.size()/2);

  compareArea(boundary, localMesh, distributed.quantizer(), "Distributed area failure");
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
     std::cout << "=== DistributedOrderedPoints passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
