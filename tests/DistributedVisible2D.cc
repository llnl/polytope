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

std::vector<std::vector<double>>
rankPoints = {{0.0563102, 0.480787,
               0.269074, 0.251452},
              {0.54725, 0.532511,
               0.631279, 0.512978,
               0.425365, 0.542825,
               0.643499, 0.843612},
              {0.187156, 0.989748},
              {0.813925, 0.024177,
               0.945232, 0.272313,
               0.971435, 0.214536},
              {0.472844, 0.87031,
               0.211005, 0.71287,
               0.0202204, 0.814865},
              {0.609406, 0.220088,
               0.564849, 0.419905,
               0.787302, 0.435326}};

void test(Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  int nranks = Communicator::getNRanks();

  const int activeRanks = 6;
  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(0);
  boundary.mClipping = false;
  std::vector<double> points;
  if (rank < activeRanks) {
    points = rankPoints[rank];
  }
  tessellator.setPadding(boundary.mPad);

  // Number of points per rank
  unsigned refCellCount = 0;
  std::vector<unsigned> finalProcN(activeRanks, 0);
  for (int i = 0; i < activeRanks; ++i) {
    finalProcN[i] = rankPoints[i].size()/2;
    refCellCount += finalProcN[i];
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
  localMesh = distributed.tessellate(points);

  unsigned totalCells = localMesh.totalPointCount(true);
  POLY_CHECK2(totalCells == refCellCount, "Incorrect total number of points");

  std::string outname = "DistributedVisible_" + tessellator.name();
  outputMesh(localMesh, outname, 0, 0.);

  // Now get the current mesh including all it's neighbor generators
  auto& quantizer = distributed.quantizer();
  QuantTessellation<2> qmesh(quantizer, points);
  distributed.tessellateQuantized(qmesh);
  if (rank < activeRanks && nranks >= activeRanks) {
    Tessellation<2, double> procMesh;
    qmesh.fillTessellation(procMesh);
    SiloWriter<2, Tessellation<2, double>> writer(procMesh);
    writer.generateTestVars();
    writer.write("ProcMesh" + tessellator.name(), "", rank, 0., 1);
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
