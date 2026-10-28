#include <cmath>
#include <exception>
#include <iostream>
#include <vector>

#include "polytope.hh"

#include "BoostTessellator.hh"
#include "PLC.hh"
#include "Tessellation.hh"
#include "polytope_test_utilities.hh"
#include "SiloUtils.hh"
#include "SiloReader.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif
#ifdef POLYTOPE_ENABLE_MPI
#include "DistributedTessellator.hh"
#endif

using namespace polytope;

namespace {

std::vector<double>
generatorPoints() {
  return {0.20, 0.20,
          0.45, 0.15,
          0.75, 0.20,
          0.15, 0.45,
          0.40, 0.40,
          0.65, 0.45,
          0.85, 0.55,
          0.20, 0.75,
          0.50, 0.70,
          0.75, 0.80,
          0.35, 0.90,
          0.60, 0.25};
}

std::vector<double>
rankLocalPoints(const std::vector<double>& allPoints,
                const int rank,
                const int size) {
  std::vector<double> result;
  const auto n = allPoints.size()/2;
  if (rank < n) {
    for (unsigned i = 0; i < n; ++i) {
      if (static_cast<int>(i % size) == rank) {
        result.push_back(allPoints[2*i]);
        result.push_back(allPoints[2*i + 1]);
      }
    }
  }
  return result;
}

void test(Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int nranks = Communicator::getNRanks();

  const auto allPoints = generatorPoints();
  const auto localPoints = rankLocalPoints(allPoints, rank, nranks);
  const auto expectedCells = allPoints.size()/2;
  Boundary2D boundary;
  boundary.mClipping = false;
  boundary.mPad = 0.04;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(0);
  tessellator.setPadding(boundary.mPad);

  std::string serial_mesh_name = "serialIO" + std::to_string(nranks) + tessellator.name();
  std::string read_mesh_name = serial_mesh_name;
  unsigned totalCells = 0;
  Tessellation<2, double> serialMesh;
  if (rank == Communicator::getRoot()) {
    serialMesh = tessellator.tessellate(allPoints);
    totalCells = serialMesh.totalPointCount(false);
    outputMesh(serialMesh, serial_mesh_name, 0, 0., 1);
    compareArea(boundary, serialMesh, tessellator.getQuantizer(),
                "Serial area failed", false);
  }
  if (nranks > 1) {
    DistributedTessellator<2> distributed(tessellator);
    Tessellation<2, double> localMesh = distributed.tessellate(localPoints);
    totalCells = localMesh.totalPointCount(true);

    std::string parallel_mesh_name = "parallelIO" + std::to_string(nranks) + tessellator.name();
    outputMesh(localMesh, parallel_mesh_name, 0, 0.);
    read_mesh_name = parallel_mesh_name;
    compareArea(boundary, localMesh, distributed.getQuantizer(), "Distributed area failed");
  }
  POLY_CHECK2(totalCells == expectedCells,
              "Distributed output has " << totalCells
              << " total cells but expected " << expectedCells);
  // Now try to open the file we just created
  Tessellation<2, double> readMesh;
  std::string masterFilename = getMasterFilename(read_mesh_name, 0);
  SiloReader<2, Tessellation<2, double>>::FieldTypeMap fields;
  SiloReader<2, Tessellation<2, double>>::read(readMesh, fields, masterFilename);
  totalCells = readMesh.totalPointCount(true);
  POLY_CHECK2(totalCells == expectedCells, "Read in mesh has " << totalCells
              << " cells but expected " << expectedCells);
}
} // anonymous namespace

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);

#ifdef POLYTOPE_ENABLE_TRIANGLE
   {
     if (Communicator::getRank() == 0) {
       cout << "\nTriangle Tessellator:\n" << endl;
     }
     TriangleTessellator tessellator;
     test(tessellator);
   }
#endif

   {
     if (Communicator::getRank() == 0) {
       cout << "\nBoost Tessellator:\n" << endl;
     }
     BoostTessellator tessellator;
     test(tessellator);
   }
   if (Communicator::getRank() == 0) {
     std::cout << "=== SiloIO2D passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
