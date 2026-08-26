#include <cmath>
#include <exception>
#include <iostream>
#include <vector>

#include "polytope.hh"
#include "mpi.h"

#include "BoostTessellator.hh"
#include "DistributedTessellator.hh"
#include "PLC.hh"
#include "Tessellation.hh"
#include "polytope_test_utilities.hh"
#include "SiloUtils.hh"
#include "SiloReader.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace polytope;

namespace {

PLC<2>
unitSquarePLC() {
  PLC<2> result;
  result.facets = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  return result;
}

std::vector<double>
unitSquarePoints() {
  return {0.0, 0.0,
          1.0, 0.0,
          1.0, 1.0,
          0.0, 1.0};
}

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
  for (unsigned i = 0; i < n; ++i) {
    if (static_cast<int>(i % size) == rank) {
      result.push_back(allPoints[2*i]);
      result.push_back(allPoints[2*i + 1]);
    }
  }
  return result;
}

void test(Tessellator<2, double>& tessellator) {
  int rank = Communicator::getRank();
  int size = Communicator::getNProcs();

  const auto allPoints = generatorPoints();
  const auto localPoints = rankLocalPoints(allPoints, rank, size);
  const auto boundaryPoints = unitSquarePoints();
  const auto boundary = unitSquarePLC();

  DistributedTessellator<2> distributed(tessellator);

  Tessellation<2, double> localMesh;
  distributed.tessellate(localPoints, localMesh);

  auto localCells = static_cast<int>(localMesh.cells.size());
  int totalCells = 0;
  MPI_Allreduce(&localCells, &totalCells, 1, MPI_INT, MPI_SUM, Communicator::communicator());

  const auto expectedCells = static_cast<int>(allPoints.size()/2);
  POLY_CHECK2(totalCells == expectedCells,
              "Distributed output has " << totalCells
              << " total cells but expected " << expectedCells);
  POLY_CHECK2(localCells == static_cast<int>(localPoints.size()/2),
              "Rank " << rank << " output " << localCells
              << " cells for " << localPoints.size()/2 << " owned generators");

  double localArea = computeTessellationArea(localMesh);
  double distributedArea = 0.0;
  MPI_Allreduce(&localArea, &distributedArea, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());

  double serialArea = 0.0;
  if (rank == 0) {
    Tessellation<2, double> serialMesh;
    tessellator.tessellate(allPoints, serialMesh);
    serialArea = computeTessellationArea(serialMesh);
  }
  MPI_Bcast(&serialArea, 1, MPI_DOUBLE, 0, Communicator::communicator());

  POLY_CHECK2(std::abs(distributedArea - serialArea) < 1.0e-8,
              "Distributed area " << distributedArea
              << " differs from serial area " << serialArea);

  std::string outname = "parallelIO_" + tessellator.name();
  outputMesh(localMesh, outname, 0, 0.);

  // Now try to open the file we just created
  Tessellation<2, double> readMesh;
  std::string masterFilename = getMasterFilename(outname, 0);
  SiloReader<2, Tessellation<2, double>>::FieldTypeMap fields;
  SiloReader<2, Tessellation<2, double>>::read(readMesh, fields, masterFilename);
  localCells = static_cast<int>(readMesh.cells.size());
  totalCells = 0;
  MPI_Allreduce(&localCells, &totalCells, 1, MPI_INT, MPI_SUM, Communicator::communicator());
  POLY_CHECK2(totalCells == expectedCells,
              "Read in mesh has " << totalCells
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
     std::cout << "=== DistributedVoronoi2D passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
