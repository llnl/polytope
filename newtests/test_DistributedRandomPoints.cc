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

  const int Ngen = 4;
  int oseed = 1049600 + 10*btype;
  int seed = oseed + rank;
  Boundary2D boundary;
  boundary.setDefaultBoundary(btype);
  Generators<2> generators(boundary);
  generators.randomPoints(Ngen, seed);

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

  double serialArea = 0.0;
  Tessellation<2, double> serialMesh;
  if (rank == root) {
    std::vector<double> allPoints;
    for (int p = 0; p < nprocs; ++p) {
      int cseed = oseed + p;
      generators.randomPoints(Ngen, cseed);
      std::copy(generators.mPoints.begin(), generators.mPoints.end(), std::back_inserter(allPoints));
    }
    tessellator.tessellate(allPoints, boundary.mPLCpoints, boundary.mPLC, serialMesh);
    serialArea = computeTessellationArea(serialMesh);
  }
  std::string serialoutname = "Serial_" + tessellator.name();
  outputMesh(serialMesh, serialoutname, btype, 0.);
  std::string outname = "DistributedRandom_" + tessellator.name();
  outputMesh(localMesh, outname, btype, 0.);
  MPI_Bcast(&serialArea, 1, MPI_DOUBLE, 0, Communicator::communicator());

  POLY_CHECK2(std::abs(distributedArea - serialArea) < 1.0e-8,
              "Distributed area " << distributedArea
              << " differs from serial area " << serialArea);
}
} // anonymous namespace

int main(int argc, char** argv) {
  auto& comm = Communicator::instance();
  comm.init(argc, argv);
  const int root = Communicator::getRoot();
  const int Nbtype = 11;

#ifdef POLYTOPE_ENABLE_TRIANGLE
   {
     if (Communicator::getRank() == root) {
       cout << "\nTriangle Tessellator:\n" << endl;
     }
     TriangleTessellator tessellator;
     for (int btype = 0; btype < Nbtype; ++btype) {
       test(btype, tessellator);
     }
   }
#endif

   {
     if (Communicator::getRank() == root) {
       cout << "\nBoost Tessellator:\n" << endl;
     }
     BoostTessellator tessellator;
     for (int btype = 0; btype < Nbtype; ++btype) {
       test(btype, tessellator);
     }
   }
   if (Communicator::getRank() == root) {
     std::cout << "=== DistributedVoronoi2D passed ===" << std::endl;
   }
  comm.finalize();
  return 0;
}
