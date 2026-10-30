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
  int nranks = Communicator::getNRanks();
  int root = Communicator::getRoot();

  const int Ngen = 100;
  int oseed = 1049600;
  int seed = oseed + rank;
  Boundary2D boundary;
  boundary.mClipping = false;
  boundary.setDefaultBoundary(0);
  tessellator.setPadding(boundary.mPad);
  Generators<2> generators(boundary);
  generators.randomPoints(Ngen, seed);

  DistributedTessellator<2> distributed(tessellator);

  Tessellation<2, double> localMesh;
  localMesh = distributed.tessellate(generators.mPoints);
  std::string outname = "parallelVoronoi_" + tessellator.name();
  outputMesh(localMesh, outname, 0, 0.);

  auto totalCells = localMesh.totalPointCount(true);
  POLY_CHECK2(totalCells = static_cast<unsigned>(Ngen*nranks),
              "Incorrect number of total cells");
  compareArea(boundary, localMesh,
              distributed.quantizer(),
              "Distributed area failure");

  Tessellation<2, double> serialMesh;
  if (rank == root) {
    std::vector<double> allPoints;
    for (int p = 0; p < nranks; ++p) {
      int cseed = oseed + p;
      generators.randomPoints(Ngen, cseed);
      std::copy(generators.mPoints.begin(), generators.mPoints.end(), std::back_inserter(allPoints));
    }
    serialMesh = tessellator.tessellate(allPoints);
    std::string serialoutname = "serialVoronoi_" + tessellator.name();
    outputMesh(serialMesh, serialoutname, 0, 0., 1);
  }
  compareArea(boundary, serialMesh,
              distributed.quantizer(),
              "Serial area failure");
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
