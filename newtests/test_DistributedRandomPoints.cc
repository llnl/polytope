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
    auto& Q = Quantizer<2>::instance();
    std::map<Point<2, double>, int> finalRanks;
    std::vector<double> allPoints;
    for (int proc = 0; proc < nprocs; ++proc) {
      int cseed = oseed + proc;
      generators.randomPoints(Ngen, cseed);
      std::copy(generators.mPoints.begin(), generators.mPoints.end(), std::back_inserter(allPoints));
      for (auto i = 0u; i < generators.nPoints; ++i) {
        Point2<double> p;
        p.x = generators.mPoints[2*i];
        p.y = generators.mPoints[2*i+1];
        auto newp = Q.dequantize(Q.quantize(p));
        finalRanks[newp] = proc;
      }
    }
    tessellator.tessellate(allPoints, boundary.mPLCpoints, boundary.mPLC, serialMesh);
    serialArea = computeTessellationArea(serialMesh);
    std::vector<double> finalRanksD;
    for (auto& p : serialMesh.points) {
      auto rankItr = finalRanks.find(p);
      POLY_CHECK2(rankItr != finalRanks.end(),
                  "Could not find final rank for serial generator " << p);
      finalRanksD.push_back(static_cast<double>(rankItr->second));
    }
    std::string serialoutname = "SerialRandom_" + tessellator.name();
    outputMesh(serialMesh, serialoutname, finalRanksD, "rank", btype, 0., 1);
  }
  std::string outname = "DistributedRandom_" + tessellator.name();
  outputMesh(localMesh, outname, btype, 0.);
  MPI_Bcast(&serialArea, 1, MPI_DOUBLE, 0, Communicator::communicator());
  POLY_CHECK2(std::abs(boundary.mArea - serialArea)/boundary.mArea < 1.0e-8,
              "Serial area does not match reference area, ref: " << boundary.mArea
              << " serial: " << serialArea);

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
