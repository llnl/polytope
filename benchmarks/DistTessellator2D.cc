// DistributedTessellationBenchmark
//
// Time the partitioned distributed 2-D tessellation with both generator
// exchange representations. Example of how to run:
//
//   mpirun -n 81 ./benchmarks/DistributedTessellationBenchmark 1000000

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "polytope_test_utilities.hh"
#include "polytope.hh"

#include "BoostTessellator.hh"
#include "Communicator.hh"
#include "DistributedTessellator.hh"
#include "Generators.hh"
#include "Partitioner.hh"
#include "LatticePartitioner.hh"
#include "QuasiVoronoiPartitioner.hh"
#include "Tessellation.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace polytope;

void
timeTessellation(const std::vector<double> allPoints,
                 Tessellator<2, double>& serialTessellator,
                 std::string distType) {
  const auto rank = Communicator::getRank();
  const auto root = Communicator::getRoot();
  const int partseed = 1042390;

  // Add other partitioners as they are implemented
  std::vector<std::unique_ptr<Partitioner<2>>> partitioners;
  partitioners.push_back(std::make_unique<LatticePartitioner<2>>());
  partitioners.push_back(std::make_unique<QuasiVoronoiPartitioner<2>>(partseed));
  double minTime = std::numeric_limits<double>::max();
  std::vector<std::string> run_info;
  for (unsigned encoding = 0; encoding < 2; ++encoding) {
    if (encoding == 0) serialTessellator.useMortonEncoding();
    else serialTessellator.usePackedEncoding();
    DistributedTessellator<2> distributed(serialTessellator);    
    const auto& quantizer = distributed.quantizer();
    if (rank == root) {
      std::cout << "Encoding: " << quantizer.keyName() << std::endl;
    }
    for (unsigned exchangePoints = 0; exchangePoints < 2; ++exchangePoints) {
      if (exchangePoints == 0) distributed.setExchangePoints(false);
      else distributed.setExchangePoints(true);
      if (rank == root) {
        std::cout << "MPI Exchange: " << ((distributed.exchangePoints()) ? "points" : "hashes") << std::endl;
      }
      for (const auto& partitioner : partitioners) {
        Tessellation<2, double> mesh;
        Communicator::Barrier();
        const auto tessellationStart = std::chrono::steady_clock::now();
        mesh = distributed.partitionAndTessellate(allPoints, *(partitioner));
        Communicator::Barrier();
        const auto tessellationTime = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - tessellationStart).count();

        if (rank == root) {
          std::cout << partitioner->name() << " " << tessellationTime << " s\n";
        }
        if (tessellationTime < minTime) {
          minTime = tessellationTime;
          run_info.clear();
          run_info.push_back(quantizer.keyName());
          run_info.push_back((distributed.exchangePoints()) ? "Points" : "Hashes");
          run_info.push_back(partitioner->name());
        }
        if (encoding == 0 && exchangePoints == 0) {
          std::string meshName = "Bench" + quantizer.keyName() + "_" + partitioner->name();
          outputMesh(mesh, meshName);
        }
      } // Partitioner loop
    } // Exchange type loop
  } // Encoding type loop
  if (rank == root) {
    std::cout << "Fastest time for " << distType << " point distribution and " << serialTessellator.name() << "\n";
    for (auto& st : run_info) {
      std::cout << st << ", ";
    }
    std::cout << std::endl << minTime << std::endl;
  }
}

int
main(int argc, char** argv) {
  auto& communicator = Communicator::instance();
  communicator.init(argc, argv);
  unsigned n = 1000000;
  if (argc > 1) n = std::strtoul(argv[1], nullptr, 10);


  const auto rank = Communicator::getRank();
  const auto root = Communicator::getRoot();
  const auto nranks = Communicator::getNRanks();
  for (int dist_type = 0; dist_type < 2; ++dist_type) {
    if (rank == root) {
      std::cout << "Generating " << n << " random points over " << nranks << " ranks\n";
    }
#ifdef POLYTOPE_ENABLE_TRIANGLE
    TriangleTessellator tessellator;
#else
    BoostTessellator tessellator;
#endif
    if (rank == root) {
      std::cout << "Using " << tessellator.name() << " tessellator\n";
    }
    // Setup the quantizer
    Boundary2D boundary;
    boundary.mCenter[0] = 0.5;
    boundary.mCenter[1] = 0.5;
    boundary.setDefaultBoundary(Boundary2D::square);
    const int genseed = 1049600;
    Generators<2> generators(boundary);
    std::string distType;
    if (dist_type == 0) {
      generators.randomPoints(n, genseed);
      distType = "Uniform";
    } else {
      generators.randomNormalPoints(n, genseed);
      distType = "Normal";
    }
    if (rank == root) {
      std::cout << distType << " point distribution\n";
    }
    timeTessellation(generators.mPoints, tessellator, distType);
  }
  communicator.finalize();
  return 0;
}
