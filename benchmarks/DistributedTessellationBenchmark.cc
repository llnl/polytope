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
#include "Tessellation.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

using namespace polytope;

struct TessTest {
public:
  std::string encodingType; // Type of key encoder used (Morton or Packed)
  std::string distType; // Distribution used for creating generators (Normal or Uniform)
  std::string exchangeType; // What gets exchanged between MPI ranks (hashes or points)
  std::string partitionerType; // Type of partitioner used to distribute generators over ranks (Lattice or QuasiVoronoi)
  std::string tessType; // Type of tessellator used (Boost or Triangle)
  double degeneracy = 0.;
  double tessTime; // Time needed to do tessellation
};

std::vector<TessTest>
timeTessellation(const std::vector<double> allPoints,
                 Tessellator<2, double>& serialTessellator,
                 std::string distType) {
  const auto rank = Communicator::getRank();
  const auto root = Communicator::getRoot();
  const auto nranks = static_cast<unsigned>(Communicator::getNProcs());
  std::vector<TessTest> output;

  // Inputs for the LatticePartitioner
  const auto r0 = static_cast<unsigned>(std::sqrt(nranks));
  const auto r1 = static_cast<unsigned>(nranks/r0);
  std::array<unsigned, 2> rarray = {r0, r1};
  const int partseed = 1042390;

  DistributedTessellator<2> distributed(serialTessellator);
  // Add other partitioners as they are implemented
  std::vector<std::unique_ptr<Partitioner<2>>> partitioners;
  partitioners.push_back(std::make_unique<LatticePartitioner<2>>(rarray));
  partitioners.push_back(std::make_unique<QuasiVoronoiPartitioner<2>>(partseed));
  for (unsigned encoding = 0; encoding < 2; ++encoding) {
    if (encoding == 0) Quantizer<2>::instance().useMortonEncoding();
    else Quantizer<2>::instance().usePackedEncoding();
    if (rank == root) {
      std::cout << "Encoding: " << Quantizer<2>::instance().keyName() << std::endl;
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
        distributed.partitionAndTessellate(allPoints, *(partitioner),
                                           mesh);
        Communicator::Barrier();
        const auto tessellationTime = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - tessellationStart).count();

        if (rank == root) {
          std::cout << partitioner->name() << " " << tessellationTime << " s\n";
        }
        TessTest tt;
        tt.distType = distType;
        tt.encodingType = Quantizer<2>::instance().keyName();
        tt.exchangeType = (distributed.exchangePoints()) ? "Points" : "Hashes";
        tt.partitionerType = partitioner->name();
        tt.tessType = serialTessellator.name();
        tt.tessTime = tessellationTime;
        output.push_back(tt);
        if (encoding == 0 && exchangePoints == 0) {
          std::string meshName = "Bench" + Quantizer<2>::instance().keyName() + "_" + partitioner->name();
          outputMesh(mesh, meshName);
        }
      } // Partitioner loop
    } // Exchange type loop
  } // Encoding type loop
  return output;
}

int
main(int argc, char** argv) {
  auto& communicator = Communicator::instance();
  communicator.init(argc, argv);
  unsigned n = 1000000;
  if (argc > 1) n = std::strtoul(argv[1], nullptr, 10);
  std::vector<TessTest> test_vals;

  // Setup the quantizer
  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(Boundary2D::square);
  Quantizer<2>::instance().init(Point<2, double>(0.0, 0.0),
                                Point<2, double>(1.0, 1.0));
  const auto rank = Communicator::getRank();
  const auto root = Communicator::getRoot();
  const auto nranks = Communicator::getNProcs();
  for (int dist_type = 0; dist_type < 2; ++dist_type) {
    if (rank == root) {
      std::cout << "Degeneracy " << Quantizer<2>::instance().degeneracy()
                << "\nGenerating " << n << " random points over " << nranks << " ranks\n";
    }
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
      std::cout << distType << " point distribution\n"
                << "-----------------------------------------\n"
                << "Using Boost tessellator\n";
    }
    BoostTessellator boost;
    auto btests = timeTessellation(generators.mPoints, boost, distType);
    std::copy(btests.begin(), btests.end(), std::back_inserter(test_vals));

#ifdef POLYTOPE_ENABLE_TRIANGLE
    if (rank == root) {
      std::cout << "Using Triangle tessellator\n";
    }
    TriangleTessellator triangle;
    auto ttests = timeTessellation(generators.mPoints, triangle, distType);
    std::copy(ttests.begin(), ttests.end(), std::back_inserter(test_vals));
#endif
  }
  communicator.finalize();
  return 0;
}
