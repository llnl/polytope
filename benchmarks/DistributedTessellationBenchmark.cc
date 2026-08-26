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

namespace {

void
timeTessellation(const unsigned n, Tessellator<2, double>& serialTessellator) {
  const auto rank = Communicator::getRank();
  const auto root = Communicator::getRoot();
  const auto nranks = static_cast<unsigned>(Communicator::getNProcs());

  // Inputs for the LatticePartitioner
  const auto r0 = static_cast<unsigned>(std::sqrt(nranks));
  const auto r1 = nranks/r0;
  POLY_VERIFY2(r1 >= 2,
               "DistributedTessellationBenchmark requires r1 >= 2");

  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(Boundary2D::square);
  Quantizer<2>::instance().init(Point<2, double>(0.0, 0.0),
                                Point<2, double>(1.0, 1.0));

  if (rank == root) {
    std::cout << "Degeneracy " << Quantizer<2>::instance().degeneracy()
              << "\nGenerating " << n << " random points\n";
  }

  const auto generationStart = std::chrono::steady_clock::now();
  Generators<2> generators(boundary);
  generators.randomPoints(n, 1049600);
  const auto generationTime = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - generationStart).count();

  if (rank == root) {
    std::cout << "Point generation: " << generationTime << " s\n"
              << "Doing " << serialTessellator.name() << " tessellation\n";
  }

  DistributedTessellator<2> distributed(serialTessellator);
  // Add other partitioners as they are implemented
  std::vector<LatticePartitioner<2>> partitioners = {LatticePartitioner<2>({r0, r1})};

  for (unsigned exchangePoints = 0; exchangePoints != 2; ++exchangePoints) {
    for (const auto& partitioner : partitioners) {
      Tessellation<2, double> mesh;
      Communicator::Barrier();
      const auto tessellationStart = std::chrono::steady_clock::now();
      distributed.partitionAndTessellate(generators.mPoints, partitioner,
                                         mesh);
      Communicator::Barrier();
      const auto tessellationTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - tessellationStart).count();

      if (rank == root) {
        std::cout << "exchange points: " << distributed.exchangePoints()
                  << ", " << partitioner.name()
                  << ": " << tessellationTime << " s\n";
      }
    }
    distributed.setExchangePoints(true);
  }
}

} // anonymous namespace

int
main(int argc, char** argv) {
  auto& communicator = Communicator::instance();
  communicator.init(argc, argv);
  const auto nranks = Communicator::getNProcs();

  unsigned n = 1000000;
  if (argc > 1) n = std::strtoul(argv[1], nullptr, 10);

  BoostTessellator boost;
  timeTessellation(n, boost);

#ifdef POLYTOPE_ENABLE_TRIANGLE
  TriangleTessellator triangle;
  timeTessellation(n, triangle);
#endif

  communicator.finalize();
  return 0;
}
