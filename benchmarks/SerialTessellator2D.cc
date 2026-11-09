// SerialTessellationBenchmark
//
// Time the serial 2-D tessellation
// Usage:
// ./benchmarks/SerialTessellator2D Nr Np
// where Nr is the number of times to run each
// tessellator and Np is the number of generator points
// By default Nr = 5 and Np = 50000

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
                 Tessellator<2, double>& tessellator,
                 unsigned nruns) {
  Tessellation<2, double> mesh;
  const auto tessellationStart = std::chrono::steady_clock::now();
  for (auto i = 0; i < nruns; ++i) {
    Tessellation<2, double> mesh = tessellator.tessellate(allPoints);
  }
  const auto tessellationTime = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - tessellationStart).count();
  std::cout << tessellator.name() << " " << tessellationTime << " s\n";
}

int
main(int argc, char** argv) {
  auto& communicator = Communicator::instance();
  communicator.init(argc, argv);
  unsigned n = 50000; // Number of generator points
  unsigned nruns = 5; // Number of times to run each tessellator
  if (argc > 1) {
    nruns = std::strtoul(argv[1], nullptr, 10);
    if (argc > 2) n = std::strtoul(argv[2], nullptr, 10);
  }
  // Create generator points
  Boundary2D boundary;
  boundary.mCenter[0] = 0.5;
  boundary.mCenter[1] = 0.5;
  boundary.setDefaultBoundary(Boundary2D::square);
  const int genseed = 1049600;
  Generators<2> generators(boundary);
  generators.randomNormalPoints(n, genseed);

  {
    BoostTessellator tessellator;
    timeTessellation(generators.mPoints, tessellator, nruns);
  }

#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    TriangleTessellator tessellator;
    timeTessellation(generators.mPoints, tessellator, nruns);
  }
#endif
  communicator.finalize();
  return 0;
}
