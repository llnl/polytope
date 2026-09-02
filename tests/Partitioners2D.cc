#include <array>
#include <cstdint>
#include <vector>

#include "Communicator.hh"
#include "Partitioner.hh"
#include "Quantizer.hh"
#include "polytope_internal.hh"

using namespace polytope;

namespace {

using PointType = QuantizedPoint<2>;

std::vector<PointType>
generators() {
  const auto& Q = Quantizer<2>::instance();
  const auto& lower = Q.minBound;
  const auto& upper = Q.maxBound;
  const auto xmid = lower.x + (upper.x - lower.x)/2;
  const auto ymid = lower.y + (upper.y - lower.y)/2;
  return {
    lower, PointType(lower.x + 100, lower.y + 100),
    PointType(xmid, lower.y + 100), PointType(upper.x, lower.y + 100),
    PointType(lower.x + 100, ymid), PointType(lower.x + 100, upper.y),
    PointType(xmid, ymid), upper
  };
}

void initializeQuantizer() {
  Quantizer<2>::instance().init(Point<2, double>(0.0, 0.0),
                                Point<2, double>(1.0, 1.0));
}

void testSerial() {
  const auto points = generators();
  RandomPartitioner<2> random(123456789ULL);
  LatticePartitioner<2> lattice({1, 1});
  POLY_CHECK(random.computeLocalPartition(points) == points);
  POLY_CHECK(lattice.computeLocalPartition(points) == points);
}

#ifdef POLYTOPE_ENABLE_MPI
void testDistributed() {
  const auto points = generators();
  const auto rank = Communicator::getRank();

  LatticePartitioner<2> lattice({2, 2});
  const std::vector<std::vector<PointType>> expected = {
    {points[0], points[1]}, {points[2], points[3]},
    {points[4], points[5]}, {points[6], points[7]}
  };
  const auto localLattice = lattice.computeLocalPartition(points);
  POLY_CHECK(localLattice == expected[rank]);

  int latticeCount = static_cast<int>(localLattice.size());
  int globalLatticeCount = 0;
  MPI_Allreduce(&latticeCount, &globalLatticeCount, 1, MPI_INT, MPI_SUM,
                Communicator::communicator());
  POLY_CHECK(globalLatticeCount == static_cast<int>(points.size()));

  RandomPartitioner<2> first(123456789ULL), second(987654321ULL);
  const auto localFirst = first.computeLocalPartition(points);
  POLY_CHECK(localFirst == first.computeLocalPartition(points));

  int randomCount = static_cast<int>(localFirst.size());
  int globalRandomCount = 0;
  MPI_Allreduce(&randomCount, &globalRandomCount, 1, MPI_INT, MPI_SUM,
                Communicator::communicator());
  POLY_CHECK(globalRandomCount == static_cast<int>(points.size()));

  const int locallyDifferent = localFirst != second.computeLocalPartition(points) ? 1 : 0;
  int anyDifferent = 0;
  MPI_Allreduce(&locallyDifferent, &anyDifferent, 1, MPI_INT, MPI_MAX,
                Communicator::communicator());
  POLY_CHECK(anyDifferent == 1);
}
#endif

} // namespace

int main(int argc, char** argv) {
  auto& communicator = Communicator::instance();
  communicator.init(argc, argv);
  initializeQuantizer();

  const auto nranks = Communicator::getNProcs();
  if (nranks == 1) {
    testSerial();
  } else {
#ifdef POLYTOPE_ENABLE_MPI
    POLY_CHECK(nranks == 4);
    testDistributed();
#else
    POLY_CHECK(false);
#endif
  }

  communicator.finalize();
  return 0;
}
