#include <array>
#include <cstdint>
#include <vector>

#include "Communicator.hh"
#include "Partitioner.hh"
#include "polytope_internal.hh"

using namespace polytope;

namespace {

using PointType = QuantizedPoint<2>;

std::vector<PointType>
generators() {
  return {
    PointType(0, 0), PointType(49, 49),
    PointType(50, 10), PointType(100, 40),
    PointType(10, 50), PointType(40, 100),
    PointType(50, 50), PointType(100, 100)
  };
}

void testSerial() {
  const auto points = generators();
  RandomPartitioner<2> random(123456789ULL);
  LatticePartitioner<2> lattice(PointType(0, 0), PointType(100, 100), {1, 1});
  POLY_CHECK(random.computePartition(points) == points);
  POLY_CHECK(lattice.computePartition(points) == points);
}

#ifdef POLYTOPE_ENABLE_MPI
void testDistributed() {
  const auto points = generators();
  const auto rank = Communicator::getRank();

  LatticePartitioner<2> lattice(PointType(0, 0), PointType(100, 100), {2, 2});
  const std::vector<std::vector<PointType>> expected = {
    {points[0], points[1]}, {points[2], points[3]},
    {points[4], points[5]}, {points[6], points[7]}
  };
  const auto localLattice = lattice.computePartition(points);
  POLY_CHECK(localLattice == expected[rank]);

  int latticeCount = static_cast<int>(localLattice.size());
  int globalLatticeCount = 0;
  MPI_Allreduce(&latticeCount, &globalLatticeCount, 1, MPI_INT, MPI_SUM,
                Communicator::communicator());
  POLY_CHECK(globalLatticeCount == static_cast<int>(points.size()));

  RandomPartitioner<2> first(123456789ULL), second(987654321ULL);
  const auto localFirst = first.computePartition(points);
  POLY_CHECK(localFirst == first.computePartition(points));

  int randomCount = static_cast<int>(localFirst.size());
  int globalRandomCount = 0;
  MPI_Allreduce(&randomCount, &globalRandomCount, 1, MPI_INT, MPI_SUM,
                Communicator::communicator());
  POLY_CHECK(globalRandomCount == static_cast<int>(points.size()));

  const int locallyDifferent = localFirst != second.computePartition(points) ? 1 : 0;
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
