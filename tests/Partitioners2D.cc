#include <array>
#include <cstdint>
#include <vector>

#include "Communicator.hh"
#include "Partitioner.hh"
#include "QuantTessellation.hh"
#include "Quantizer.hh"
#include "Tessellation.hh"
#include "polytope_internal.hh"

using namespace polytope;

namespace {

using PointType = QuantizedPoint<2>;

std::vector<PointType>
generators() {
  const auto& Q = Quantizer<2>::instance();
  const auto& lower = Q.minBound;
  const auto& upper = Q.maxBound;
  const auto xmid = lower.x + (upper.x - lower.x + 1)/2;
  const auto ymid = lower.y + (upper.y - lower.y + 1)/2;
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
  QuasiVoronoiPartitioner<2> quasiVoronoi(123456789U);
  LatticePartitioner<2> lattice({1, 1});
  POLY_CHECK(random.computeLocalPartition(points) == points);
  POLY_CHECK(quasiVoronoi.computeLocalPartition(points) == points);
  POLY_CHECK(lattice.computeLocalPartition(points) == points);
  std::vector<Point<2, double>> realGeneratorPoints;
  realGeneratorPoints.reserve(points.size());
  const auto& Q = Quantizer<2>::instance();
  for (const auto& point : points) realGeneratorPoints.push_back(Q.dequantize(point));
  POLY_CHECK(lattice.computeLocalPartition(realGeneratorPoints) == realGeneratorPoints);

  const auto randomResult = random.computePartitionResult(points);
  POLY_CHECK(randomResult.ownerByInput.size() == points.size());
  POLY_CHECK(randomResult.generatorsByPartition.size() == 1u);
  POLY_CHECK(randomResult.generatorsByPartition[0] == points);

  // Serial output can label cells using an arbitrary number of logical
  // partitions, independent of the number of MPI ranks.
  QuasiVoronoiPartitioner<2> virtualRanks(123456789U, 4);
  const auto virtualResult = virtualRanks.computePartitionResult(points);
  POLY_CHECK(virtualResult.generatorsByPartition.size() == 4u);
  POLY_CHECK(virtualResult.ownerByInput.size() == points.size());
  for (const auto owner : virtualResult.ownerByInput) {
    POLY_CHECK(owner < 4u);
  }
  QuantTessellation<2> quantmesh(points);
  virtualRanks.assignCellRanks(quantmesh);
  const auto meshOwners = virtualRanks.computeOwners(quantmesh.points);
  POLY_CHECK(quantmesh.cellRank.size() == meshOwners.size());
  for (std::size_t i = 0; i < meshOwners.size(); ++i) {
    POLY_CHECK(quantmesh.cellRank[i] == static_cast<int>(meshOwners[i]));
  }
  Tessellation<2, double> serialMesh;
  serialMesh.points.resize(quantmesh.points.size());
  for (std::size_t i = 0; i < quantmesh.points.size(); ++i) {
    serialMesh.points[i] = Q.dequantize(quantmesh.points[i]);
  }
  virtualRanks.assignCellRanks(serialMesh);//, points);
  POLY_CHECK(serialMesh.cellRank == quantmesh.cellRank);
  Tessellation<2, double> clippedMesh;
  clippedMesh.points = {serialMesh.points[0]};
  virtualRanks.assignCellRanks(clippedMesh);//, points);
  POLY_CHECK(clippedMesh.cellRank == std::vector<int>{serialMesh.cellRank[0]});

  // Exercise the real-valued convenience interface as well as the quantized
  // one above.  With one rank, all generators must remain local.
  const std::vector<double> realPoints = {0.25, 0.25, 0.75, 0.75};
  POLY_CHECK(random.computeLocalPartition(realPoints).size() == 2u);
  POLY_CHECK(quasiVoronoi.computeLocalPartition(realPoints).size() == 2u);
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

  // Quasi-Voronoi must assign every point to exactly one rank.  Also check
  // the all-ranks interface agrees with the local convenience interface.
  QuasiVoronoiPartitioner<2> quasiVoronoi(123456789U);
  const auto allQuasiVoronoi = quasiVoronoi.computePartition(points);
  POLY_CHECK(allQuasiVoronoi.size() == static_cast<std::size_t>(Communicator::getNProcs()));
  const auto localQuasiVoronoi = quasiVoronoi.computeLocalPartition(points);
  POLY_CHECK(localQuasiVoronoi == allQuasiVoronoi[rank]);
  const auto quasiVoronoiResult = quasiVoronoi.computePartitionResult(points);
  POLY_CHECK(quasiVoronoiResult.generatorsByPartition == allQuasiVoronoi);
  POLY_CHECK(quasiVoronoiResult.ownerByInput.size() == points.size());
  int quasiVoronoiCount = static_cast<int>(localQuasiVoronoi.size());
  int globalQuasiVoronoiCount = 0;
  MPI_Allreduce(&quasiVoronoiCount, &globalQuasiVoronoiCount, 1, MPI_INT, MPI_SUM,
                Communicator::communicator());
  POLY_CHECK(globalQuasiVoronoiCount == static_cast<int>(points.size()));

  // When there are fewer points than ranks, the short-input path gives the
  // first two ranks one point each and leaves the remaining ranks empty.
  const std::vector<PointType> shortPoints = {points[0], points[1]};
  const auto allShortQuasiVoronoi = quasiVoronoi.computePartition(shortPoints);
  POLY_CHECK(allShortQuasiVoronoi.size() == static_cast<std::size_t>(Communicator::getNProcs()));
  POLY_CHECK(allShortQuasiVoronoi[0] == std::vector<PointType>{points[0]});
  POLY_CHECK(allShortQuasiVoronoi[1] == std::vector<PointType>{points[1]});
  for (std::size_t i = 2; i < allShortQuasiVoronoi.size(); ++i) {
    POLY_CHECK(allShortQuasiVoronoi[i].empty());
  }
  const auto localShortQuasiVoronoi = quasiVoronoi.computeLocalPartition(shortPoints);
  POLY_CHECK(localShortQuasiVoronoi == allShortQuasiVoronoi[rank]);

  // Empty input must be safe for every local rank.
  const std::vector<PointType> noPoints;
  POLY_CHECK(quasiVoronoi.computePartition(noPoints).size() ==
             static_cast<std::size_t>(Communicator::getNProcs()));
  POLY_CHECK(quasiVoronoi.computeLocalPartition(noPoints).empty());
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

  if (Communicator::getRank() == Communicator::getRoot()) {
    std::cout << "Partitioner tests passed!" << std::endl;
  }
  communicator.finalize();
  return 0;
}
