//----------------------------------------------------------------------------//
// LatticePartitioner
//
// The complete partitioning domain is Quantizer<Dimension>::minBound through
// Quantizer<Dimension>::maxBound, inclusive. The Quantizer must be initialized
// before computePartition is called. This domain is divided uniformly into
// ranksPerAxis[d] tiles along each axis; a point on the global upper bound
// belongs to the final tile on that axis.
//----------------------------------------------------------------------------//
#ifndef __Polytope_LatticePartitioner__
#define __Polytope_LatticePartitioner__

#include "Partitioner.hh"

namespace polytope {

template<int Dimension>
class LatticePartitioner: public Partitioner<Dimension> {
public:
  using PointType = QuantizedPoint<Dimension>;
  using OwnerType = typename Partitioner<Dimension>::OwnerType;
  using RealPoint = typename Partitioner<Dimension>::RealPoint;
  using QuantPoint = typename Partitioner<Dimension>::QuantPoint;
  using Partitioner<Dimension>::computeOwners;
  using RanksPerAxis = std::array<unsigned, Dimension>;

  virtual std::string name() const override {
    std::string dist = " [";
    for (int d = 0; d < Dimension; ++d) {
      dist += std::to_string(m_ranksPerAxis[d]) + ", ";
    }
    return "LatticePartitioner" + dist + "]";
  }

  explicit LatticePartitioner(const unsigned numPartitions = Communicator::getNRanks()):
    Partitioner<Dimension>(numPartitions) {
    // Try to distribute evenly across each dimension
    if constexpr (Dimension == 2) {
      optimalLattice2D();
    } else {
      optimalLattice3D();
    }
    init();
  }

  explicit LatticePartitioner(const RanksPerAxis& ranksPerAxis,
                              const unsigned numPartitions = Communicator::getNRanks()):
    Partitioner<Dimension>(numPartitions),
    m_ranksPerAxis(ranksPerAxis) {
    init();
  }

  void init() {
    std::size_t expectedRanks = 1;
    for (int d = 0; d < Dimension; ++d) {
      POLY_VERIFY2(m_ranksPerAxis[d] > 0,
                   "Each lattice axis must have at least one rank");
      POLY_VERIFY2(expectedRanks <= std::numeric_limits<std::size_t>::max()/m_ranksPerAxis[d],
                   "Lattice rank count overflow");
      expectedRanks *= m_ranksPerAxis[d];
    }
    POLY_VERIFY2(expectedRanks <= static_cast<std::size_t>(m_numPartitions),
                 "Product of lattice ranks per axis must not exceed the partition count");
  }

  std::vector<OwnerType>
  computeOwners(const std::vector<RealPoint>& globalPoints) const override {
    const auto& Q = Quantizer<Dimension>::instance();
    POLY_VERIFY2(Q.m_init, "The Quantizer must be initialized before lattice partitioning");
    std::vector<OwnerType> result(globalPoints.size());
    for (std::size_t i = 0; i < globalPoints.size(); ++i) {
      result[i] = owner(Q.quantize(globalPoints[i]), Q.minBound, Q.maxBound);
    }
    return result;
  }

  std::vector<OwnerType>
  computeOwners(const std::vector<QuantPoint>& globalPoints) const override {
    const auto& Q = Quantizer<Dimension>::instance();
    POLY_VERIFY2(Q.m_init, "The Quantizer must be initialized before lattice partitioning");
    std::vector<OwnerType> result(globalPoints.size());
    for (std::size_t i = 0; i < globalPoints.size(); ++i) {
      result[i] = owner(globalPoints[i], Q.minBound, Q.maxBound);
    }
    return result;
  }

  RanksPerAxis m_ranksPerAxis;
  using Partitioner<Dimension>::m_numPartitions;
private:
  std::size_t owner(const PointType& point,
                    const PointType& lower,
                    const PointType& upper) const {
    std::size_t rank = 0;
    std::size_t stride = 1;
    for (int d = 0; d < Dimension; ++d) {
      POLY_VERIFY2(point[d] >= lower[d] && point[d] <= upper[d],
                   "Generator is outside the lattice domain");

      const auto offset = static_cast<std::uint64_t>(point[d]) -
                          static_cast<std::uint64_t>(lower[d]);
      const auto extent = static_cast<std::uint64_t>(upper[d]) -
                          static_cast<std::uint64_t>(lower[d]);
      const auto tile = point[d] == upper[d] ?
        m_ranksPerAxis[d] - 1 :
        std::min(static_cast<unsigned>((static_cast<unsigned __int128>(offset)*m_ranksPerAxis[d])/extent),
                 m_ranksPerAxis[d] - 1);
      rank += stride*tile;
      stride *= m_ranksPerAxis[d];
    }
    return rank;
  }

  // Determine the optimal ranks per axis for a given number of partitions
  // without exceeding the number of partitions
  void optimalLattice2D() {
    unsigned bestx = 1;
    unsigned besty = m_numPartitions;
    for (unsigned x = 1; x <= std::sqrt(m_numPartitions); ++x) {
      unsigned y = m_numPartitions/x;
      unsigned bestTotal = bestx*besty;
      if (x*y > bestTotal ||
          (x*y == bestTotal && (y - x) < (besty - bestx))) {
        bestx = x;
        besty = y;
      }
    }
    m_ranksPerAxis = {bestx, besty};
  }

  void optimalLattice3D() {
    unsigned bestx = 1;
    unsigned besty = 1;
    unsigned bestz = m_numPartitions;
    for (unsigned x = 1; x*x*x <= m_numPartitions + x*x; ++x) {
      for (unsigned y = x; x*y*y <= m_numPartitions; ++y) {
        unsigned z = unsigned(m_numPartitions/(x*y));
        if (z < y) {
          continue;
        }
        unsigned total = x*y*z;
        unsigned bestTotal = bestx*besty*bestz;
        if (total > bestTotal ||
            (total == bestTotal && (z - x) < (bestz - bestx))) {
          bestx = x;
          besty = y;
          bestz = z;
        }
      }
    }
    m_ranksPerAxis = {bestx, besty, bestz};
  }
};

} // namespace polytope

#endif
