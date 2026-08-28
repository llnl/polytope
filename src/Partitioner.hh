//----------------------------------------------------------------------------//
// Partitioner
//
// Deterministic MPI-domain partitioners for replicated quantized generator
// point sets.
//----------------------------------------------------------------------------//
#ifndef __Polytope_Partitioner__
#define __Polytope_Partitioner__

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>
#include <random>

#include "Communicator.hh"
#include "Point.hh"
#include "Quantizer.hh"
#include "polytope_internal.hh"
#include "GeomUtils.hh"

namespace polytope {

template<int Dimension>
class Partitioner {
public:
  using PointType = QuantizedPoint<Dimension>;

  virtual ~Partitioner() = default;

  virtual std::string name() const = 0;

  //! Return the generators owned by this MPI rank.  Inputs must have been
  //! quantized with a globally consistent Quantizer instance.
  virtual std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const = 0;
};

//----------------------------------------------------------------------------//
// RandomPartitioner
//
// Ownership is a deterministic hash of the seed, quantized coordinates, and
// input ordinal. Every rank must receive the same, identically ordered input.
//----------------------------------------------------------------------------//
template<int Dimension>
class RandomPartitioner: public Partitioner<Dimension> {
public:
  using PointType = typename Partitioner<Dimension>::PointType;

  explicit RandomPartitioner(const std::uint64_t seed,
                             const unsigned maxNRank = Communicator::getNProcs()):
    m_seed(seed),
    m_maxNRank(maxNRank) {
    POLY_VERIFY(m_maxNRank > 0 && m_maxNRank <= Communicator::getNProcs());
  }

  virtual std::string name() const override { return "RandomPartitioner"; }

  std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const override {
    const auto rank = static_cast<std::uint64_t>(Communicator::getRank());
    const auto nranks = m_maxNRank;

    std::vector<PointType> result;
    result.reserve(globalPoints.size()/nranks + 1);
    for (std::size_t i = 0; i < globalPoints.size(); ++i) {
      const auto& point = globalPoints[i];
      if (owner(point, i, nranks) == rank) {
        result.push_back(point);
      }
    }
    return result;
  }

protected:
  static std::uint64_t mix(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }

  std::uint64_t owner(const PointType& point,
                      const std::size_t ordinal,
                      const std::uint64_t nranks) const {
    auto hash = mix(m_seed);
    for (int d = 0; d < Dimension; ++d) {
      hash = mix(hash ^ mix(static_cast<std::uint64_t>(point[d])));
    }
    return mix(hash ^ mix(static_cast<std::uint64_t>(ordinal))) % nranks;
  }

  std::uint64_t m_seed;
  unsigned m_maxNRank; // Max number of ranks to use
};

//----------------------------------------------------------------------------//
// QuasiVoronoiPartitioner
//
// Assigns a random point to each rank and gathers spatially nearby points to
// that rank.
//----------------------------------------------------------------------------//
template<int Dimension>
class QuasiVoronoiPartitioner: public Partitioner<Dimension> {
public:
  using PointType = typename RandomPartitioner<Dimension>::PointType;

  explicit QuasiVoronoiPartitioner(const std::uint64_t seed,
                                   const unsigned maxNRank = Communicator::getNProcs()):
    m_seed(seed),
    m_maxNRank(maxNRank) {
  }

  virtual std::string name() const override { return "QuasiVoronoiPartitioner"; }

  std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const override {
    const auto rank = Communicator::getRank();
    const auto nranks = m_maxNRank;
    POLY_VERIFY(nranks > 0);
    std::mt19937 gen(m_seed);
    const auto N = globalPoints.size();
    std::uniform_int_distribution<unsigned> distrib(0, N);
    std::set<unsigned> procPointIndices;
    std::vector<PointType> procPoints;
    procPoints.reserve(nranks);
    // Assign each rank a random generator point
    for (int rank = 0; rank < nranks; ++rank) {
      auto i = distrib(gen);
      while (procPointIndices.find(i) != procPointIndices.end()) {
        i = distrib(gen);
      }
      procPointIndices.insert(i);
      procPoints.push_back(globalPoints[i]);
    }

    std::vector<PointType> result;
    result.reserve(globalPoints.size()/nranks + 1);
    // Iterate over each point and determine which proc seed it closest
    for (std::size_t i = 0; i < N; ++i) {
      int proc_owner = 0;
      const auto& point = globalPoints[i];
      auto diff = point - procPoints[0];
      auto minDist = qmagnitude2(diff);
      for (int ip = 1; ip < nranks; ++ip) {
        diff = point - procPoints[ip];
        auto dist = qmagnitude2(diff);
        if (dist < minDist) {
          proc_owner = ip;
          minDist = dist;
        }
      }
      if (proc_owner == rank) {
        result.push_back(point);
      }
    }
    return result;
  }
  std::uint64_t m_seed;
  unsigned m_maxNRank; // Max number of ranks to use
};

//----------------------------------------------------------------------------//
// LatticePartitioner
//
// The complete partitioning domain is Quantizer<Dimension>::minBound through
// Quantizer<Dimension>::maxBound, inclusive. The Quantizer must be initialized
// before computePartition is called. This domain is divided uniformly into
// ranksPerAxis[d] tiles along each axis; a point on the global upper bound
// belongs to the final tile on that axis.
//----------------------------------------------------------------------------//
template<int Dimension>
class LatticePartitioner: public Partitioner<Dimension> {
public:
  using PointType = typename Partitioner<Dimension>::PointType;
  using RanksPerAxis = std::array<unsigned, Dimension>;

  virtual std::string name() const override { return "LatticePartitioner"; }

  explicit LatticePartitioner(const RanksPerAxis& ranksPerAxis):
    m_ranksPerAxis(ranksPerAxis) {
    std::size_t expectedRanks = 1;
    for (int d = 0; d < Dimension; ++d) {
      POLY_VERIFY2(m_ranksPerAxis[d] > 0,
                   "Each lattice axis must have at least one rank");
      POLY_VERIFY2(expectedRanks <= std::numeric_limits<std::size_t>::max()/m_ranksPerAxis[d],
                   "Lattice rank count overflow");
      expectedRanks *= m_ranksPerAxis[d];
    }
    POLY_VERIFY2(expectedRanks <= static_cast<std::size_t>(Communicator::getNProcs()),
                 "Product of lattice ranks per axis must less than or equal the MPI rank count");
  }

  std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const override {
    const auto rank = static_cast<std::size_t>(Communicator::getRank());
    const auto& Q = Quantizer<Dimension>::instance();
    POLY_VERIFY2(Q.m_init, "The Quantizer must be initialized before lattice partitioning");
    std::vector<PointType> result;
    for (const auto& point : globalPoints) {
      if (owner(point, Q.minBound, Q.maxBound) == rank) {
        result.push_back(point);
      }
    }
    return result;
  }

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

  RanksPerAxis m_ranksPerAxis;
};

} // namespace polytope

#endif
