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

#include "Communicator.hh"
#include "MortonKeyTraits.hh"
#include "Point.hh"
#include "polytope_internal.hh"

namespace polytope {

template<int Dimension>
class Partitioner {
public:
  using PointType = QuantizedPoint<Dimension>;

  virtual ~Partitioner() = default;

  //! Return the generators owned by this MPI rank.  Inputs must have been
  //! quantized with a globally consistent Quantizer instance.
  virtual std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const = 0;
};

//----------------------------------------------------------------------------//
// RandomPartitioner
//
// Ownership is a deterministic hash of the seed, quantized coordinates, and
// input ordinal.  Every rank must receive the same, identically ordered input.
//----------------------------------------------------------------------------//
template<int Dimension>
class RandomPartitioner: public Partitioner<Dimension> {
public:
  using PointType = typename Partitioner<Dimension>::PointType;

  explicit RandomPartitioner(const std::uint64_t seed):
    m_seed(seed) {
  }

  std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const override {
    const auto rank = static_cast<std::uint64_t>(Communicator::getRank());
    const auto nranks = static_cast<std::uint64_t>(Communicator::getNProcs());
    POLY_VERIFY(nranks > 0);

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

private:
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
};

//----------------------------------------------------------------------------//
// LatticePartitioner
//----------------------------------------------------------------------------//
template<int Dimension>
class LatticePartitioner: public Partitioner<Dimension> {
public:
  using PointType = typename Partitioner<Dimension>::PointType;
  using RanksPerAxis = std::array<unsigned, Dimension>;

  LatticePartitioner(const PointType& lower,
                     const PointType& upper,
                     const RanksPerAxis& ranksPerAxis):
    m_lower(lower),
    m_upper(upper),
    m_ranksPerAxis(ranksPerAxis) {
    std::size_t expectedRanks = 1;
    for (int d = 0; d < Dimension; ++d) {
      POLY_VERIFY2(m_lower[d] < m_upper[d],
                   "Each lattice lower bound must be less than its upper bound");
      POLY_VERIFY2(m_ranksPerAxis[d] > 0,
                   "Each lattice axis must have at least one rank");
      POLY_VERIFY2(expectedRanks <= std::numeric_limits<std::size_t>::max()/m_ranksPerAxis[d],
                   "Lattice rank count overflow");
      expectedRanks *= m_ranksPerAxis[d];
    }
    POLY_VERIFY2(expectedRanks == static_cast<std::size_t>(Communicator::getNProcs()),
                 "Product of lattice ranks per axis must equal the MPI rank count");
  }

  std::vector<PointType>
  computePartition(const std::vector<PointType>& globalPoints) const override {
    const auto rank = static_cast<std::size_t>(Communicator::getRank());
    std::vector<PointType> result;
    for (const auto& point : globalPoints) {
      if (owner(point) == rank) {
        result.push_back(point);
      }
    }
    return result;
  }

private:
  std::size_t owner(const PointType& point) const {
    std::size_t rank = 0;
    std::size_t stride = 1;
    for (int d = 0; d < Dimension; ++d) {
      POLY_VERIFY2(point[d] >= m_lower[d] && point[d] <= m_upper[d],
                   "Generator is outside the lattice domain");

      const auto offset = static_cast<std::uint64_t>(point[d]) -
                          static_cast<std::uint64_t>(m_lower[d]);
      const auto extent = static_cast<std::uint64_t>(m_upper[d]) -
                          static_cast<std::uint64_t>(m_lower[d]);
      const auto tile = point[d] == m_upper[d] ?
        m_ranksPerAxis[d] - 1 :
        std::min(static_cast<unsigned>((static_cast<unsigned __int128>(offset)*m_ranksPerAxis[d])/extent),
                 m_ranksPerAxis[d] - 1);
      rank += stride*tile;
      stride *= m_ranksPerAxis[d];
    }
    return rank;
  }

  PointType m_lower, m_upper;
  RanksPerAxis m_ranksPerAxis;
};

} // namespace polytope

#endif
