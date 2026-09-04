//----------------------------------------------------------------------------//
// Partitioner
//
// Deterministic MPI-domain partitioners for replicated generator point sets.
//----------------------------------------------------------------------------//
#ifndef __Polytope_Partitioner__
#define __Polytope_Partitioner__

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <type_traits>
#include <vector>
#include <random>

#include "Communicator.hh"
#include "Point.hh"
#include "Quantizer.hh"
#include "polytope_internal.hh"
#include "GeomUtils.hh"
#include "QuantTessellation.hh"

namespace polytope {

template<int Dimension>
class Partitioner {
public:
  using OwnerType = unsigned;
  using RealPoint = Point<Dimension, double>;
  using QuantPoint = QuantizedPoint<Dimension>;

  //! A partitioning together with the owner of each input generator.
  template<typename CoordType>
  struct PartitionResult {
    std::vector<std::vector<Point<Dimension, CoordType>>> generatorsByPartition;
    std::vector<OwnerType> ownerByInput;
  };

  //! The number of logical partitions.  These are MPI ranks only when used
  //! by DistributedTessellator.
  Partitioner(const unsigned numPartitions = Communicator::getNProcs()) :
    m_numPartitions(numPartitions) {
    POLY_VERIFY(m_numPartitions > 0);
  }

  virtual ~Partitioner() = default;

  virtual std::string name() const = 0;

  void setNumPartitions(const unsigned numPartitions) {
    m_numPartitions = numPartitions;
    POLY_VERIFY(m_numPartitions > 0);
  }

  unsigned numPartitions() const { return m_numPartitions; }

  //! Determine a logical owner for every real-valued input generator.  The
  //! returned vector is in the same order as globalPoints.
  virtual std::vector<OwnerType>
  computeOwners(const std::vector<RealPoint>& globalPoints) const = 0;

  //! Determine a logical owner for every quantized input generator.  This is
  //! kept separate from the real-valued overload so partitioners that use
  //! integer coordinates can preserve them exactly.
  virtual std::vector<OwnerType>
  computeOwners(const std::vector<QuantPoint>& globalPoints) const = 0;

  //! Determine owners for another coordinate representation by treating it
  //! as real-valued. The real and quantized overloads above avoid conversion.
  template<typename CoordType,
           typename std::enable_if_t<!std::is_same_v<CoordType, double> &&
                                     !std::is_same_v<CoordType, QuantizedCoordinate<Dimension>>, int> = 0>
  std::vector<OwnerType>
  computeOwners(const std::vector<Point<Dimension, CoordType>>& globalPoints) const {
    return computeOwners(realPoints(globalPoints));
  }

  //! Return the complete partitioning and the owner of every input point.
  template<typename CoordType>
  PartitionResult<CoordType>
  computePartitionResult(const std::vector<Point<Dimension, CoordType>>& globalPoints) const {
    // To ensure consistency, hash and sort the points
    QuantTessellation<Dimension> qmesh(globalPoints);
    std::vector<Point<Dimension, CoordType>> sortedPoints;
    if constexpr (std::is_same_v<CoordType, double>) {
      sortedPoints = qmesh.getRealPoints();
    } else {
      sortedPoints = qmesh.getQuantizedPoints();
    }
    PartitionResult<CoordType> result;
    result.ownerByInput = computeOwners(sortedPoints);
    POLY_VERIFY(result.ownerByInput.size() == sortedPoints.size());
    result.generatorsByPartition.resize(m_numPartitions);
    for (std::size_t i = 0; i < sortedPoints.size(); ++i) {
      const auto owner = result.ownerByInput[i];
      POLY_VERIFY(owner < m_numPartitions);
      result.generatorsByPartition[owner].push_back(sortedPoints[i]);
    }
    return result;
  }

  //! Return a vector of vector of generators arranged by result[rank][point].
  template<typename CoordType>
  std::vector<std::vector<Point<Dimension, CoordType>>>
  computePartition(const std::vector<Point<Dimension, CoordType>>& globalPoints) const {
    return computePartitionResult(globalPoints).generatorsByPartition;
  }

  //! Label every generator in a mesh with its logical partition. This supports
  //! serial tessellation followed by output through mesh.cellRank.
  template<typename CoordType>
  void assignCellRanks(Tessellation<Dimension, CoordType>& mesh) const {
    const auto owners = computeOwners(mesh.points);
    POLY_VERIFY(owners.size() == mesh.points.size());
    mesh.cellRank.resize(owners.size());
    for (std::size_t i = 0; i < owners.size(); ++i) {
      mesh.cellRank[i] = static_cast<int>(owners[i]);
    }
  }

  //! As above, using a flat real-valued generator array.
  void assignCellRanks(Tessellation<Dimension, double>& mesh,
                       const std::vector<double>& flatRealPoints) const {
    assignCellRanks(mesh, extractCoords<Dimension, double>(flatRealPoints));
  }

  //! Return the generators owned by this MPI rank.
  template<typename CoordType>
  std::vector<Point<Dimension, CoordType>>
  computeLocalPartition(const std::vector<Point<Dimension, CoordType>>& globalPoints) const {
    const auto rank = Communicator::getRank();
    auto result = computePartition(globalPoints);
    if (rank < result.size()) {
      return result[rank];
    } else {
      return std::vector<Point<Dimension, CoordType>>();
    }
  }

  //! Return the generators owned by this MPI rank.
  template<typename CoordType>
  std::vector<Point<Dimension, CoordType>>
  computeLocalPartition(const std::vector<CoordType>& flatRealPoints) const {
    auto globalPoints = extractCoords<Dimension, CoordType>(flatRealPoints);
    return computeLocalPartition(globalPoints);
  }

protected:
  template<typename CoordType>
  static std::vector<RealPoint>
  realPoints(const std::vector<Point<Dimension, CoordType>>& points) {
    std::vector<RealPoint> result;
    result.reserve(points.size());
    if constexpr (std::is_same_v<CoordType, QuantizedCoordinate<Dimension>>) {
      const auto& Q = Quantizer<Dimension>::instance();
      for (const auto& point : points) result.push_back(Q.dequantize(point));
    } else {
      for (const auto& point : points) result.push_back(point.template type_cast<double>());
    }
    return result;
  }

  unsigned m_numPartitions;
};

//----------------------------------------------------------------------------//
// RandomPartitioner
//
// Ownership is a deterministic hash of the seed, generator coordinates, and
// input ordinal. Every rank must receive the same, identically ordered input.
//----------------------------------------------------------------------------//
template<int Dimension>
class RandomPartitioner: public Partitioner<Dimension> {
public:
  using OwnerType = typename Partitioner<Dimension>::OwnerType;
  using RealPoint = typename Partitioner<Dimension>::RealPoint;
  using QuantPoint = typename Partitioner<Dimension>::QuantPoint;
  using Partitioner<Dimension>::computeOwners;

  explicit RandomPartitioner(const std::uint64_t seed,
                             const unsigned numPartitions = Communicator::getNProcs()):
    Partitioner<Dimension>(numPartitions),
    m_seed(seed){ }

  virtual std::string name() const override { return "RandomPartitioner"; }

  std::vector<OwnerType>
  computeOwners(const std::vector<RealPoint>& globalPoints) const override {
    return computeOwnersImpl(globalPoints);
  }

  std::vector<OwnerType>
  computeOwners(const std::vector<QuantPoint>& globalPoints) const override {
    return computeOwnersImpl(globalPoints);
  }

protected:
  template<typename CoordType>
  std::vector<OwnerType>
  computeOwnersImpl(const std::vector<Point<Dimension, CoordType>>& globalPoints) const {
    const auto nranks = m_numPartitions;
    std::vector<OwnerType> result(globalPoints.size());
    for (std::size_t i = 0; i < globalPoints.size(); ++i) {
      const auto& point = globalPoints[i];
      result[i] = owner(point, i, nranks);
    }
    return result;
  }

  static std::uint64_t mix(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }

  template<typename CoordType>
  std::uint64_t owner(const Point<Dimension, CoordType>& point,
                      const std::size_t ordinal,
                      const std::uint64_t nranks) const {
    auto hash = mix(m_seed);
    for (int d = 0; d < Dimension; ++d) {
      hash = mix(hash ^ mix(coordinateHash(point[d])));
    }
    return mix(hash ^ mix(static_cast<std::uint64_t>(ordinal))) % nranks;
  }

  template<typename CoordType>
  static std::uint64_t coordinateHash(const CoordType coordinate) {
    if constexpr (std::is_floating_point_v<CoordType>) {
      static_assert(sizeof(CoordType) <= sizeof(std::uint64_t));
      std::uint64_t result = 0;
      std::memcpy(&result, &coordinate, sizeof(CoordType));
      return result;
    } else {
      return static_cast<std::uint64_t>(coordinate);
    }
  }

  std::uint64_t m_seed;
  using Partitioner<Dimension>::m_numPartitions;
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
  using OwnerType = typename Partitioner<Dimension>::OwnerType;
  using RealPoint = typename Partitioner<Dimension>::RealPoint;
  using QuantPoint = typename Partitioner<Dimension>::QuantPoint;
  using Partitioner<Dimension>::computeOwners;

  explicit QuasiVoronoiPartitioner(const unsigned seed,
                                   const unsigned numPartitions = Communicator::getNProcs()):
    Partitioner<Dimension>(numPartitions),
    m_seed(seed) { }

  virtual std::string name() const override { return "QuasiVoronoiPartitioner"; }

  std::vector<OwnerType>
  computeOwners(const std::vector<RealPoint>& globalPoints) const override {
    return computeOwnersImpl(globalPoints);
  }

  std::vector<OwnerType>
  computeOwners(const std::vector<QuantPoint>& globalPoints) const override {
    return computeOwnersImpl(globalPoints);
  }

private:
  template<typename CoordType>
  std::vector<OwnerType>
  computeOwnersImpl(const std::vector<Point<Dimension, CoordType>>& globalPoints) const {
    if (globalPoints.size() == 0) {
      return std::vector<OwnerType>();
    }
    const auto nranks = m_numPartitions;
    std::mt19937 gen(m_seed);
    const auto N = static_cast<unsigned>(globalPoints.size());
    std::uniform_int_distribution<unsigned> distrib(0, N-1);
    const int Nmin = std::min(nranks, N);
    std::vector<OwnerType> result(N);
    if (N <= nranks) {
      for (int i = 0; i < Nmin; ++i) {
        result[i] = i;
      }
      return result;
    }
    std::set<unsigned> procPointIndices;
    std::vector<Point<Dimension, CoordType>> procPoints;
    procPoints.reserve(nranks);
    // Assign each rank a random generator point
    for (int rank = 0; rank < nranks; ++rank) {
      auto i = distrib(gen);
      // Ensure generator point is not assigned to another rank
      while (procPointIndices.find(i) != procPointIndices.end()) {
        i = distrib(gen);
      }
      procPointIndices.insert(i);
      procPoints.push_back(globalPoints[i]);
    }

    // Iterate over each point and determine which proc seed is closest
    for (std::size_t i = 0; i < N; ++i) {
      const auto& point = globalPoints[i];
      auto diff = point - procPoints[0];
      auto minDist = magnitude2(diff);
      int proc_owner = 0;
      for (int ip = 1; ip < nranks; ++ip) {
        diff = point - procPoints[ip];
        auto dist = magnitude2(diff);
        if (dist < minDist) {
          proc_owner = ip;
          minDist = dist;
        }
      }
      result[i] = proc_owner;
    }
    return result;
  }
  unsigned m_seed;
  using Partitioner<Dimension>::m_numPartitions;
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
  using PointType = QuantizedPoint<Dimension>;
  using OwnerType = typename Partitioner<Dimension>::OwnerType;
  using RealPoint = typename Partitioner<Dimension>::RealPoint;
  using QuantPoint = typename Partitioner<Dimension>::QuantPoint;
  using Partitioner<Dimension>::computeOwners;
  using RanksPerAxis = std::array<unsigned, Dimension>;

  virtual std::string name() const override { return "LatticePartitioner"; }

  explicit LatticePartitioner(const RanksPerAxis& ranksPerAxis,
                              const unsigned numPartitions = Communicator::getNProcs()):
    Partitioner<Dimension>(numPartitions),
    m_ranksPerAxis(ranksPerAxis) {
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
protected:
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

  RanksPerAxis m_ranksPerAxis;
};

} // namespace polytope

#endif
