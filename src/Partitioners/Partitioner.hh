//----------------------------------------------------------------------------//
// Partitioner
//
// Common interface for MPI-domain partitioners.
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

  virtual ~Partitioner() = default;

  virtual std::string name() const = 0;

  //! Collectively produce this rank's local generators.  Implementations may
  //! accept replicated input or arbitrary rank-local input; see their class
  //! documentation for the corresponding input contract.
  virtual std::vector<RealPoint>
  partition(const std::vector<RealPoint>& localPoints) const = 0;

  //! Number of logical output partitions.
  virtual unsigned getNumPartitions() const = 0;
};

//----------------------------------------------------------------------------//
// ReplicatedPartitioner
//
// Base for partitioners that calculate owners from a generator list replicated
// on every rank.  It retains the owner-query API used by serial callers.
//----------------------------------------------------------------------------//
template<int Dimension>
class ReplicatedPartitioner: public Partitioner<Dimension> {
public:
  using OwnerType = typename Partitioner<Dimension>::OwnerType;
  using RealPoint = typename Partitioner<Dimension>::RealPoint;
  using QuantPoint = typename Partitioner<Dimension>::QuantPoint;

  //! A partitioning together with the owner of each input generator.
  template<typename CoordType>
  struct PartitionResult {
    std::vector<std::vector<Point<Dimension, CoordType>>> generatorsByPartition;
    std::vector<OwnerType> ownerByInput;
  };

  //! The number of logical partitions.  These are MPI ranks only when used
  //! by DistributedTessellator.
  ReplicatedPartitioner(const unsigned numPartitions = Communicator::getNRanks()) :
    m_numPartitions(numPartitions) {
    POLY_VERIFY(m_numPartitions > 0);
  }

  void setNumPartitions(const unsigned numPartitions) {
    m_numPartitions = numPartitions;
    POLY_VERIFY(m_numPartitions > 0);
  }

  unsigned getNumPartitions() const override { return m_numPartitions; }

  //! Return this rank's subset of an identically replicated generator list.
  std::vector<RealPoint>
  partition(const std::vector<RealPoint>& globalPoints) const override {
    return computeLocalPartition(globalPoints);
  }

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

  unsigned m_numPartitions;
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
};

} // namespace polytope

#endif
