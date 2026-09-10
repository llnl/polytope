//----------------------------------------------------------------------------//
// QuasiVoronoiPartitioner
//
// Assigns a random point to each rank and gathers spatially nearby points to
// that rank. It also uses Lloyd's algorithm to improve load balancing.
//----------------------------------------------------------------------------//

#ifndef __Polytope_QuasiVoronoiPartitioner__
#define __Polytope_QuasiVoronoiPartitioner__

#include "Partitioner.hh"

namespace polytope {

template<int Dimension>
class QuasiVoronoiPartitioner: public Partitioner<Dimension> {
public:
  using OwnerType = typename Partitioner<Dimension>::OwnerType;
  using RealPoint = typename Partitioner<Dimension>::RealPoint;
  using QuantPoint = typename Partitioner<Dimension>::QuantPoint;
  using Partitioner<Dimension>::computeOwners;

  explicit QuasiVoronoiPartitioner(const unsigned seed,
                                   const unsigned numPartitions = Communicator::getNRanks(),
                                   const unsigned Niter = 100):
    Partitioner<Dimension>(numPartitions),
    m_seed(seed),
    m_Niter(Niter) { }

  virtual std::string name() const override { return "QuasiVoronoiPartitioner"; }

  std::vector<OwnerType>
  computeOwners(const std::vector<RealPoint>& globalPoints) const override {
    return computeOwnersImpl(globalPoints);
  }

  std::vector<OwnerType>
  computeOwners(const std::vector<QuantPoint>& globalPoints) const override {
    return computeOwnersImpl(globalPoints);
  }

  void setNumIter(const unsigned Niter) { m_Niter = Niter; }
  unsigned getNumIter() { return m_Niter; }

  unsigned m_seed;
  using Partitioner<Dimension>::m_numPartitions;
  unsigned m_Niter; // Number of Lloyd's iterations to run

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
    // Set of which generator point indices are assigned to a rank
    std::set<unsigned> procPointIndices;
    std::vector<Point<Dimension, CoordType>> rankOrigins;
    rankOrigins.reserve(nranks);
    // Assign each rank a random generator point
    for (int rank = 0; rank < nranks; ++rank) {
      auto i = distrib(gen);
      // Ensure generator point is not assigned to another rank
      while (procPointIndices.find(i) != procPointIndices.end()) {
        i = distrib(gen);
      }
      procPointIndices.insert(i);
      rankOrigins.push_back(globalPoints[i]);
    }

    for (int i = 0; i < m_Niter; ++i) {
      updateRankSites(globalPoints, rankOrigins);
    }

    // Iterate over each point and determine which proc seed is closest
    for (std::size_t i = 0; i < N; ++i) {
      const auto& point = globalPoints[i];
      auto diff = point - rankOrigins[0];
      auto minDist = magnitude2(diff);
      int proc_owner = 0;
      for (int ip = 1; ip < nranks; ++ip) {
        diff = point - rankOrigins[ip];
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

  // Update the rank sites based the centroid of its point cloud
  template<typename CoordType>
  void
  updateRankSites(const std::vector<Point<Dimension, CoordType>>& globalPoints,
                  std::vector<Point<Dimension, CoordType>>& rankOrigins) const {
    using Wide = typename WideIntHelper<Dimension, CoordType>::type;
    const auto nranks = m_numPartitions;
    const auto N = globalPoints.size();
    std::vector<std::vector<Point<Dimension, CoordType>>> procPoints(nranks);
    std::vector<std::vector<unsigned>> rankGens(nranks);
    // Gather generators closest to each rank seed site
    for (std::size_t i = 0; i < N; ++i) {
      const auto& point = globalPoints[i];
      auto diff = point - rankOrigins[0];
      auto minDist = magnitude2(diff);
      int proc_owner = 0;
      for (int ip = 1; ip < nranks; ++ip) {
        diff = point - rankOrigins[ip];
        auto dist = magnitude2(diff);
        if (dist < minDist) {
          proc_owner = ip;
          minDist = dist;
        }
      }
      procPoints[proc_owner].push_back(point);
      rankGens[proc_owner].push_back(i);
    }
    for (std::size_t i = 0; i < nranks; ++i) {
      rankOrigins[i] = pointCentroid(procPoints[i]);
    }
  }
};

} // namespace polytope

#endif
