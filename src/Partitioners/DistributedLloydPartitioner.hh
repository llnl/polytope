//----------------------------------------------------------------------------//
// DistributedLloydPartitioner
//
// Collectively redistributes rank-local generators.  A lattice decomposition
// supplies deterministic geometric Lloyd sites; subsequent Lloyd iterations
// improve spatial locality without replicating the generator set.
//----------------------------------------------------------------------------//
#ifndef __Polytope_DistributedLloydPartitioner__
#define __Polytope_DistributedLloydPartitioner__

#include <array>
#include <limits>
#include <string>
#include <vector>

#include "LatticePartitioner.hh"
#ifdef POLYTOPE_ENABLE_MPI
#include "ParallelUtils.hh"
#endif

namespace polytope {

template<int Dimension>
class DistributedLloydPartitioner {
public:
  using RealPoint = Point<Dimension, double>;

  explicit DistributedLloydPartitioner(const unsigned Niter = 100):
    m_Niter(Niter) { }

  std::string name() const { return "DistributedLloydPartitioner"; }

  void setNumIter(const unsigned Niter) { m_Niter = Niter; }
  unsigned getNumIter() const { return m_Niter; }

  //! Collectively redistribute arbitrary rank-local input.  Point coordinates
  //! and indices are preserved exactly during every exchange.
  std::vector<RealPoint>
  partition(const std::vector<RealPoint>& localPoints) const {
    auto& Q = Quantizer<Dimension>::instance();
    if (!Q.m_init && !localPoints.empty()) Q.init(localPoints);
#ifdef POLYTOPE_ENABLE_MPI
    const auto nranks = Communicator::getNRanks();
    LatticePartitioner<Dimension> lattice(static_cast<unsigned>(nranks));
    auto points = redistribute(localPoints, lattice.computeOwners(localPoints));

    // Start with the geometric centers of the lattice tiles.  This leaves the
    // first assignment equal to the lattice split and makes lattice-then-Lloyd
    // explicit even for initially empty tiles.
    auto sites = latticeSites(lattice);
    for (unsigned iteration = 0; iteration < m_Niter; ++iteration) {
      std::vector<unsigned> owners(points.size());
      for (std::size_t i = 0; i < points.size(); ++i) {
        owners[i] = nearestSite(points[i], sites);
      }
      points = redistribute(points, owners);

      auto localSite = sites[Communicator::getRank()];
      if (!points.empty()) {
        localSite = centroid(points);
      }
      sites = gatherSites(localSite);
    }
    return points;
#else
    // In a non-MPI build there is one logical partition, so redistribution is
    // a no-op.  Initialize the Quantizer consistently with other serial APIs.
    return localPoints;
#endif
  }

private:
#ifdef POLYTOPE_ENABLE_MPI
  std::vector<RealPoint>
  redistribute(const std::vector<RealPoint>& points,
               const std::vector<unsigned>& owners) const {
    const auto nranks = Communicator::getNRanks();
    POLY_VERIFY(points.size() == owners.size());
    std::vector<std::vector<RealPoint>> byDestination(nranks);
    for (std::size_t i = 0; i < points.size(); ++i) {
      POLY_VERIFY2(owners[i] < static_cast<unsigned>(nranks),
                   "Generator owner is outside the MPI communicator");
      byDestination[owners[i]].push_back(points[i]);
    }
    return redistributeGenerators<Dimension>(byDestination);
  }

  std::vector<RealPoint>
  latticeSites(const LatticePartitioner<Dimension>& lattice) const {
    const auto& Q = Quantizer<Dimension>::instance();
    const auto lower = Q.dequantize(Q.minBound);
    const auto upper = Q.dequantize(Q.maxBound);
    const auto nranks = lattice.m_numPartitions;
    std::vector<RealPoint> result(nranks);
    for (int rank = 0; rank < nranks; ++rank) {
      auto tileIndex = static_cast<unsigned>(rank);
      RealPoint site;
      for (int d = 0; d < Dimension; ++d) {
        const auto tile = tileIndex % lattice.m_ranksPerAxis[d];
        tileIndex /= lattice.m_ranksPerAxis[d];
        site[d] = lower[d] + (upper[d] - lower[d]) *
          (static_cast<double>(tile) + 0.5) / lattice.m_ranksPerAxis[d];
      }
      result[rank] = site;
    }
    return result;
  }

  unsigned
  nearestSite(const RealPoint& point,
              const std::vector<RealPoint>& sites) const {
    POLY_VERIFY(!sites.empty());
    auto owner = 0u;
    auto minimumDistance = magnitude2(point - sites[0]);
    for (std::size_t rank = 1; rank < sites.size(); ++rank) {
      const auto distance = magnitude2(point - sites[rank]);
      // Strict comparison gives rank-order tie breaking.
      if (distance < minimumDistance) {
        owner = static_cast<unsigned>(rank);
        minimumDistance = distance;
      }
    }
    return owner;
  }

  RealPoint centroid(const std::vector<RealPoint>& points) const {
    POLY_VERIFY(!points.empty());
    RealPoint result;
    for (const auto& point : points) result += point;
    return result/static_cast<double>(points.size());
  }

  std::vector<RealPoint>
  gatherSites(const RealPoint& localSite) const {
    const auto gathered = allGatherGenerators<Dimension>(std::vector<RealPoint>{localSite});
    POLY_VERIFY(gathered.size() == static_cast<std::size_t>(Communicator::getNRanks()));
    std::vector<RealPoint> result(gathered.size());
    for (std::size_t rank = 0; rank < gathered.size(); ++rank) {
      POLY_VERIFY(gathered[rank].size() == 1);
      result[rank] = gathered[rank][0];
    }
    return result;
  }
#endif

  unsigned m_Niter;
};

} // namespace polytope

#endif
