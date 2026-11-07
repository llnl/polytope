#include "DistributedTessellator.hh"

#include <algorithm>
#include <set>
#include <unordered_map>

#include "ParallelUtils.hh"

namespace polytope {

namespace {

// All ranks must agree before an empty input can skip the collective
// quantization and generator-exchange operations below.
inline
bool
hasAnyGenerators(const bool localHasGenerators) {
#ifdef POLYTOPE_ENABLE_MPI
  int localValue = localHasGenerators ? 1 : 0;
  int globalValue = 0;
  MPI_Allreduce(&localValue, &globalValue, 1, MPI_INT, MPI_MAX,
                Communicator::communicator());
  return globalValue != 0;
#else
  return localHasGenerators;
#endif
}

template<int Dimension, typename Generator>
void
appendVisibleGenerators(const QuantTessellation<Dimension>& visibleMesh,
                        const std::vector<Generator>& visibleRecords,
                        const std::set<int>& neighborRanks,
                        std::vector<std::vector<Generator>>& rankGenerators) {
  const int rank = Communicator::getRank();
  POLY_ASSERT2(visibleRecords.size() == visibleMesh.cellRank.size(),
               "Incorrect cell ranks in visible mesh");
  for (auto i = 0u; i < visibleRecords.size(); ++i) {
    const int source = visibleMesh.cellRank[i];
    if (source != rank && neighborRanks.find(source) == neighborRanks.end()) {
      POLY_ASSERT2(source >= 0 && source < int(rankGenerators.size()),
                   "Invalid visible generator rank");
      rankGenerators[source].push_back(visibleRecords[i]);
    }
  }
}

}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<int Dimension>
DistributedTessellator<Dimension>::
DistributedTessellator(Base& serialTessellator):
  Tessellator<Dimension, double>(serialTessellator),
  m_serialTessellator(serialTessellator) {
  // Make quantizer parallel
  m_Q.setParallel(true);
}

//------------------------------------------------------------------------------
// Synchronize generator exchange representation from the configured root rank.
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
synchronizeExchangePoints() {
  int value = m_exchangePoints ? 1 : 0;
  MPI_Bcast(&value, 1, MPI_INT,
            Communicator::getRoot(), Communicator::communicator());
  m_exchangePoints = value != 0;
}

//------------------------------------------------------------------------------
// Name
//------------------------------------------------------------------------------
template<int Dimension>
std::string
DistributedTessellator<Dimension>::
name() const {
  return "Distributed" + m_serialTessellator.name();
}

//------------------------------------------------------------------------------
// Unbounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
Tessellation<Dimension, double>
DistributedTessellator<Dimension>::
tessellate(const std::vector<Point<Dimension, RealType>>& points) {
  m_clipping = false;
  TessellationType mesh;
  if (!hasAnyGenerators(!points.empty())) return mesh;
  initQ(points);

  QuantizedTessellation quantmesh(m_Q, points);
  quantmesh.tessellatorName = this->name();
  this->tessellateQuantized(quantmesh);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  return mesh;
}

//------------------------------------------------------------------------------
// PLC-bounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
Tessellation<Dimension, double>
DistributedTessellator<Dimension>::
tessellate(const std::vector<Point<Dimension, RealType>>& points,
           const std::vector<Point<Dimension, RealType>>& PLCpoints,
           const PLC<Dimension>& geometry) {
  m_clipping = true;
  TessellationType mesh;
  if (!hasAnyGenerators(!points.empty())) return mesh;
  initQ(points, PLCpoints);

  QuantizedTessellation quantmesh(m_Q, points);
  quantmesh.tessellatorName = this->name();
  m_QPLCpoints = QuantPLC<Dimension>::quantizePoints(m_Q, PLCpoints);
  m_QPLC.emplace(m_Q, m_QPLCpoints, geometry);
  quantmesh.cullExternalPoints(*m_QPLC, m_QPLCpoints);
  this->tessellateQuantized(quantmesh);
  quantmesh.clipTessellation(*m_QPLC, m_QPLCpoints, m_serialTessellator);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  return mesh;
}

//------------------------------------------------------------------------------
// Partitioned unbounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
Tessellation<Dimension, double>
DistributedTessellator<Dimension>::
partitionAndTessellate(const std::vector<Point<Dimension, RealType>>& points,
                       const Partitioner<Dimension>& partitioner) {
  m_clipping = false;
  TessellationType mesh;
  if (!hasAnyGenerators(!points.empty())) return mesh;
  initQ(points);
  POLY_VERIFY2(partitioner.getNumPartitions() <= Communicator::getNRanks(),
               "Distributed partition count must not exceed the MPI rank count");

  const auto localPoints = partitioner.partition(points);
  QuantizedTessellation quantmesh(m_Q, localPoints);
  quantmesh.tessellatorName = this->name();
  this->tessellateQuantized(quantmesh);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  return mesh;
}

//------------------------------------------------------------------------------
// Partitioned PLC-bounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
Tessellation<Dimension, double>
DistributedTessellator<Dimension>::
partitionAndTessellate(const std::vector<Point<Dimension, RealType>>& points,
                       const std::vector<Point<Dimension, RealType>>& PLCpoints,
                       const PLC<Dimension>& geometry,
                       const Partitioner<Dimension>& partitioner) {
  m_clipping = true;
  TessellationType mesh;
  if (!hasAnyGenerators(!points.empty())) return mesh;
  initQ(points, PLCpoints);
  POLY_VERIFY2(partitioner.getNumPartitions() <= Communicator::getNRanks(),
               "Distributed partition count must not exceed the MPI rank count");

  m_QPLCpoints = QuantPLC<Dimension>::quantizePoints(m_Q, PLCpoints);
  m_QPLC.emplace(m_Q, m_QPLCpoints, geometry);
  const auto localPoints = partitioner.partition(points);
  QuantizedTessellation quantmesh(m_Q, localPoints);
  quantmesh.tessellatorName = this->name();
  quantmesh.cullExternalPoints(*m_QPLC, m_QPLCpoints);
  this->tessellateQuantized(quantmesh);
  quantmesh.clipTessellation(*m_QPLC, m_QPLCpoints, m_serialTessellator);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  return mesh;
}

//------------------------------------------------------------------------------
// Quantized distributed tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
tessellateQuantized(QuantizedTessellation& quantmesh) {
  int rank = Communicator::getRank();
  int nranks = Communicator::getNRanks();

  // Get the local visible generators and gather them on all processors
  auto visibleMesh = generateVisibleMesh(quantmesh);
  std::set<int> neighborRanks;
  if (visibleMesh.points.size() > 1) {
    neighborRanks = visibleMesh.neighboringRanks();
    // If clipping, holes can modify which ranks are neighbors
    if (m_clipping) {
      visibleMesh.clipTessellation(*m_QPLC, m_QPLCpoints, m_serialTessellator);
      std::set<int> clipped_neighbors = visibleMesh.neighboringRanks();
      neighborRanks.insert(clipped_neighbors.begin(), clipped_neighbors.end());
    }
  }

  // Check if any there is any overlap with other convex hulls
  if (!quantmesh.convexHull.m_convex) {
    quantmesh.makeConvexHull();
  }
  const auto localHull = quantmesh.convexHull;
  const bool validLocalHull = localHull.isValid(quantmesh.points) &&
                              localHull.m_convex &&
                              localHull.facets.size() >= Dimension + 1;
  const auto allHulls = allGatherHulls(validLocalHull, localHull, quantmesh.points);
  POLY_ASSERT2(int(allHulls.size()) == nranks,
               "Incorrect number of hulls after communication");
  if (validLocalHull) {
    for (int source = 0; source < nranks; ++source) {
      if (source == rank || !allHulls[source].first) continue;
      const auto& remoteHull = allHulls[source].second.first;
      bool boundsOverlap = true;
      for (int axis = 0; axis < Dimension; ++axis) {
        if (localHull.m_hiBounds[axis] < remoteHull.m_loBounds[axis] ||
            remoteHull.m_hiBounds[axis] < localHull.m_loBounds[axis]) {
          boundsOverlap = false;
          break;
        }
      }
      if (boundsOverlap &&
          QuantPLC<Dimension>::convexPLCIntersection(localHull, quantmesh.points,
                                                    remoteHull,
                                                    allHulls[source].second.second)) {
        neighborRanks.insert(source);
      }
    }
  }

  // Exchange complete neighbor generator sets using the selected representation.
  // Clipping also needs visible generators from non-neighbor ranks (see
  // DistributedHoleTests.cc). Add them before extending the mesh.
  bool addedGenerators = false;
  if (m_exchangePoints) {
    auto neighborGenerators =
      exchangeNeighborGenerators<Dimension, QuantizedPoint<Dimension>>(quantmesh.points, neighborRanks);
    if (quantmesh.points.empty()) {
      for (auto& generators : neighborGenerators) generators.clear();
    }
    if (m_clipping) {
      appendVisibleGenerators(visibleMesh, visibleMesh.points, neighborRanks,
                              neighborGenerators);
    }
    addedGenerators = std::any_of(neighborGenerators.begin(), neighborGenerators.end(),
                                  [](const auto& generators) { return !generators.empty(); });
    if (addedGenerators) {
      quantmesh.extend(neighborGenerators);
    }
  } else {
    auto neighborGenerators =
      exchangeNeighborGenerators<Dimension, QuantizedKey<Dimension>>(quantmesh.hashes, neighborRanks);
    if (quantmesh.points.empty()) {
      for (auto& generators : neighborGenerators) generators.clear();
    }
    if (m_clipping) {
      appendVisibleGenerators(visibleMesh, visibleMesh.hashes, neighborRanks,
                              neighborGenerators);
    }
    addedGenerators = std::any_of(neighborGenerators.begin(), neighborGenerators.end(),
                                  [](const auto& generators) { return !generators.empty(); });
    if (addedGenerators) {
      quantmesh.extend(neighborGenerators);
    }
  }
  if (addedGenerators) {
    m_serialTessellator.tessellateQuantized(quantmesh);
  } else if (quantmesh.cellRank.size() != quantmesh.points.size()) {
    // extend() normally initializes the local ranks for later filtering.
    quantmesh.cellRank.assign(quantmesh.points.size(), rank);
  }
}

//------------------------------------------------------------------------------
// Generate the visible mesh
//------------------------------------------------------------------------------
template<int Dimension>
QuantTessellation<Dimension>
DistributedTessellator<Dimension>::
generateVisibleMesh(QuantizedTessellation& quantmesh) {
  if (!quantmesh.points.empty()) {
    m_serialTessellator.tessellateQuantized(quantmesh);
  }
  QuantizedTessellation visibleMesh(quantmesh.quantizer());
  visibleMesh.tessellatorName = this->name();
  if (m_exchangePoints) {
    const auto visiblePoints = quantmesh.visibleGeneratorPoints();
    auto allVisibleRecords =
      allGatherGenerators<Dimension, QuantizedPoint<Dimension>>(visiblePoints);
    size_t nvisible = 0;
    for (const auto& rankPoints : allVisibleRecords) {
      nvisible += rankPoints.size();
    }
    if (nvisible > 0) {
      visibleMesh.extend(allVisibleRecords);
      m_serialTessellator.tessellateQuantized(visibleMesh);
    }
  } else {
    const auto visibleHashes = quantmesh.visibleGeneratorKeys();
    auto allVisibleRecords =
      allGatherGenerators<Dimension, QuantizedKey<Dimension>>(visibleHashes);
    size_t nvisible = 0;
    for (const auto& rankHashes : allVisibleRecords) {
      nvisible += rankHashes.size();
    }
    if (nvisible > 0) {
      visibleMesh.extend(allVisibleRecords);
      m_serialTessellator.tessellateQuantized(visibleMesh);
    }
  }
  return visibleMesh;
}

//------------------------------------------------------------------------------
// Add visible generator points to quantized mesh
// TODO: This might be overkill. It might be sufficient to only use the visible
// generators from neighbors of neighbors
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
addVisibleMesh(const QuantizedTessellation& visibleMesh,
               const std::set<int>& neighborRanks,
               QuantizedTessellation& quantmesh) {
  int nranks = Communicator::getNRanks();
  if (m_exchangePoints) {
    std::vector<std::vector<QuantizedPoint<Dimension>>> visibleGenerators(nranks);
    appendVisibleGenerators(visibleMesh, visibleMesh.points, neighborRanks,
                            visibleGenerators);
    quantmesh.extend(visibleGenerators);
  } else {
    std::vector<std::vector<QuantizedKey<Dimension>>> visibleGenerators(nranks);
    appendVisibleGenerators(visibleMesh, visibleMesh.hashes, neighborRanks,
                            visibleGenerators);
    quantmesh.extend(visibleGenerators);
  }
}

template class DistributedTessellator<2>;
template class DistributedTessellator<3>;

} // namespace polytope
