#include "DistributedTessellator.hh"

#include <algorithm>
#include <set>
#include <unordered_map>

#include "ParallelUtils.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<int Dimension>
DistributedTessellator<Dimension>::
DistributedTessellator(Base& serialTessellator):
  m_serialTessellator(serialTessellator),
  m_keyEncode(Quantizer<Dimension>::instance().keyEncoding()) {
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
  return "DistributedTessellator_" + m_serialTessellator.name();
}

//------------------------------------------------------------------------------
// Unbounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
tessellate(const std::vector<Point<Dimension, RealType>>& points,
           TessellationType& mesh) {
  POLY_ASSERT(mesh.empty());
  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    Q.init(points);
  }

  m_keyEncode = Q.keyEncoding();
  QuantizedTessellation quantmesh(points);
  this->tessellateQuantized(quantmesh);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  checkEncoding();
}

//------------------------------------------------------------------------------
// PLC-bounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
tessellate(const std::vector<Point<Dimension, RealType>>& points,
           const std::vector<RealType>& PLCpoints,
           const PLC<Dimension>& geometry,
           TessellationType& mesh) {
  m_clipping = true;
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(PLCpoints.size() % Dimension == 0);
  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    if (PLCpoints.empty()) {
      Q.init(points);
    } else {
      Q.init(PLCpoints);
    }
  }

  m_keyEncode = Q.keyEncoding();
  QuantizedTessellation qmesh(points);
  m_QPLC.init(geometry, PLCpoints);
  qmesh.cullExternalPoints(m_QPLC);
  this->tessellateQuantized(qmesh);
  qmesh.clipTessellation(m_QPLC, m_serialTessellator);
  qmesh.filterToLocalGenerators();
  qmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  checkEncoding();
}

//------------------------------------------------------------------------------
// Partitioned unbounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
partitionAndTessellate(const std::vector<Point<Dimension, RealType>>& points,
                       const Partitioner<Dimension>& partitioner,
                       TessellationType& mesh) {
  POLY_ASSERT(mesh.empty());
  POLY_VERIFY2(partitioner.numPartitions() <= Communicator::getNRanks(),
               "Distributed partition count must not exceed the MPI rank count");

  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    Q.init(points);
  }

  m_keyEncode = Q.keyEncoding();
  const auto localPoints = partitioner.computeLocalPartition(points);
  QuantizedTessellation quantmesh(localPoints);
  this->tessellateQuantized(quantmesh);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  checkEncoding();
}

//------------------------------------------------------------------------------
// Partitioned PLC-bounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
partitionAndTessellate(const std::vector<Point<Dimension, RealType>>& points,
                       const std::vector<RealType>& PLCpoints,
                       const PLC<Dimension>& geometry,
                       const Partitioner<Dimension>& partitioner,
                       TessellationType& mesh) {
  m_clipping = true;
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(PLCpoints.size() % Dimension == 0);
  POLY_VERIFY2(partitioner.numPartitions() <= Communicator::getNRanks(),
               "Distributed partition count must not exceed the MPI rank count");

  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    if (PLCpoints.empty()) {
      Q.init(points);
    } else {
      Q.init(PLCpoints);
    }
  }

  m_keyEncode = Q.keyEncoding();
  m_QPLC.init(geometry, PLCpoints);
  const auto localPoints = partitioner.computeLocalPartition(points);
  QuantizedTessellation quantmesh(localPoints);
  quantmesh.cullExternalPoints(m_QPLC);
  this->tessellateQuantized(quantmesh);
  quantmesh.clipTessellation(m_QPLC, m_serialTessellator);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
  checkEncoding();
}

//------------------------------------------------------------------------------
// Quantized distributed tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
tessellateQuantizedImpl(QuantizedTessellation& qmesh) {
  int rank = Communicator::getRank();
  int nranks = Communicator::getNRanks();

  // Get the local visible generators and gather them on all processors
  auto visibleMesh = generateVisibleMesh(qmesh);
  std::set<int> neighborRanks;
  if (visibleMesh.points.size() > 1) {
    neighborRanks = visibleMesh.neighboringRanks();
    // If clipping, holes can modify which ranks are neighbors
    if (m_clipping) {
      visibleMesh.clipTessellation(m_QPLC, m_serialTessellator);
      std::set<int> clipped_neighbors = visibleMesh.neighboringRanks();
      neighborRanks.insert(clipped_neighbors.begin(), clipped_neighbors.end());
    }
  }

  // Check if any there is any overlap with other convex hulls
  if (!qmesh.convexHull.m_convex) {
    qmesh.makeConvexHull();
  }
  const auto localHull = qmesh.convexHull;
  const bool validLocalHull = localHull.isValid() &&
                              localHull.m_convex &&
                              localHull.facets.size() >= Dimension + 1;
  const auto allHulls = allGatherHulls(validLocalHull, localHull);
  POLY_ASSERT2(int(allHulls.size()) == nranks,
               "Incorrect number of hulls after communication");
  if (validLocalHull) {
    for (int source = 0; source < nranks; ++source) {
      if (source != rank &&
          allHulls[source].first &&
          QuantPLC<Dimension>::convexPLCIntersection(localHull, allHulls[source].second)) {
        neighborRanks.insert(source);
      }
    }
  }

  // Exchange complete neighbor generator sets using the selected representation.
  if (m_exchangePoints) {
    auto neighborGenerators =
      exchangeNeighborGenerators<Dimension, QuantizedPoint<Dimension>>(qmesh.points, neighborRanks);
    if (!qmesh.points.empty()) {
      qmesh.extend(neighborGenerators);
    }
  } else {
    auto neighborGenerators =
      exchangeNeighborGenerators<Dimension, QuantizedKey<Dimension>>(qmesh.hashes, neighborRanks);
    if (!qmesh.points.empty()) {
      qmesh.extend(neighborGenerators);
    }
  }
  // If we are clipping, we must include all visible generator points to avoid issues
  // See the DistributedHoleTests.cc for examples why this is needed
  if (m_clipping) {
    addVisibleMesh(visibleMesh, neighborRanks, qmesh);
  }
  m_serialTessellator.tessellateQuantized(qmesh);
}

//------------------------------------------------------------------------------
// Generate the visible mesh
//------------------------------------------------------------------------------
template<int Dimension>
QuantTessellation<Dimension>
DistributedTessellator<Dimension>::
generateVisibleMesh(QuantizedTessellation& qmesh) {
  if (!qmesh.points.empty()) {
    m_serialTessellator.tessellateQuantized(qmesh);
  }
  QuantizedTessellation visibleMesh;
  if (m_exchangePoints) {
    const auto visiblePoints = qmesh.visibleGeneratorPoints();
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
    const auto visibleHashes = qmesh.visibleGeneratorKeys();
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
               QuantizedTessellation& qmesh) {
  int rank = Communicator::getRank();
  int nranks = Communicator::getNRanks();
  if (m_exchangePoints) {
    std::vector<std::vector<QuantizedPoint<Dimension>>> visibleGenerators(nranks);
    const auto& vpoints = visibleMesh.points;
    POLY_ASSERT2(vpoints.size() == visibleMesh.cellRank.size(), "Incorrect cell ranks in visible mesh");
    for (auto i = 0u; i < vpoints.size(); ++i) {
      auto source = visibleMesh.cellRank[i];
      // Make sure the points are not part of the current or neighbor mesh to avoid repeats
      if (source != rank && neighborRanks.find(source) == neighborRanks.end()) {
        visibleGenerators[source].push_back(vpoints[i]);
      }
    }
    qmesh.extend(visibleGenerators);
  } else {
    std::vector<std::vector<QuantizedKey<Dimension>>> visibleGenerators(nranks);
    const auto& vhashes = visibleMesh.hashes;
    POLY_ASSERT2(vhashes.size() == visibleMesh.cellRank.size(), "Incorrect cell ranks in visible mesh");
    for (auto i = 0u; i < vhashes.size(); ++i) {
      auto source = visibleMesh.cellRank[i];
      // Make sure the points are not part of the current or neighbor mesh to avoid repeats
      if (source != rank && neighborRanks.find(source) == neighborRanks.end()) {
        visibleGenerators[source].push_back(vhashes[i]);
      }
    }
    qmesh.extend(visibleGenerators);
  }
}

template class DistributedTessellator<2>;
template class DistributedTessellator<3>;

} // namespace polytope
