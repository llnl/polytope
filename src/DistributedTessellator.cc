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
  m_serialTessellator(serialTessellator) {
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
tessellate(const std::vector<RealType>& points,
           TessellationType& mesh) {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() % Dimension == 0);

  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    typename Quantizer<Dimension>::RealPoint globalMin, globalMax;
    findGlobalBounds<Dimension>(points, globalMin, globalMax);
    Q.init(globalMin, globalMax);
  }

  QuantizedTessellation quantmesh(points);
  this->tessellateQuantized(quantmesh);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
}

//------------------------------------------------------------------------------
// PLC-bounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
tessellate(const std::vector<RealType>& points,
           const std::vector<RealType>& PLCpoints,
           const PLC<Dimension>& geometry,
           TessellationType& mesh) {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() % Dimension == 0);
  POLY_ASSERT(PLCpoints.size() % Dimension == 0);

  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    const auto& boundsPoints = PLCpoints.empty() ? points : PLCpoints;
    typename Quantizer<Dimension>::RealPoint globalMin, globalMax;
    findGlobalBounds<Dimension>(boundsPoints, globalMin, globalMax);
    Q.init(globalMin, globalMax);
  }

  QuantizedTessellation qmesh(points);
  QuantPLC<Dimension> qplc(geometry, PLCpoints);
  qmesh.cullExternalPoints(qplc);
  this->tessellateQuantized(qmesh);
  qmesh.clipTessellation(qplc, m_serialTessellator);
  qmesh.filterToLocalGenerators();
  qmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
}

//------------------------------------------------------------------------------
// Partitioned unbounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
partitionAndTessellate(const std::vector<RealType>& points,
                       const Partitioner<Dimension>& partitioner,
                       TessellationType& mesh) {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() % Dimension == 0);

  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    typename Quantizer<Dimension>::RealPoint globalMin, globalMax;
    findGlobalBounds<Dimension>(points, globalMin, globalMax);
    Q.init(globalMin, globalMax);
  }

  QuantizedTessellation quantmesh;
  {
    QuantizedTessellation replicatedMesh(points);
    quantmesh.init(partitioner.computePartition(replicatedMesh.getQuantizedPoints()));
  }
  this->tessellateQuantized(quantmesh);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
}

//------------------------------------------------------------------------------
// Partitioned PLC-bounded tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
partitionAndTessellate(const std::vector<RealType>& points,
                       const std::vector<RealType>& PLCpoints,
                       const PLC<Dimension>& geometry,
                       const Partitioner<Dimension>& partitioner,
                       TessellationType& mesh) {
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() % Dimension == 0);
  POLY_ASSERT(PLCpoints.size() % Dimension == 0);

  auto& Q = Quantizer<Dimension>::instance();
  if (!Q.m_init) {
    const auto& boundsPoints = PLCpoints.empty() ? points : PLCpoints;
    typename Quantizer<Dimension>::RealPoint globalMin, globalMax;
    findGlobalBounds<Dimension>(boundsPoints, globalMin, globalMax);
    Q.init(globalMin, globalMax);
  }

  QuantizedTessellation quantmesh;
  QuantPLC<Dimension> qplc(geometry, PLCpoints);
  {
    QuantizedTessellation replicatedMesh(points);
    replicatedMesh.cullExternalPoints(qplc);
    quantmesh.init(partitioner.computePartition(replicatedMesh.getQuantizedPoints()));
  }
  this->tessellateQuantized(quantmesh);
  quantmesh.clipTessellation(qplc, m_serialTessellator);
  quantmesh.filterToLocalGenerators();
  quantmesh.fillTessellation(mesh);
  findBoundaryElements(mesh, mesh.boundaryFaces, mesh.boundaryNodes);
}

//------------------------------------------------------------------------------
// Quantized distributed tessellation
//------------------------------------------------------------------------------
template<int Dimension>
void
DistributedTessellator<Dimension>::
tessellateQuantizedImpl(QuantizedTessellation& qmesh) {
  int rank = Communicator::getRank();
  int nranks = Communicator::getNProcs();

  // Get the local visible generators and gather them on all processors
  auto visibleMesh = generateVisibleMesh(qmesh);
  std::set<int> neighborRanks;
  if (visibleMesh.points.size() > 1) {
    neighborRanks = visibleMesh.neighboringRanks();
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
      qmesh.init(neighborGenerators);
      m_serialTessellator.tessellateQuantized(qmesh);
    }
  } else {
    auto neighborGenerators =
      exchangeNeighborGenerators<Dimension, MortonKey<Dimension>>(qmesh.hashes, neighborRanks);
    if (!qmesh.points.empty()) {
      qmesh.init(neighborGenerators);
      m_serialTessellator.tessellateQuantized(qmesh);
    }
  }
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
      visibleMesh.init(allVisibleRecords);
      m_serialTessellator.tessellateQuantized(visibleMesh);
    }
  } else {
    const auto visibleHashes = qmesh.visibleGeneratorKeys();
    auto allVisibleRecords =
      allGatherGenerators<Dimension, MortonKey<Dimension>>(visibleHashes);
    size_t nvisible = 0;
    for (const auto& rankHashes : allVisibleRecords) {
      nvisible += rankHashes.size();
    }
    if (nvisible > 0) {
      visibleMesh.init(allVisibleRecords);
      m_serialTessellator.tessellateQuantized(visibleMesh);
    }
  }
  return visibleMesh;
}

template class DistributedTessellator<2>;
template class DistributedTessellator<3>;

} // namespace polytope
