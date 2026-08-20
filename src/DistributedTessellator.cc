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

  auto neighborGenerators =
    exchangeNeighborGenerators<Dimension, MortonKey<Dimension>>(qmesh.hashes, neighborRanks);

  // Extend the QuantTessellation by it's neighbor generators and retessellate
  if (!qmesh.points.empty()) {
    qmesh.init(neighborGenerators);
    m_serialTessellator.tessellateQuantized(qmesh);
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
  auto visibleHashes = qmesh.visibleGenerators();
  auto allVisibleRecords = allGatherGenerators<Dimension, MortonKey<Dimension>>(visibleHashes);
  size_t nvisible = 0;
  for (const auto& rankHashes : allVisibleRecords) {
    nvisible += rankHashes.size();
  }
  QuantizedTessellation visibleMesh;
  if (nvisible > 0) {
    visibleMesh.init(allVisibleRecords);
    m_serialTessellator.tessellateQuantized(visibleMesh);
  }
  return visibleMesh;
}

template class DistributedTessellator<2>;
template class DistributedTessellator<3>;

} // namespace polytope
