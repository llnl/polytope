#include "DistributedTessellator.hh"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <fstream>

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
  std::ostringstream os;
  os << "DistributedTessellator_"
     << m_serialTessellator.name();
  return os.str();
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

  typename Quantizer<Dimension>::RealPoint globalMin, globalMax;
  findGlobalBounds<Dimension>(points, globalMin, globalMax);
  Quantizer<Dimension>::instance().init(globalMin, globalMax);

  QuantizedTessellation quantmesh(points);
  this->tessellateQuantized(quantmesh);
  filterToLocalGenerators(quantmesh, m_localRecords);
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

  const auto& boundsPoints = PLCpoints.empty() ? points : PLCpoints;
  typename Quantizer<Dimension>::RealPoint globalMin, globalMax;
  findGlobalBounds<Dimension>(boundsPoints, globalMin, globalMax);
  Quantizer<Dimension>::instance().init(globalMin, globalMax);

  QuantizedTessellation quantmesh(points);
  QuantPLC<Dimension> qplc(geometry, PLCpoints);
  quantmesh.cullExternalPoints(qplc);

  this->tessellateQuantized(quantmesh);
  quantmesh.clipTessellation(qplc, m_serialTessellator);
  filterToLocalGenerators(quantmesh, m_localRecords);
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
  auto rank = Communicator::getRank();
  auto size = Communicator::getNProcs();
  auto localRecords = recordsFromQuantTessellation(qmesh);
  checkUniqueGeneratorHashes(localRecords, "local generator set");

  if (!localRecords.empty()) {
    m_serialTessellator.tessellateQuantized(qmesh);
  }

  QuantPLC<Dimension> localHull;
  const bool validLocalHull = makeConvexHull(qmesh, localHull);

  auto visibleRecords =
    visibleLocalGenerators(qmesh, localRecords, validLocalHull, localHull);
  auto allVisibleRecords =
    allGatherGenerators(visibleRecords);
  checkUniqueGeneratorHashes(allVisibleRecords, "visible generator set");

  std::set<int> neighbors;
  const auto allHulls = allGatherHulls(validLocalHull, localHull);
  if (validLocalHull) {
    for (int r = 0; r < size; ++r) {
      if (r != rank && allHulls[r].first &&
          QuantPLC<Dimension>::convexPLCIntersection(localHull, allHulls[r].second)) {
        neighbors.insert(r);
      }
    }
  }

  if (allVisibleRecords.size() > 1) {
    auto visiblePoints = pointsFromRecords(allVisibleRecords);
    QuantizedTessellation visibleMesh(visiblePoints);
    m_serialTessellator.tessellateQuantized(visibleMesh);
    auto visibleNeighbors =
      neighborRanksFromVisibleVoronoi(visibleMesh, allVisibleRecords);
    neighbors.insert(visibleNeighbors.begin(), visibleNeighbors.end());
  }

  auto neighborRecords =
    exchangeNeighborGenerators(localRecords, neighbors);

  std::vector<GeneratorRecord<Dimension>> finalRecords(localRecords);
  finalRecords.insert(finalRecords.end(), neighborRecords.begin(), neighborRecords.end());
  checkUniqueGeneratorHashes(finalRecords, "final local plus neighbor generator set");

  if (finalRecords.empty()) {
    return;
  }

  if (finalRecords.size() == localRecords.size()) {
    return;
  }

  auto finalPoints = pointsFromRecords(finalRecords);
  QuantizedTessellation finalMesh(finalPoints);
  m_serialTessellator.tessellateQuantized(finalMesh);
  // Retain the records to filter to only local generators after clipping
  m_localRecords = std::move(finalRecords);
  qmesh = std::move(finalMesh);
}

template class DistributedTessellator<2>;
template class DistributedTessellator<3>;

} // namespace polytope

