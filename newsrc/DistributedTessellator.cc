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

#include "Intersections.hh"
#include "QuantPLC.hh"
#include "findBoundaryElements.hh"
#include "Serializer.hh"

namespace polytope {

namespace {

template<int Dimension>
void
packGenerators(const std::vector<GeneratorRecord<Dimension>>& generators,
               std::vector<char>& buffer) {
  const auto n = static_cast<unsigned>(generators.size());
  serialize(n, buffer);
  for (const auto& g : generators) {
    serialize(g.rank, buffer);
    serialize(g.ordinal, buffer);
    serialize(g.point, buffer);
    serialize(g.hash, buffer);
  }
}

template<int Dimension>
std::vector<GeneratorRecord<Dimension>>
unpackGenerators(const std::vector<char>& buffer) {
  std::vector<GeneratorRecord<Dimension>> result;
  auto itr = buffer.begin();
  unsigned n = 0;
  deserialize(n, itr, buffer.end());
  result.resize(n);
  for (auto& g : result) {
    deserialize(g.rank, itr, buffer.end());
    deserialize(g.ordinal, itr, buffer.end());
    deserialize(g.point, itr, buffer.end());
    deserialize(g.hash, itr, buffer.end());
  }
  return result;
}

template<int Dimension>
std::vector<std::vector<char>>
allGatherBuffers(const std::vector<char>& localBuffer) {
  auto& comm = Communicator::communicator();
  auto size = Communicator::getNProcs();
  const auto localSize = static_cast<int>(localBuffer.size());
  std::vector<int> recvSizes(size, 0);
  MPI_Allgather(&localSize, 1, MPI_INT,
                recvSizes.data(), 1, MPI_INT,
                comm);

  std::vector<int> displs(size, 0);
  int totalSize = 0;
  for (int i = 0; i < size; ++i) {
    displs[i] = totalSize;
    totalSize += recvSizes[i];
  }

  std::vector<char> recvBuffer(totalSize);
  MPI_Allgatherv(localBuffer.data(), localSize, MPI_BYTE,
                 recvBuffer.data(), recvSizes.data(), displs.data(), MPI_BYTE,
                 comm);

  std::vector<std::vector<char>> result(size);
  for (int i = 0; i < size; ++i) {
    result[i].assign(recvBuffer.begin() + displs[i],
                     recvBuffer.begin() + displs[i] + recvSizes[i]);
  }
  return result;
}

template<int Dimension>
std::vector<GeneratorRecord<Dimension>>
allGatherGenerators(const std::vector<GeneratorRecord<Dimension>>& localGenerators) {
  std::vector<char> localBuffer;
  packGenerators(localGenerators, localBuffer);

  std::vector<GeneratorRecord<Dimension>> result;
  auto buffers = allGatherBuffers<Dimension>(localBuffer);
  for (const auto& buffer : buffers) {
    auto generators = unpackGenerators<Dimension>(buffer);
    result.insert(result.end(), generators.begin(), generators.end());
  }
  return result;
}

template<int Dimension>
std::vector<GeneratorRecord<Dimension>>
exchangeNeighborGenerators(const std::vector<GeneratorRecord<Dimension>>& localGenerators,
                           const std::set<int>& neighbors) {
  auto& comm = Communicator::communicator();
  auto rank = Communicator::getRank();
  auto size = Communicator::getNProcs();
  std::vector<char> localBuffer;
  packGenerators(localGenerators, localBuffer);

  std::vector<int> sendSizes(size, 0), recvSizes(size, 0);
  for (const auto neighbor : neighbors) {
    if (neighbor != rank) {
      sendSizes[neighbor] = static_cast<int>(localBuffer.size());
    }
  }
  MPI_Alltoall(sendSizes.data(), 1, MPI_INT,
               recvSizes.data(), 1, MPI_INT,
               comm);

  std::vector<std::vector<char>> recvBuffers(size);
  std::vector<MPI_Request> requests;
  for (int r = 0; r < size; ++r) {
    if (recvSizes[r] > 0) {
      recvBuffers[r].resize(recvSizes[r]);
      requests.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(recvBuffers[r].data(), recvSizes[r], MPI_BYTE,
                r, 9721, comm, &requests.back());
    }
  }
  for (int r = 0; r < size; ++r) {
    if (sendSizes[r] > 0) {
      requests.push_back(MPI_REQUEST_NULL);
      MPI_Isend(localBuffer.data(), sendSizes[r], MPI_BYTE,
                r, 9721, comm, &requests.back());
    }
  }
  if (!requests.empty()) {
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  }

  std::vector<GeneratorRecord<Dimension>> result;
  for (const auto& buffer : recvBuffers) {
    if (!buffer.empty()) {
      auto generators = unpackGenerators<Dimension>(buffer);
      result.insert(result.end(), generators.begin(), generators.end());
    }
  }
  return result;
}

template<int Dimension>
std::vector<GeneratorRecord<Dimension>>
recordsFromQuantTessellation(const QuantTessellation<Dimension>& qmesh) {
  auto rank = Communicator::getRank();
  std::vector<GeneratorRecord<Dimension>> result;
  result.reserve(qmesh.points.size());
  for (unsigned i = 0; i < qmesh.points.size(); ++i) {
    GeneratorRecord<Dimension> record;
    record.rank = rank;
    record.ordinal = i;
    record.point = qmesh.points[i];
    record.hash = qmesh.hashes[i];
    result.push_back(record);
  }
  return result;
}

template<int Dimension>
std::vector<typename Quantizer<Dimension>::IntPoint>
pointsFromRecords(const std::vector<GeneratorRecord<Dimension>>& records) {
  std::vector<typename Quantizer<Dimension>::IntPoint> result;
  result.reserve(records.size());
  for (auto record : records) {
    record.point.index = result.size();
    result.push_back(record.point);
  }
  return result;
}

template<int Dimension>
void
checkUniqueGeneratorHashes(const std::vector<GeneratorRecord<Dimension>>& records,
                           const std::string& context) {
  using CoordHash = typename Quantizer<Dimension>::CoordHash;
  using HashType = typename HashKey<Dimension>::HashType;

  std::unordered_map<CoordHash, std::pair<int, unsigned>, HashType> owners;
  for (const auto& record : records) {
    const auto [itr, inserted] =
      owners.emplace(record.hash, std::make_pair(record.rank, record.ordinal));
    if (!inserted) {
      const auto& owner = itr->second;
      POLY_CHECK2(owner.first == record.rank && owner.second == record.ordinal,
                  "Duplicate quantized generator hash in " << context
                  << " between rank " << owner.first << " ordinal " << owner.second
                  << " and rank " << record.rank << " ordinal " << record.ordinal);
    }
  }
}

template<int Dimension>
bool
makeConvexHull(const QuantTessellation<Dimension>& qmesh,
               QuantPLC<Dimension>& hull) {
  if (qmesh.points.size() < Dimension + 1) return false;

  PLC<Dimension> emptyPLC;
  hull.init(emptyPLC, qmesh.points);
  hull.makeConvex();
  return hull.m_convex && hull.facets.size() >= Dimension + 1;
}

template<int Dimension>
void
packHull(const bool valid,
         const QuantPLC<Dimension>& hull,
         std::vector<char>& buffer) {
  const int validInt = valid ? 1 : 0;
  serialize(validInt, buffer);
  if (!valid) return;

  serialize(hull.facets, buffer);
  serialize(hull.holes, buffer);
  serialize(hull.points, buffer);
  serialize(hull.hashes, buffer);
  serialize(hull.m_normals, buffer);
}

template<int Dimension>
bool
unpackHull(const std::vector<char>& buffer,
           QuantPLC<Dimension>& hull) {
  auto itr = buffer.begin();
  int validInt = 0;
  deserialize(validInt, itr, buffer.end());
  if (validInt == 0) return false;

  deserialize(hull.facets, itr, buffer.end());
  deserialize(hull.holes, itr, buffer.end());
  deserialize(hull.points, itr, buffer.end());
  deserialize(hull.hashes, itr, buffer.end());
  deserialize(hull.m_normals, itr, buffer.end());
  hull.m_convex = true;
  hull.m_reduced = true;
  return true;
}

template<int Dimension>
std::vector<std::pair<bool, QuantPLC<Dimension>>>
allGatherHulls(const bool localValid,
               const QuantPLC<Dimension>& localHull) {
  auto size = Communicator::getNProcs();
  std::vector<char> localBuffer;
  packHull(localValid, localHull, localBuffer);
  auto buffers = allGatherBuffers<Dimension>(localBuffer);

  std::vector<std::pair<bool, QuantPLC<Dimension>>> result(size);
  for (int i = 0; i < size; ++i) {
    result[i].first = unpackHull(buffers[i], result[i].second);
  }
  return result;
}

template<int Dimension>
void
findGlobalBounds(const std::vector<double>& coords,
                 typename Quantizer<Dimension>::RealPoint& globalMin,
                 typename Quantizer<Dimension>::RealPoint& globalMax) {
  auto& comm = Communicator::communicator();
  using RealPoint = typename Quantizer<Dimension>::RealPoint;
  using RealType = typename Quantizer<Dimension>::RealType;

  std::array<RealType, Dimension> localMin, localMax, reducedMin, reducedMax;
  localMin.fill(std::numeric_limits<RealType>::max());
  localMax.fill(-std::numeric_limits<RealType>::max());

  const auto points = extractCoords<Dimension, RealType>(coords);
  for (const auto& p : points) {
    for (int d = 0; d < Dimension; ++d) {
      localMin[d] = std::min(localMin[d], p[d]);
      localMax[d] = std::max(localMax[d], p[d]);
    }
  }

  MPI_Allreduce(localMin.data(), reducedMin.data(), Dimension, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(localMax.data(), reducedMax.data(), Dimension, MPI_DOUBLE, MPI_MAX, comm);

  globalMin = RealPoint();
  globalMax = RealPoint();
  for (int d = 0; d < Dimension; ++d) {
    globalMin[d] = reducedMin[d];
    globalMax[d] = reducedMax[d];
  }
}

template<int Dimension>
std::vector<std::vector<unsigned>>
computeFaceCells(const QuantTessellation<Dimension>& qmesh) {
  std::vector<std::vector<unsigned>> faceCells(qmesh.faces.size());
  for (unsigned cellID = 0; cellID < qmesh.cells.size(); ++cellID) {
    for (auto faceID : qmesh.cells[cellID]) {
      const auto absFaceID = faceID < 0 ? ~faceID : faceID;
      POLY_ASSERT(absFaceID < faceCells.size());
      faceCells[absFaceID].push_back(cellID);
    }
  }
  return faceCells;
}

template<int Dimension>
std::set<int>
neighborRanksFromVisibleVoronoi(const QuantTessellation<Dimension>& visibleMesh,
                                const std::vector<GeneratorRecord<Dimension>>& visibleRecords) {
  auto rank = Communicator::getRank();
  std::set<int> neighbors;
  const auto faceCells = computeFaceCells(visibleMesh);
  for (const auto& cells : faceCells) {
    if (cells.size() < 2) continue;
    for (unsigned i = 0; i < cells.size(); ++i) {
      const auto ranki = visibleRecords[cells[i]].rank;
      for (unsigned j = i + 1; j < cells.size(); ++j) {
        const auto rankj = visibleRecords[cells[j]].rank;
        if (ranki == rank && rankj != rank) neighbors.insert(rankj);
        if (rankj == rank && ranki != rank) neighbors.insert(ranki);
      }
    }
  }
  return neighbors;
}

template<int Dimension>
void
filterToLocalGenerators(QuantTessellation<Dimension>& qmesh,
                        const std::vector<GeneratorRecord<Dimension>>& records) {
  auto rank = Communicator::getRank();
  POLY_ASSERT(records.size() == qmesh.points.size());
  POLY_ASSERT(qmesh.cells.size() == qmesh.points.size());

  std::vector<typename QuantTessellation<Dimension>::IntPoint> newPoints;
  std::vector<typename QuantTessellation<Dimension>::CoordHash> newHashes;
  std::vector<std::vector<int>> newCells;
  newPoints.reserve(qmesh.points.size());
  newHashes.reserve(qmesh.hashes.size());
  newCells.reserve(qmesh.cells.size());

  for (unsigned i = 0; i < records.size(); ++i) {
    if (records[i].rank == rank) {
      auto point = qmesh.points[i];
      point.index = newPoints.size();
      newPoints.push_back(point);
      newHashes.push_back(qmesh.hashes[i]);
      newCells.push_back(qmesh.cells[i]);
    }
  }

  qmesh.points = std::move(newPoints);
  qmesh.hashes = std::move(newHashes);
  qmesh.cells = std::move(newCells);
}

template<int Dimension>
std::vector<GeneratorRecord<Dimension>>
visibleLocalGenerators(const QuantTessellation<Dimension>& localMesh,
                       const std::vector<GeneratorRecord<Dimension>>& localRecords,
                       const bool validHull,
                       const QuantPLC<Dimension>& localHull) {
  if constexpr (Dimension == 2) {
    std::vector<GeneratorRecord<Dimension>> result;
    result.reserve(localRecords.size());
    for (unsigned i = 0; i < localRecords.size(); ++i) {
      if (!validHull || localMesh.cellIntersectsHull(localHull, i)) {
        result.push_back(localRecords[i]);
      }
    }
    return result;
  } else {
    POLY_CHECK2(false, "DistributedTessellator visibility is implemented for 2D only");
    return std::vector<GeneratorRecord<Dimension>>();
  }
}

} // anonymous namespace

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

