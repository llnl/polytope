//----------------------------------------------------------------------------//
// ParallelUtils
//
// Provides utilities for doing tessellation in parallel
//----------------------------------------------------------------------------//
#ifndef __Polytope_ParallelUtils__
#define __Polytope_ParallelUtils__

#include <utility>

#include "Communicator.hh"
#include "Intersections.hh"
#include "QuantPLC.hh"
#include "findBoundaryElements.hh"
#include "Serializer.hh"

namespace polytope {

inline
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

// Gather a vector of vectors of the coordinate hashes.
// It is ordered result[rank][point index]
template<int Dimension, typename CoordHash>
std::vector<std::vector<CoordHash>>
allGatherGenerators(const std::vector<CoordHash>& localGenerators) {
  auto size = Communicator::getNProcs();
  std::vector<char> localBuffer;
  serialize(localGenerators, localBuffer);
  auto buffers = allGatherBuffers(localBuffer);
  std::vector<std::vector<CoordHash>> result(buffers.size());
  // Loop over source ranks so we can track which points came from where
  for (int source = 0; source < size; ++source) {
    auto itr = buffers[source].cbegin();
    std::vector<CoordHash> recvGens;
    deserialize<std::vector<CoordHash>>(recvGens, itr, buffers[source].end());
    result[source] = std::move(recvGens);
  }
  return result;
}

template<int Dimension, typename CoordHash>
std::vector<std::vector<CoordHash>>
exchangeNeighborGenerators(const std::vector<CoordHash>& localGenerators,
                           const std::set<int>& neighbors) {
  auto& comm = Communicator::communicator();
  auto rank = Communicator::getRank();
  auto size = Communicator::getNProcs();
  std::vector<char> localBuffer;
  serialize(localGenerators, localBuffer);

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

  std::vector<std::vector<CoordHash>> result(size);
  for (int source = 0; source < size; ++source) {
    auto itr = recvBuffers[source].cbegin();
    std::vector<CoordHash> received;
    if (recvBuffers[source].size() > 0) {
      deserialize<std::vector<CoordHash>>(received, itr, recvBuffers[source].end());
    }
    result[source] = std::move(received);
  }
  return result;
}

template<int Dimension>
std::vector<std::pair<bool, QuantPLC<Dimension>>>
allGatherHulls(const bool localValid,
               const QuantPLC<Dimension>& localHull) {
  auto size = Communicator::getNProcs();
  std::vector<char> localBuffer;
  serialize(localValid, localBuffer);
  if (localValid) {
    serialize(localHull, localBuffer);
  }
  auto buffers = allGatherBuffers(localBuffer);
  std::vector<std::pair<bool, QuantPLC<Dimension>>> result(size);
  for (int i = 0; i < size; ++i) {
    auto itr = buffers[i].cbegin();
    deserialize(result[i].first, itr, buffers[i].end());
    if (result[i].first) {
      deserialize<QuantPLC<Dimension>>(result[i].second, itr, buffers[i].end());
    }
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

}
#endif
