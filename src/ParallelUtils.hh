//----------------------------------------------------------------------------//
// ParallelUtils
//
// Provides utilities for doing tessellation in parallel
//----------------------------------------------------------------------------//
#ifndef __Polytope_ParallelUtils__
#define __Polytope_ParallelUtils__

#include <limits>
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
  auto size = Communicator::getNRanks();
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

// Exchange rank-indexed byte buffers.  A null send-buffer pointer denotes an
// empty message for that rank; this permits multiple destinations to share a
// single serialized buffer without copying it.
inline
std::vector<std::vector<char>>
exchangeBuffersByRank(const std::vector<const std::vector<char>*>& sendBuffers,
                      const int mpiTag) {
  auto& comm = Communicator::communicator();
  const auto nranks = Communicator::getNRanks();
  POLY_VERIFY2(sendBuffers.size() == static_cast<std::size_t>(nranks),
               "Expected one send buffer per MPI rank");

  std::vector<int> sendSizes(nranks, 0), recvSizes(nranks, 0);
  for (int destination = 0; destination < nranks; ++destination) {
    const auto* sendBuffer = sendBuffers[destination];
    if (sendBuffer != nullptr) {
      POLY_VERIFY2(sendBuffer->size() <=
                     static_cast<std::size_t>(std::numeric_limits<int>::max()),
                   "Generator exchange message exceeds MPI int count limit");
      sendSizes[destination] = static_cast<int>(sendBuffer->size());
    }
  }

  MPI_Alltoall(sendSizes.data(), 1, MPI_INT,
               recvSizes.data(), 1, MPI_INT,
               comm);

  std::vector<std::vector<char>> recvBuffers(nranks);
  std::vector<MPI_Request> requests;
  for (int source = 0; source < nranks; ++source) {
    if (recvSizes[source] > 0) {
      recvBuffers[source].resize(recvSizes[source]);
      requests.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(recvBuffers[source].data(), recvSizes[source], MPI_BYTE,
                source, mpiTag, comm, &requests.back());
    }
  }
  for (int destination = 0; destination < nranks; ++destination) {
    if (sendSizes[destination] > 0) {
      const auto& sendBuffer = *sendBuffers[destination];
      requests.push_back(MPI_REQUEST_NULL);
      MPI_Isend(sendBuffer.data(), sendSizes[destination], MPI_BYTE,
                destination, mpiTag, comm, &requests.back());
    }
  }
  if (!requests.empty()) {
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  }
  return recvBuffers;
}

// Gather serialized generator records, ordered result[rank][point index].
// Records may be Morton keys or quantized points.
template<int Dimension, typename Generator>
std::vector<std::vector<Generator>>
allGatherGenerators(const std::vector<Generator>& localGenerators) {
  auto size = Communicator::getNRanks();
  std::vector<char> localBuffer;
  serialize(localGenerators, localBuffer);
  auto buffers = allGatherBuffers(localBuffer);
  std::vector<std::vector<Generator>> result(buffers.size());
  // Loop over source ranks so we can track which points came from where
  for (int source = 0; source < size; ++source) {
    auto itr = buffers[source].cbegin();
    std::vector<Generator> recvGens;
    deserialize<std::vector<Generator>>(recvGens, itr, buffers[source].end());
    result[source] = std::move(recvGens);
  }
  return result;
}

// Exchange serialized generator records with the selected neighboring ranks.
// Records may be Morton keys or quantized points.
template<int Dimension, typename Generator>
std::vector<std::vector<Generator>>
exchangeNeighborGenerators(const std::vector<Generator>& localGenerators,
                           const std::set<int>& neighbors) {
  auto rank = Communicator::getRank();
  auto size = Communicator::getNRanks();
  std::vector<char> localBuffer;
  serialize(localGenerators, localBuffer);

  std::vector<const std::vector<char>*> sendBuffers(size, nullptr);
  for (const auto neighbor : neighbors) {
    if (neighbor != rank) {
      sendBuffers[neighbor] = &localBuffer;
    }
  }
  auto recvBuffers = exchangeBuffersByRank(sendBuffers, 9721);

  std::vector<std::vector<Generator>> result(size);
  for (int source = 0; source < size; ++source) {
    auto itr = recvBuffers[source].cbegin();
    std::vector<Generator> received;
    if (recvBuffers[source].size() > 0) {
      deserialize<std::vector<Generator>>(received, itr, recvBuffers[source].end());
    }
    result[source] = std::move(received);
  }
  return result;
}

// Redistribute generator records to their destination ranks.  Every rank must
// provide one record vector for every rank in the communicator; the returned
// vector contains all records sent to this rank, grouped by source rank.
// Unlike exchangeNeighborGenerators, each destination can receive a different
// subset of the local generators.
template<int Dimension, typename Generator>
std::vector<Generator>
redistributeGenerators(const std::vector<std::vector<Generator>>& generatorsByDestination) {
  const auto nranks = Communicator::getNRanks();
  POLY_VERIFY2(generatorsByDestination.size() == static_cast<std::size_t>(nranks),
               "Expected one generator list per MPI rank");

  std::vector<std::vector<char>> sendBuffers(nranks);
  for (int destination = 0; destination < nranks; ++destination) {
    const auto& generators = generatorsByDestination[destination];
    if (!generators.empty()) {
      serialize(generators, sendBuffers[destination]);
    }
  }

  std::vector<const std::vector<char>*> sendBufferPointers(nranks);
  for (int destination = 0; destination < nranks; ++destination) {
    sendBufferPointers[destination] = &sendBuffers[destination];
  }
  auto recvBuffers = exchangeBuffersByRank(sendBufferPointers, 9722);

  std::vector<Generator> result;
  for (int source = 0; source < nranks; ++source) {
    if (!recvBuffers[source].empty()) {
      auto itr = recvBuffers[source].cbegin();
      deserialize(result, itr, recvBuffers[source].cend());
    }
  }
  return result;
}

template<int Dimension>
std::vector<std::pair<bool, QuantPLC<Dimension>>>
allGatherHulls(const bool localValid,
               const QuantPLC<Dimension>& localHull) {
  auto size = Communicator::getNRanks();
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
