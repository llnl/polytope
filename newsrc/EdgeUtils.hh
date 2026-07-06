#ifndef POLYTOPE_EDGEUTILS_HH
#define POLYTOPE_EDGEUTILS_HH

// Set of logic for handling indices to edge vertices
// Useful for face merging

#include "Point.hh"
#include "HashKey.hh"

namespace polytope {
namespace edge {
//------------------------------------------------------------------------------
// Utilities for edges specifically
//------------------------------------------------------------------------------
using Edge = std::pair<int, int>;
//------------------------------------------------------------------------------
// Custom hashing function for pairs
//------------------------------------------------------------------------------
struct EdgeHash {
  std::size_t operator()(const Edge& p) const {
    auto x = p.first;
    auto y = p.second;
    if (x > y) {
      x = p.second;
      y = p.first;
    }
    auto h1 = std::hash<int>{}(x);
    auto h2 = std::hash<int>{}(y);
    return h1 ^ (h2 << 1);
  }
};

using EdgeToFaceMap = std::unordered_map<Edge, int, EdgeHash>;

//------------------------------------------------------------------------------
// Make ordered pair for edges
//------------------------------------------------------------------------------
inline Edge orderEdge(const int v0, const int v1) {
  return v0 < v1 ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
}

inline Edge orderEdge(const Edge edge) {
  return edge.first < edge.second ? edge : std::make_pair(edge.second, edge.first);
}

//------------------------------------------------------------------------------
// Order a loop of edges to form a connected chain
// Ensures edges[i][1] connects to edges[i+1][0]
//------------------------------------------------------------------------------
inline void orderEdgeLoop(std::vector<edge::Edge>& edges,
                          std::map<int, int>& oldToNew) {
  if (edges.empty()) return;

  std::vector<edge::Edge> ordered;
  ordered.reserve(edges.size());
  oldToNew.clear();

  // Build map: start vertex -> list of edge indices starting at that vertex
  std::map<int, std::vector<int>> startMap;
  for (size_t i = 0; i < edges.size(); ++i) {
    startMap[edges[i].first].push_back(i);
  }

  // Follow the chain starting from first edge
  std::set<int> used;
  int current = 0;

  while (used.size() < edges.size()) {
    // Add current edge to ordered list
    oldToNew[current] = ordered.size();
    ordered.push_back(edges[current]);
    used.insert(current);

    // Find next edge: one that starts where this one ends
    int nextVertex = edges[current].second;
    bool foundNext = false;

    if (startMap.count(nextVertex)) {
      for (int candidate : startMap[nextVertex]) {
        if (!used.count(candidate)) {
          current = candidate;
          foundNext = true;
          break;
        }
      }
    }

    // If chain is broken but we haven't used all edges, find an unused edge to continue
    if (!foundNext && used.size() < edges.size()) {
      for (size_t i = 0; i < edges.size(); ++i) {
        if (!used.count(i)) {
          current = i;
          foundNext = true;
          break;
        }
      }
    }

    // If we still can't find an edge, we're done
    if (!foundNext) break;
  }

  edges = ordered;
}

inline void orderEdgeLoop(std::vector<std::vector<int>>& edges,
                          std::map<int, int>& oldToNew) {
  if (edges.empty()) return;

  std::vector<std::vector<int>> ordered;
  ordered.reserve(edges.size());
  oldToNew.clear();
  

  // Build map: start vertex -> edge index
  std::map<int, int> startMap;
  for (size_t i = 0; i < edges.size(); ++i) {
    startMap[edges[i][0]] = i;
  }

  // Follow the chain starting from first edge
  std::set<int> used;
  int current = 0;
  while (used.size() < edges.size()) {
    oldToNew[current] = ordered.size();
    ordered.push_back(edges[current]);
    used.insert(current);

    int nextVertex = edges[current][1];
    if (startMap.count(nextVertex) && !used.count(startMap[nextVertex])) {
      current = startMap[nextVertex];
    } else {
      break;  // Chain broken
    }
  }

  edges = ordered;
}

inline void orderEdgeLoop(std::vector<std::vector<int>>& edges) {
  std::map<int, int> dummy;
  orderEdgeLoop(edges, dummy);
}

//------------------------------------------------------------------------------
// Convert an EdgeToFaceMap into an ordered facet.
// Uses deterministic starting point and sorted neighbors for consistent ordering
//------------------------------------------------------------------------------
inline std::vector<int> traceBoundary(const EdgeToFaceMap& uniqueEdges) {
  // Determine boundary faces since they only appear once
  std::unordered_map<int, std::vector<int>> adjacency;
  for (const auto& [edge, count] : uniqueEdges) {
    if (count == 1) {
      adjacency[edge.first].push_back(edge.second);
      adjacency[edge.second].push_back(edge.first);
    }
  }
  std::vector<int> boundary;
  if (adjacency.empty()) return boundary;

  // Sort neighbors for deterministic traversal
  for (auto& [vertex, neighbors] : adjacency) {
    std::sort(neighbors.begin(), neighbors.end());
  }
  // Pick smallest vertex as start for deterministic ordering
  int start = std::min_element(adjacency.begin(), adjacency.end(),
                                [](const auto& a, const auto& b) {
                                  return a.first < b.first;
                                })->first;
  int current = start;
  int prev = -1;

  do {
    boundary.push_back(current);
    const auto& neighbors = adjacency.at(current);
    auto it = std::find_if(neighbors.begin(), neighbors.end(),
                           [prev](int n) { return n != prev; });
    // Check if we found a valid neighbor
    if (it == neighbors.end()) {
      // Boundary is not a simple closed loop - return what we have
      break;
    }
    prev = current;
    current = *it;
  } while (current != start && boundary.size() < adjacency.size() + 1);
  return boundary;
}

//------------------------------------------------------------------------------
// Adds any unique edges to a set.
//------------------------------------------------------------------------------
inline void addUniqueEdges(const std::vector<int>& facet,
                           EdgeToFaceMap& uniqueEdges) {
  const auto N = facet.size();
  for (size_t i = 0; i < N; ++i) {
    int v0 = facet[i];
    int v1 = facet[(i + 1) % N];
    auto edge = orderEdge(v0, v1);
    uniqueEdges[edge]++;
  }
}

//------------------------------------------------------------------------------
// Check if facet shares edges with a unique set.
//------------------------------------------------------------------------------
inline bool sharedEdges(const std::vector<int>& facet,
                        EdgeToFaceMap& uniqueEdges) {
  const auto N = facet.size();
  bool shared = false;
  for (size_t i = 0; i < N; ++i) {
    int v0 = facet[i];
    int v1 = facet[(i + 1) % N];
    auto edge = orderEdge(v0, v1);
    if (uniqueEdges.count(edge)) {
      uniqueEdges.erase(edge);
      shared = true;
    } else {
      uniqueEdges[edge]++;
    }
  }
  return shared;
}

//------------------------------------------------------------------------------
// Edge storage and orientation tracking
//------------------------------------------------------------------------------

// Map from canonical edge to its face/edge index
using EdgeToFaceMap = std::unordered_map<Edge, int, EdgeHash>;

//------------------------------------------------------------------------------
// Add an oriented edge to the edge map
// Returns the signed face index:
//   - Positive if edge orientation matches canonical form
//   - Negative (bitwise NOT) if edge orientation is reversed
//------------------------------------------------------------------------------
inline int addOrientedEdge(int n0, int n1,
                           std::vector<std::vector<int>>& faces,
                           EdgeToFaceMap& edgeToFace) {
  Edge canonical = orderEdge(n0, n1);

  auto it = edgeToFace.find(canonical);
  int faceIndex;

  if (it == edgeToFace.end()) {
    // New edge - add to faces in canonical form
    faceIndex = faces.size();
    edgeToFace[canonical] = faceIndex;
    faces.push_back({canonical.first, canonical.second});
  } else {
    faceIndex = it->second;
  }

  // Return signed index based on whether orientation matches canonical
  return (canonical.first == n0) ? faceIndex : ~faceIndex;
}

//------------------------------------------------------------------------------
// Get the absolute (unsigned) face index from a signed index
//------------------------------------------------------------------------------
inline int unsignedIndex(int signedIndex) {
  return (signedIndex < 0) ? ~signedIndex : signedIndex;
}

//------------------------------------------------------------------------------
// Check if a signed index indicates reversed orientation
//------------------------------------------------------------------------------
inline bool isReversed(int signedIndex) {
  return signedIndex < 0;
}

//------------------------------------------------------------------------------
// Get the node indices for an edge, respecting signed orientation
// If signedIndex < 0, returns nodes in reverse order
//------------------------------------------------------------------------------
inline std::pair<int, int> getOrientedNodes(int signedIndex,
                                            const std::vector<std::vector<int>>& faces) {
  int faceIndex = unsignedIndex(signedIndex);
  const auto& face = faces[faceIndex];

  if (isReversed(signedIndex)) {
    return {face[1], face[0]};
  } else {
    return {face[0], face[1]};
  }
}

//------------------------------------------------------------------------------
// Reverse the orientation of a signed edge index
//------------------------------------------------------------------------------
inline int reverseOrientation(int signedIndex) {
  return ~signedIndex;
}

//------------------------------------------------------------------------------
// Check if there are nearly duplicate nodes
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
inline bool hasNearDuplicates(const Point<Dimension, CoordType>& point,
                              std::map<Point<Dimension, CoordType>, int>& node2id) {
  for (int offset : {-1, 1}) {
    for (int dim = 0; dim < Dimension; ++dim) {
      Point<Dimension, CoordType> pp(point);
      pp[dim] += offset;
      auto it = node2id.find(pp);
      if (it != node2id.end()) {
        return true;
      }
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Modify the nodes list if points do not exist in a given node id map
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
inline edge::Edge updateNodeMap(const Point<Dimension, CoordType>& p0,
                                const Point<Dimension, CoordType>& p1,
                                std::map<Point<Dimension, CoordType>, int>& node2id,
                                std::vector<Point<Dimension, CoordType>>& nodes) {
  auto it0 = node2id.find(p0);
  int n0;
  if (it0 == node2id.end()) {
    n0 = nodes.size();
    node2id[p0] = n0;
    nodes.push_back(p0);
  } else {
    n0 = it0->second;
  }
  auto it1 = node2id.find(p1);
  int n1;
  if (it1 == node2id.end()) {
    n1 = nodes.size();
    node2id[p1] = n1;
    nodes.push_back(p1);
  } else {
    n1 = it1->second;
  }
  return edge::Edge(std::make_pair(n0, n1));
}

//------------------------------------------------------------------------------
// Utilities for edge data, meaning edges paired with generator points
// This allows us to keep track of which edges belong to which generators
//------------------------------------------------------------------------------
using GenPair = std::pair<int, int>;
inline GenPair orderPair(const int a, const int b) {
  return orderEdge(a, b);
}

struct EdgeData {
  edge::Edge curEdge;
  int startSide, endSide;
};

using GenPairToEdgeDataMap = std::map<GenPair, EdgeData>;
}
}
#endif
