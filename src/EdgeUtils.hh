#ifndef __Polytope_EdgeUtils__
#define __Polytope_EdgeUtils__

// Set of logic for handling indices to edge vertices
// Useful for face merging

#include "Point.hh"
#include "QuantizedKeyTraits.hh"
#include "GeomUtils.hh"
#include "Shapes.hh"

namespace polytope {
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
// A directed edge and the clipping-box sides associated with its start and end
// nodes.  Keep these together whenever the edge loop is reordered.
//------------------------------------------------------------------------------
struct ClippedEdge {
  Edge curEdge;
  std::pair<int, int> clippedSides = std::make_pair(-1, -1);
};

inline void walkBoxEdges(const BoxSide& startSide,
                         const BoxSide& endSide,
                         const unsigned& startPoint,
                         const unsigned& endPoint,
                         const std::map<BoxSide, unsigned>& cornerIndices,
                         std::vector<Edge>& edges) {
  BoxSides sides;
  BoxSide thisSide = startSide;
  unsigned curPoint = startPoint;
  POLY_ASSERT(static_cast<int>(thisSide) >= 0);
  POLY_ASSERT(static_cast<int>(endSide) >= 0);
  while (thisSide != endSide) {
    if (isCorner(thisSide)) {
      unsigned nextPoint = cornerIndices.at(thisSide);
      if (curPoint != nextPoint) {
        edges.push_back(std::make_pair(curPoint, nextPoint));
        curPoint = nextPoint;
      }
    }
    thisSide = sides.next(thisSide);
  }
  edges.push_back(std::make_pair(curPoint, endPoint));
}

//------------------------------------------------------------------------------
// Close an ordered loop of clipped edges.  Consecutive edges either meet at a
// node or are joined by the CCW box path between their clipped endpoints.
//------------------------------------------------------------------------------
inline std::vector<Edge>
closeClippedEdges(const std::vector<ClippedEdge>& clippedEdges,
                  const std::map<BoxSide, unsigned>& cornerIndices) {
  const auto N = clippedEdges.size();
  POLY_ASSERT2(N > 0, "Must have at least 1 edge");
  std::vector<Edge> out;
  out.reserve(2 * N + 4);

  for (auto i = 0u; i < N; ++i) {
    const auto& cur = clippedEdges[i];
    const auto& next = clippedEdges[(i + 1) % N];
    const auto& curEdge = cur.curEdge;
    out.push_back(curEdge);

    if (curEdge.second == next.curEdge.first) continue;

    POLY_ASSERT2(cur.clippedSides.second >= 0 &&
                 next.clippedSides.first >= 0,
                 "Disconnected clipped edges without a box connection");
    walkBoxEdges(static_cast<BoxSide>(cur.clippedSides.second),
                 static_cast<BoxSide>(next.clippedSides.first),
                 curEdge.second, next.curEdge.first, cornerIndices, out);
  }
  return out;
}

//------------------------------------------------------------------------------
// Order any clipped edges based on the generator location. This must be done
// before calling orderClippedEdges
//------------------------------------------------------------------------------
template<typename CoordType>
void orientClippedEdges(std::vector<ClippedEdge>& clippedEdges,
                        const Point<2, CoordType>& generator,
                        const std::vector<Point<2, CoordType>>& nodes) {
  for (auto& clipped : clippedEdges) {
    const auto& p0 = nodes[clipped.curEdge.first];
    const auto& p1 = nodes[clipped.curEdge.second];
    if (qcross<CoordType>(p1 - p0, generator - p0) < 0) {
      std::swap(clipped.curEdge.first, clipped.curEdge.second);
      std::swap(clipped.clippedSides.first, clipped.clippedSides.second);
    }
  }
}

//------------------------------------------------------------------------------
// Order a loop of edges to form a connected chain
// Ensures edges[i][1] connects to edges[i+1][0] when possible and relies
// on clipping-box sides when not possible. Assumes 2D.
//------------------------------------------------------------------------------
inline void orderClippedEdges(std::vector<ClippedEdge>& clippedEdges) {
  if (clippedEdges.empty()) return;
  // Remove any degenerate edges
  clippedEdges.erase(
    std::remove_if(clippedEdges.begin(), clippedEdges.end(),
                   [](const ClippedEdge& edge) {
                     return edge.curEdge.first == edge.curEdge.second;
                   }),
    clippedEdges.end());
  if (clippedEdges.empty()) return;
  const auto N = clippedEdges.size();
  auto sideDistance = [](int fromSide, int toSide) {
                        constexpr int numBoxSides = 8;
                        return (toSide - fromSide + numBoxSides) % numBoxSides;
                      };
  std::vector<ClippedEdge> orderedEdges;
  orderedEdges.reserve(N);
  std::vector<bool> used(N, false);
  size_t current = 0;

  for (size_t count = 0; count < N; ++count) {
    orderedEdges.push_back(clippedEdges[current]);
    used[current] = true;
    if (count + 1 == N) break;
    int next = -1;
    // Prefer exact Voronoi-edge adjacency.
    for (size_t candidate = 0; candidate < N; ++candidate) {
      if (!used[candidate] &&
          clippedEdges[current].curEdge.second ==
          clippedEdges[candidate].curEdge.first) {
        next = static_cast<int>(candidate);
        break;
      }
    }
    // If the current edge ends on the clipping box, connect to the next
    // edge that starts on the box by walking the box in CCW side order.
    if (next == -1 && clippedEdges[current].clippedSides.second >= 0) {
      int bestDistance = 8;
      for (size_t candidate = 0; candidate < N; ++candidate) {
        if (!used[candidate] &&
            clippedEdges[candidate].clippedSides.first >= 0) {
          const int distance = sideDistance(
              clippedEdges[current].clippedSides.second,
              clippedEdges[candidate].clippedSides.first);
          if (distance < bestDistance) {
            bestDistance = distance;
            next = static_cast<int>(candidate);
          }
        }
      }
    }
    POLY_ASSERT2(next != -1,
                 "Unable to order clipped edges by node or CCW box adjacency");
    current = static_cast<size_t>(next);
  }
  clippedEdges = std::move(orderedEdges);
}

inline void orderEdgeLoop(std::vector<std::vector<unsigned>>& edges) {
  if (edges.empty()) return;

  std::vector<std::vector<unsigned>> ordered;
  ordered.reserve(edges.size());

  // Build map: start vertex -> edge index
  std::map<int, int> startMap;
  for (size_t i = 0; i < edges.size(); ++i) {
    startMap[edges[i][0]] = i;
  }

  // Follow the chain starting from first edge
  std::set<int> used;
  int current = 0;
  while (used.size() < edges.size()) {
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

//------------------------------------------------------------------------------
// Edge storage and orientation tracking
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Utilities for edge data, meaning edges paired with generator points
// This allows us to keep track of which edges belong to which generators
//------------------------------------------------------------------------------
using GenPair = std::pair<int, int>;
inline GenPair orderGenPair(const int a, const int b) {
  return orderEdge(a, b);
}

inline ClippedEdge flipEdge(const ClippedEdge& clippedEdge) {
  ClippedEdge out(clippedEdge);
  std::swap(out.curEdge.first, out.curEdge.second);
  std::swap(out.clippedSides.first, out.clippedSides.second);
  return out;
}

using GenPairToClippedEdgeMap = std::map<GenPair, ClippedEdge>;

//------------------------------------------------------------------------------
// Add an oriented edge to the edge map
// Returns the signed face index:
//   - Positive if edge orientation matches canonical form
//   - Negative (bitwise NOT) if edge orientation is reversed
//------------------------------------------------------------------------------
inline int addOrientedEdge(int n0, int n1,
                           std::vector<std::vector<unsigned>>& faces,
                           EdgeToFaceMap& edgeToFace) {
  Edge canonical = orderEdge(n0, n1);

  auto it = edgeToFace.find(canonical);
  int faceIndex;

  if (it == edgeToFace.end()) {
    // New edge - add to faces in canonical form
    faceIndex = faces.size();
    edgeToFace[canonical] = faceIndex;
    faces.push_back({static_cast<unsigned>(canonical.first),
                     static_cast<unsigned>(canonical.second)});
  } else {
    faceIndex = it->second;
  }

  // Return signed index based on whether orientation matches canonical
  return (canonical.first == n0) ? faceIndex : ~faceIndex;
}

//------------------------------------------------------------------------------
// Reverse the order of the edges, both the order of each edge and
// the order of the edges
//------------------------------------------------------------------------------
inline void reverseEdgeLoop(std::vector<Edge>& edges) {
  // Reverse every directed edge.
  for (auto& edge : edges) {
    std::swap(edge.first, edge.second);
  }
  // Reverse traversal order so the edges remain a connected loop.
  std::reverse(edges.begin(), edges.end());
}

//------------------------------------------------------------------------------
// Assemble an unordered set of clipped Voronoi edges into one CCW cell.
// Each edge is first directed with the owning generator on its left.  The
// routine then follows direct node connections, filling only box-boundary gaps
// with CCW box edges, before creating the signed face references.
//------------------------------------------------------------------------------
template<typename CoordType>
inline std::vector<int>
makeCCWCellFromClippedEdges(std::vector<ClippedEdge> clippedEdges,
                            const Point<2, CoordType>& generator,
                            const std::map<BoxSide, unsigned>& cornerIndices,
                            const std::vector<Point<2, CoordType>>& nodes,
                            std::vector<std::vector<unsigned>>& faces,
                            EdgeToFaceMap& edgeToFace) {
  POLY_ASSERT2(!clippedEdges.empty(), "Cannot construct a cell without edges");

  // Direct each Voronoi edge so the cell interior is on its left.
  orientClippedEdges(clippedEdges, generator, nodes);

  // Establish edge order, then fill every non-node-connected transition with
  // its CCW clipping-box path.
  orderClippedEdges(clippedEdges);
  POLY_ASSERT2(!clippedEdges.empty(), "All clipped edges collapsed to zero length");
  auto edges = closeClippedEdges(clippedEdges, cornerIndices);

  POLY_ASSERT2(edges.size() >= 3, "Degenerate cell after clipping");

  // Store each canonical face once; retain per-cell direction in its
  // signed face index.
  std::vector<int> cell;
  cell.reserve(edges.size());
  for (const auto& edge : edges) {
    cell.push_back(addOrientedEdge(edge.first, edge.second, faces, edgeToFace));
  }

  return cell;
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
                                            const std::vector<std::vector<unsigned>>& faces) {
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
// Modify the nodes list if points do not exist in a given node id map
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
inline Edge updateNodeMap(const Point<Dimension, CoordType>& p0,
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
  return Edge(std::make_pair(n0, n1));
}

}
#endif
