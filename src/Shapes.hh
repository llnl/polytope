#ifndef __Polytope_Shapes__
#define __Polytope_Shapes__

//------------------------------------------------------------------------------
// Simple 2D box routines.
//------------------------------------------------------------------------------

#include "Point.hh"
#include "MortonKeyTraits.hh"
#include "EdgeUtils.hh"

namespace polytope {

enum BoxSide {
    L, // Left
    LL,
    B, // Bottom
    LR,
    R, // Right
    UR,
    T, // Top
    UL
};

struct BoxSides {
  std::array<BoxSide, 8> sides;
  std::array<BoxSide, 4> corners;
  BoxSides() {
    for (int i = 0; i < 8; ++i) {
      sides[i] = static_cast<BoxSide>(i);
    }
    for (int i = 0; i < 4; ++i) {
      corners[i] = static_cast<BoxSide>(2*i+1);
    }
  }
  BoxSide next(const BoxSide& side) {
    auto i = static_cast<int>(side);
    return sides[(i+1)%8];
  }
  BoxSide prev(const BoxSide& side) {
    auto i = static_cast<int>(side);
    return sides[(i+8-1)%8];
  }
  BoxSide corner(const int i) {
    return corners[i];
  }
};

inline bool isCorner(const BoxSide& side) {
  return (static_cast<int>(side)%2 == 1);
}

// 2D specialization with explicit CCW ordering
template<typename CoordType>
inline std::vector<Point2<CoordType>>
createSquarePoints(const Point2<CoordType>& min,
                   const Point2<CoordType>& max) {
  std::vector<Point2<CoordType>> out;
  out.reserve(4);

  // Explicit CCW order: LL → LR → UR → UL
  out.push_back(Point2<CoordType>(min.x, min.y)); // LL (index 0)
  out.push_back(Point2<CoordType>(max.x, min.y)); // LR (index 1)
  out.push_back(Point2<CoordType>(max.x, max.y)); // UR (index 2)
  out.push_back(Point2<CoordType>(min.x, max.y)); // UL (index 3)

  for (unsigned i = 0; i < 4; ++i) {
    out[i].index = i;
  }
  return out;
}

inline std::vector<std::vector<unsigned>> createSquareFaces() {
  std::vector<std::vector<unsigned>> coords(4);
  for (auto f = 0; f < 4; ++f) {
    coords[f].resize(2);
    coords[f][0] = f;
    coords[f][1] = (f+1)%4;
  }
  return coords;
}

// Return the corner based on two sides
inline BoxSide getBoxCorner(const BoxSide& s1, const BoxSide& s2) {
  if ((s1 == BoxSide::L && s2 == BoxSide::B) ||
      (s1 == BoxSide::B && s2 == BoxSide::L)) {
    return BoxSide::LL;
  } else if ((s1 == BoxSide::R && s2 == BoxSide::B) ||
             (s1 == BoxSide::B && s2 == BoxSide::R)) {
    return BoxSide::LR;
  } else if ((s1 == BoxSide::R && s2 == BoxSide::T) ||
             (s1 == BoxSide::T && s2 == BoxSide::R)) {
    return BoxSide::UR;
  } else if ((s1 == BoxSide::L && s2 == BoxSide::T) ||
             (s1 == BoxSide::T && s2 == BoxSide::L)) {
    return BoxSide::UL;
  }
  return BoxSide::LL; // Should never get here
}

template<typename CoordType>
inline std::map<BoxSide, unsigned>
addBoxPoints(const Quantizer<2>& Q,
             std::map<Point<2, CoordType>, int>& node2id,
             std::vector<Point<2, CoordType>>& nodes) {
  std::map<BoxSide, unsigned> cornerIndices; // Ordered lower left and CCW
  std::vector<Point<2, CoordType>> box = createSquarePoints(Q.minBound, Q.maxBound);
  BoxSides sides;
  for (unsigned i = 0; i < 4; i++) {
    const auto n = nodes.size();
    cornerIndices[sides.corner(i)] = n;
    node2id[box[i]] = n;
    nodes.push_back(box[i]);
  }
  return cornerIndices;
}

// Walk box edges only in CCW direction
inline void walkBoxEdges(const BoxSide& startSide,
                         const BoxSide& endSide,
                         const unsigned& startPoint,
                         const unsigned& endPoint,
                         const std::map<BoxSide, unsigned>& cornerIndices,
                         std::vector<edge::Edge>& edges) {
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

// Close any clipped edges
inline std::vector<edge::Edge> closeClippedEdges(const std::vector<edge::Edge>& origEdges,
                                                 const std::vector<std::pair<int, int>>& clippedNodeSides,
                                                 const std::map<BoxSide, unsigned>& cornerIndices) {
  auto N = origEdges.size();
  POLY_ASSERT2(N > 0, "Must have at least 1 edge");
  std::vector<edge::Edge> out;
  out.reserve(N);
  if (N == 1) {
    BoxSide endSide = static_cast<BoxSide>(clippedNodeSides[0].first);
    BoxSide startSide = static_cast<BoxSide>(clippedNodeSides[0].second);
    out.push_back(origEdges[0]);
    walkBoxEdges(startSide, endSide, origEdges[0].second, origEdges[0].first, cornerIndices, out);
    return out;
  }
  for (auto i = 0u; i < N; ++i) {
    auto curEdge = origEdges[i];
    auto side1 = clippedNodeSides[i].second;
    out.push_back(curEdge);
    if (side1 >= 0) {
      auto ip = (i+1)%N;
      auto nextEdge = origEdges[ip];
      BoxSide endSide = static_cast<BoxSide>(clippedNodeSides[ip].first);
      BoxSide startSide = static_cast<BoxSide>(side1);
      walkBoxEdges(startSide, endSide, curEdge.second, nextEdge.first, cornerIndices, out);
    }
  }
  return out;
}
}
#endif
