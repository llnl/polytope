#ifndef __Polytope_Shapes__
#define __Polytope_Shapes__

//------------------------------------------------------------------------------
// Simple 2D box routines.
//------------------------------------------------------------------------------

#include "Point.hh"
#include "QuantizedKeyTraits.hh"

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
addBoxPoints(std::map<Point<2, CoordType>, int>& node2id,
             std::vector<Point<2, CoordType>>& nodes) {
  const auto& Q = Quantizer<2>::instance();
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

}
#endif
