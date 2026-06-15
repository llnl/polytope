#ifndef POLYTOPE_SHAPES_HH
#define POLYTOPE_SHAPES_HH

//------------------------------------------------------------------------------
// Simple shape routines for 2D and 3D.
//------------------------------------------------------------------------------

#include "Point.hh"
#include "HashKey.hh"
#include "EdgeUtils.hh"

namespace polytope {
namespace shapes {

enum BoxCorner {
    LL, // Lower left
    LR, // Lower right
    UR, // Upper right
    UL  // Upper left
};

enum BoxSide {
    L, // Left
    B, // Bottom
    R, // Right
    T, // Top
};

struct BoxSides {
  std::array<BoxSide, 4> sides;
  BoxSides() {
    for (int i = 0; i < 4; ++i) {
      sides[i] = static_cast<BoxSide>(i);
    }
  }
  BoxSide next(const BoxSide& side) {
    auto i = static_cast<int>(side);
    return sides[(i+1)%4];
  }
  BoxSide prev(const BoxSide& side) {
    auto i = static_cast<int>(side);
    return sides[(i+4-1)%4];
  }
  static BoxSide opposite(const BoxSide& side) {
    switch (side) {
    case BoxSide::T:
      return BoxSide::B;
    case BoxSide::B:
      return BoxSide::T;
    case BoxSide::L:
      return BoxSide::R;
    case BoxSide::R:
      return BoxSide::L;
    }
  }
};

// 2D specialization with explicit CCW ordering
template<typename CoordType>
inline std::vector<Point2<CoordType>>
createBoxPoints(const Point2<CoordType>& min,
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

// 3D version using bit pattern (no specific ordering required yet)
template<typename CoordType>
inline std::vector<Point3<CoordType>>
createBoxPoints(const Point3<CoordType>& min,
                const Point3<CoordType>& max) {
  std::vector<Point3<CoordType>> out;
  int count = 8;
  unsigned kk = 0;
  for (int i = 0; i < count; ++i) {
    Point3<CoordType> corner;
    for (int d = 0; d < 3; ++d) {
      corner[d] = min[d] + max[d]*((i >> d) & 1);
    }
    corner.index = kk++;
    out.push_back(corner);
  }
  return out;
}

template<typename CoordType>
inline std::vector<Point2<CoordType>>
createSquarePoints(const Point2<CoordType>& min,
                   const Point2<CoordType>& max) {
  return createBoxPoints<CoordType>(min, max);
}

template<typename CoordType>
inline std::vector<Point3<CoordType>>
createCubePoints(const Point3<CoordType>& min,
                 const Point3<CoordType>& max) {
  return createBoxPoints<CoordType>(min, max);
}

inline std::vector<std::vector<int>> createCubeFaces() {
  // Faces with outward-pointing normals (right-hand rule)
  return {
          {0, 1, 2, 3},  // bottom (-z)
          {4, 7, 6, 5},  // top (+z)
          {0, 4, 5, 1},  // front (-y)
          {2, 6, 7, 3},  // back (+y)
          {0, 3, 7, 4},  // left (-x)
          {1, 5, 6, 2}   // right (+x)
  };
}

inline std::vector<std::vector<int>> createSquareFaces() {
  std::vector<std::vector<int>> coords(4);
  for (auto f = 0; f < 4; ++f) {
    coords[f].resize(2);
    coords[f][0] = f;
    coords[f][1] = (f+1)%4;
  }
  return coords;
}

// Return the corner based on two sides
inline BoxCorner getBoxCorner(const BoxSide& s1, const BoxSide& s2) {
  if ((s1 == BoxSide::L && s2 == BoxSide::B) ||
      (s1 == BoxSide::B && s2 == BoxSide::L)) {
    return BoxCorner::LL;
  } else if ((s1 == BoxSide::R && s2 == BoxSide::B) ||
             (s1 == BoxSide::B && s2 == BoxSide::R)) {
    return BoxCorner::LR;
  } else if ((s1 == BoxSide::R && s2 == BoxSide::T) ||
             (s1 == BoxSide::T && s2 == BoxSide::R)) {
    return BoxCorner::UR;
  } else if ((s1 == BoxSide::L && s2 == BoxSide::T) ||
             (s1 == BoxSide::T && s2 == BoxSide::L)) {
    return BoxCorner::UL;
  }
}

// If curEdge is (inf, valid vertex), endAtCurEdge = true
inline void walkBoxEdges(const BoxSide& origSide,
                         const BoxSide& curSide,
                         const std::map<BoxCorner, unsigned>& cornerIndices,
                         const edge::Edge& curEdge,
                         const unsigned& origPoint,
                         const bool endAtCurEdge,
                         std::vector<edge::Edge>& edges) {
  BoxSide thisSide = (endAtCurEdge) ? origSide : curSide;
  BoxSide endSide = (endAtCurEdge) ? curSide : origSide;
  unsigned curIndx = (endAtCurEdge) ? origPoint : curEdge.second;
  unsigned nextIndx = (endAtCurEdge) ? curEdge.first : origPoint;
  BoxSides sides;
  if (endAtCurEdge) {
    while (thisSide != endSide) {
      auto nextSide = sides.next(thisSide);
      BoxCorner corner = getBoxCorner(thisSide, nextSide);
      nextIndx = cornerIndices.at(corner);
      if (curIndx != nextIndx) {
        edges.push_back(std::make_pair(curIndx, nextIndx));
        curIndx = nextIndx;
      }
      thisSide = nextSide;
    }
    edges.push_back(std::make_pair(curIndx, curEdge.first));
    edges.push_back(curEdge);
  } else {
    edges.push_back(curEdge);
    while (thisSide != endSide) {
      auto nextSide = sides.next(thisSide);
      BoxCorner corner = getBoxCorner(thisSide, nextSide);
      nextIndx = cornerIndices.at(corner);
      if (curIndx != nextIndx) {
        edges.push_back(std::make_pair(curIndx, nextIndx));
        curIndx = nextIndx;
      }
      thisSide = nextSide;
    }
    edges.push_back(std::make_pair(curIndx, origPoint));
  }
}
}
}
#endif
