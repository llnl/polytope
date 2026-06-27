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
}

// Walk box edges from origSide to curSide, adding corner vertices along the way
//
// When clipping against a box boundary, cells may need to traverse multiple
// box sides. This function adds edges that walk along the box perimeter from
// one side to another, visiting corners in counter-clockwise order.
//
// Parameters:
//   origSide, curSide: Box sides to walk between (L, B, R, or T - NOT corners)
//   cornerIndices: Map from corner BoxSide (LL, LR, UR, UL) to vertex indices
//   curEdge: The clipped edge that intersects the boundary
//   origPoint: Starting vertex index
//   endAtCurEdge: Walk direction
//     - true: Walk origSide -> curSide, then add curEdge
//     - false: Add curEdge first, then walk curSide -> origSide
//   edges: Output edge list
//
// Example (endAtCurEdge=true):
//   origSide=L, curSide=B, origPoint=3 (UL), curEdge=(1,4) (LR->mid)
//   Walks L -> LL (add 3->0) -> B (add 0->1) -> adds curEdge (1->4)
//   Result: [(3,0), (0,1), (1,4)]
inline void walkBoxEdges(const BoxSide& origSide,
                         const BoxSide& curSide,
                         const std::map<BoxSide, unsigned>& cornerIndices,
                         const edge::Edge& curEdge,
                         const unsigned& origPoint,
                         const bool endAtCurEdge,
                         std::vector<edge::Edge>& edges) {
  BoxSides sides;

  if (endAtCurEdge) {
    // Walk origSide -> curSide, visiting corners, then add curEdge
    BoxSide thisSide = origSide;
    BoxSide endSide = curSide;
    unsigned curIndx = origPoint;

    while (thisSide != endSide) {
      if (isCorner(thisSide)) {
        unsigned nextIndx = cornerIndices.at(thisSide);
        if (curIndx != nextIndx) {
          edges.push_back(std::make_pair(curIndx, nextIndx));
          curIndx = nextIndx;
        }
      }
      thisSide = sides.next(thisSide);
    }

    // Connect to curEdge start if not already there
    if (curIndx != curEdge.first) {
      edges.push_back(std::make_pair(curIndx, curEdge.first));
    }
    edges.push_back(curEdge);

  } else {
    // Add curEdge first, then walk curSide -> origSide
    edges.push_back(curEdge);

    BoxSide thisSide = curSide;
    BoxSide endSide = origSide;
    unsigned curIndx = curEdge.second;

    while (thisSide != endSide) {
      if (isCorner(thisSide)) {
        unsigned nextIndx = cornerIndices.at(thisSide);
        if (curIndx != nextIndx) {
          edges.push_back(std::make_pair(curIndx, nextIndx));
          curIndx = nextIndx;
        }
      }
      thisSide = sides.next(thisSide);
    }

    // Close back to origPoint if needed
    if (curIndx != origPoint) {
      edges.push_back(std::make_pair(curIndx, origPoint));
    }
  }
}
}
}
#endif
