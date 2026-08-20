#ifndef __Polytope_Intersections__
#define __Polytope_Intersections__

//------------------------------------------------------------------------------
// Intersection routines with integer overflow protection
//
// These routines work with quantized (integer) coordinates and prevent overflow
// by casting to larger integer types before performing multiplications.
//
// Type system:
//   - CoordType: Input coordinate type (int for 2D, int64_t for 3D)
//   - WideInt: Larger type for intermediate calculations (int64_t for 2D, __int128 for 3D)
//   - WideInt is automatically deduced from dimension via MortonKeyTraits<Dim>
//
// Strategy:
//   - Cast to WideInt before multiplication to prevent overflow
//   - Delay divisions by carrying denominators through calculations
//------------------------------------------------------------------------------
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>

#include "GeomUtils.hh"
#include "Cell.hh"
#include "Quantizer.hh"
#include "Shapes.hh"

#ifdef POLYTOPE_ENABLE_BOOST
#include "RegisterBoostPolygonTypes.hh"
#endif

namespace polytope {

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Convex hull intersection methods
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Test if points are contained in a convex hull using half-space intersection
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolygon_convex(const std::vector<std::vector<unsigned>>& facets,
                           const std::vector<Point2<CoordType>>& vertices,
                           const std::vector<Point2<CoordType>>& points) {
  for (const auto& f : facets) {
    auto vi = vertices[f[0]];
    auto vj = vertices[f[1]];
    if (aboveBelow(vi, vj, points)) {
      return true;
    }
  }
  return false;
}

template<typename CoordType>
bool pointInPolygon_convex(const typename Cell<2, CoordType>::CellType& vertices,
                           const std::vector<Point2<CoordType>>& points) {
  auto N = vertices.size();
  for (int i = 0; i < N; ++i) {
    auto vi = vertices[i];
    auto vj = vertices[(i+1)%N];
    if (aboveBelow(vi, vj, points)) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Tests if two convex hulls overlap.
//------------------------------------------------------------------------------
template<typename CoordType>
bool convexIntersect(const typename Cell<2, CoordType>::CellType& pointsA,
                     const typename Cell<2, CoordType>::CellType& pointsB) {
  if (pointsA.empty() || pointsB.empty()) return false;

  using Wide = WideInt<2>;
  auto separated = [&pointsA, &pointsB](const typename Cell<2, CoordType>::CellType& axesFrom) {
    const auto N = axesFrom.size();
    for (auto i = 0u; i < N; ++i) {
      const auto& pi = axesFrom[i];
      const auto& pj = axesFrom[(i + 1)%N];
      const Wide dx = static_cast<Wide>(pj.x) - static_cast<Wide>(pi.x);
      const Wide dy = static_cast<Wide>(pj.y) - static_cast<Wide>(pi.y);
      if (dx == 0 && dy == 0) continue;
      if (SAT(pointsA, pointsB, Point2<Wide>(-dy, dx))) {
        return true;
      }
    }
    return false;
  };

  if (separated(pointsA)) return false;
  if (separated(pointsB)) return false;
  return true;
}

template<typename CoordType>
bool convexIntersect(const std::vector<Point2<CoordType>>& pointsA,
                     const std::vector<std::vector<unsigned>>& facetsA,
                     const std::vector<Point2<CoordType>>& pointsB,
                     const std::vector<std::vector<unsigned>>& facetsB) {
  const auto cellA = Cell<2, CoordType>::extractCell(pointsA, facetsA);
  const auto cellB = Cell<2, CoordType>::extractCell(pointsB, facetsB);
  return convexIntersect<CoordType>(cellA, cellB);
}

template<typename CoordType>
bool convexIntersect(const std::vector<Point3<CoordType>>& /*pointsA*/,
                     const std::vector<std::vector<unsigned>>& /*facetsA*/,
                     const std::vector<Point3<CoordType>>& /*normalsA*/,
                     const std::vector<Point3<CoordType>>& /*pointsB*/,
                     const std::vector<std::vector<unsigned>>& /*facetsB*/,
                     const std::vector<Point3<CoordType>>& /*normalsB*/) {
  // Implement this
  return false;
}

//------------------------------------------------------------------------------
// Tests if the boundaries of two polygons intersect. Only returns true if any
// edges intersect but not if one completely contains the other.
//------------------------------------------------------------------------------

template<typename CoordType>
bool convexBoundaryIntersect(const typename Cell<2, CoordType>::CellType& pointsA,
                             const typename Cell<2, CoordType>::CellType& pointsB) {
  Point2<CoordType> result;
  auto Na = pointsA.size();
  auto Nb = pointsB.size();
  for (auto fa = 0u; fa < Na; ++fa) {
    auto vai = pointsA[fa];
    auto vaj = pointsA[(fa+1)%Na];
    for (auto fb = 0u; fb < Nb; ++fb) {
      auto vbi = pointsB[fb];
      auto vbj = pointsB[(fb+1)%Nb];
      if (segmentIntersection2D(vai, vaj, vbi, vbj, result)) {
        return true;
      }
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Tests if the boundaries of two polyhedra intersect but not if one is contained
// in the other.
//------------------------------------------------------------------------------
template<typename CoordType>
bool convexBoundaryIntersect(const typename Cell<3, CoordType>::CellType& /*pointsA*/,
                             const typename Cell<3, CoordType>::CellType& /*pointsB*/) {
  // Implement this
  return false;
}

//------------------------------------------------------------------------------
// Point-in-polyhedron test using half-space intersection
//
// Tests if ANY of the test points is inside the convex polyhedron.
// A point is inside if it's on the negative (inside) side of ALL face planes.
// Assumes facets are oriented with outward-pointing normals.
// Uses precomputed normals for efficiency.
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolyhedron_convex(const std::vector<std::vector<unsigned>>& facets,
                              const std::vector<Point3<CoordType>>& vertices,
                              const std::vector<Point3<CoordType>>& normals,
                              const std::vector<Point3<CoordType>>& points) {
  // For each facet, use precomputed normal and test all points
  for (unsigned ifacet = 0; ifacet < facets.size(); ++ifacet) {
    const auto& v0 = vertices[facets[ifacet][0]];
    const auto& normal = normals[ifacet];

    if (aboveBelow(v0, normal, points)) {
      return true;
    }
  }

  // If we get here, at least one point is inside (not marked outside)
  return false;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// General (potentially non-convex) intersection methods
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef POLYTOPE_ENABLE_BOOST
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Boost Polygon methods
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace bp = boost::polygon;
using namespace boost::polygon::operators;

using Polygon = bp::polygon_data<QuantizedCoordinate<2>>;
using PolygonWithHoles = bp::polygon_with_holes_data<QuantizedCoordinate<2>>;
using PolygonSet = bp::polygon_set_data<QuantizedCoordinate<2>>;


inline std::vector<PolygonWithHoles>
boostUnion(const PolygonWithHoles& p1,
           const PolygonWithHoles& p2) {
  std::vector<PolygonWithHoles> intersect;
  bp::assign(intersect, p1 | p2);
  return intersect;
}

inline std::vector<PolygonWithHoles>
boostIntersect(const PolygonWithHoles& p1,
               const PolygonWithHoles& p2) {
  std::vector<PolygonWithHoles> out;
  bp::assign(out, p1 & p2);
  return out;
}

// Clip a Polytope cell against a boost Polygon and return a vector of polygons
inline std::vector<PolygonWithHoles>
boostIntersect(const Cell<2, QuantizedCoordinate<2>>::CellType& pcell,
               const PolygonWithHoles& boundary) {
  PolygonWithHoles cell = bp::polytopeToBoost(pcell);
  std::vector<PolygonWithHoles> out;
  bp::assign(out, cell & boundary);
  return out;
}

// Check if two polygons can be properly intersected.
// If they can, overwrite outPoly with the intersection and return true
inline bool
validUnion(const PolygonWithHoles& inPoly,
           PolygonWithHoles& outPoly) {
  auto trialUnion = boostUnion(outPoly, inPoly);
  if (trialUnion.size() == 1) {
    outPoly = trialUnion[0];
    return true;
  }
  return false;
}
#endif // POLYTOPE_ENABLE_BOOST

//------------------------------------------------------------------------------
// Remove collinear points and combine edges where necessary
//------------------------------------------------------------------------------
template<typename CoordType>
void removeCollinear(std::vector<Point2<CoordType>>& vertices) {
  const auto N = vertices.size();
  std::vector<Point2<CoordType>> result;
  for (unsigned i = 0; i < N; ++i) {
    auto p1 = vertices[(i+N-1)%N];
    auto p2 = vertices[i];
    auto p3 = vertices[(i+1)%N];
    if (!collinear(p1, p3, p2)) {
      result.push_back(p2);
    }
  }
  vertices = std::move(result);
}

//------------------------------------------------------------------------------
// Remove collinear points from an edge loop
// Edges should form a connected chain where edges[i][1] == edges[i+1][0]
//------------------------------------------------------------------------------
template<typename CoordType>
void removeCollinear(std::vector<edge::Edge>& edges,
                     const std::vector<Point2<CoordType>>& vertices) {
  const auto N = edges.size();
  if (N < 3) return;  // Need at least 3 edges to have collinear points

  std::vector<edge::Edge> result;
  result.reserve(N);

  for (unsigned i = 0; i < N; ++i) {
    // Get three consecutive vertices from the edge loop
    auto p1_idx = edges[(i+N-1)%N].first;   // Previous edge start
    auto p2_idx = edges[i].first;           // Current edge start
    auto p3_idx = edges[i].second;          // Current edge end

    const auto& p1 = vertices[p1_idx];
    const auto& p2 = vertices[p2_idx];
    const auto& p3 = vertices[p3_idx];

    // Keep this edge only if the middle point is not collinear
    if (!collinear(p1, p3, p2)) {
      result.push_back(edges[i]);
    } else {
      // Collinear: skip this edge and update the next edge to span across
      // The next edge should start from p1 instead of p2
      if (!result.empty()) {
        // Update the last added edge to connect to p3
        result.back().second = p3_idx;
      } else {
        // Special case: first edges are collinear, need to fix at the end
        // We'll handle this by checking if the first and last edges merge
      }
    }
  }

  // Handle wraparound: if we removed the first edge, the last edge might need adjustment
  if (result.size() < N && !result.empty()) {
    // Check if first and last result edges are now collinear
    auto p1_idx = result.back().first;
    auto p2_idx = result.back().second;
    auto p3_idx = result.front().second;

    const auto& p1 = vertices[p1_idx];
    const auto& p2 = vertices[p2_idx];
    const auto& p3 = vertices[p3_idx];

    if (collinear(p1, p3, p2)) {
      // Merge last and first edges
      result.back().second = result.front().second;
      result.erase(result.begin());
    }
  }

  edges = std::move(result);
}

//------------------------------------------------------------------------------
// Clip infinite ray provided by tessellator.
//
// Given the start of the ray and the direction, determine the boundary
// intersection point and the box side being intersected.
//------------------------------------------------------------------------------

template<typename CoordType>
bool
clipInfiniteRay(const Point2<CoordType>& validVertex,
                const Point2<CoordType>& normdiffg,
                Point2<CoordType>& result,
                BoxSide& side) {
  auto& Q = Quantizer<2>::instance();
  CoordType x_lim = (normdiffg.x > 0) ? Q.maxBound.x : Q.minBound.x;
  CoordType y_lim = (normdiffg.y > 0) ? Q.maxBound.y : Q.minBound.y;
  BoxSide LR = (normdiffg.x > 0) ? BoxSide::R : BoxSide::L;
  BoxSide TB = (normdiffg.y > 0) ? BoxSide::T : BoxSide::B;
  Point2<CoordType> planey1(Q.minBound.x, y_lim);
  Point2<CoordType> planey2(Q.maxBound.x, y_lim);
  Point2<CoordType> intersectionx, intersectiony;
  bool yint = segmentRayIntersection2D(planey1, planey2, validVertex, normdiffg, intersectiony);
  Point2<CoordType> planex1(x_lim, Q.minBound.y);
  Point2<CoordType> planex2(x_lim, Q.maxBound.y);
  bool xint = segmentRayIntersection2D(planex1, planex2, validVertex, normdiffg, intersectionx);
  if (!xint && !yint) {
    return false;
  }
  bool hitX = true;
  if (xint && yint) {
    if (intersectionx == intersectiony) {
      result = intersectionx;
      side = getBoxCorner(LR, TB);
      return true;
    }
    // Assume it intersects both planes, check if ||p-x|| is longer than ||p-y||
    if (magComparison(validVertex - intersectionx, validVertex - intersectiony)) {
      hitX = false;
    }
  } else if (yint) {
    hitX = false;
  }
  if (hitX) {
    result = intersectionx;
    side = LR;
  } else {
    result = intersectiony;
    side = TB;
  }
  return true;
}

//------------------------------------------------------------------------------
// Point-in-polygon test using ray casting algorithm
//
// Casts a horizontal ray from the point to the right and counts edge crossings.
// Odd count = inside, even count = outside.
//
// Handles vertex coincidence by shifting test point up slightly when it lies
// exactly on a horizontal edge endpoint.
// Since 2D routines are also used by 3D routines, the types are more variable
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolygon(const Point2<CoordType>& point,
                    const std::vector<Point2<CoordType>>& vertices) {
  using Wide = WideInt<2>;
  const auto N = vertices.size();
  bool inside = false;

  // Cast test point once
  const Wide px = static_cast<Wide>(point.x);
  Wide py = static_cast<Wide>(point.y);

  for (size_t i = 0; i < N; ++i) {
    // Get edge vertices (cast to Wide for all operations)
    auto vi = vertices[i].template type_cast<Wide>();
    auto vj = vertices[(i+1)%N].template type_cast<Wide>();

    // Ensure vi.y <= vj.y
    if (vi.y > vj.y) std::swap(vi, vj);

    // Shift test point if coincident with vertex to avoid double-counting
    Wide py_test = py;
    if (py_test == vi.y || py_test == vj.y) py_test++;

    // Check if ray at py_test height could intersect this edge
    if (py_test <= vi.y || py_test > vj.y) continue;

    // Check if point is definitely to the right of the edge
    if (px >= std::max(vi.x, vj.x)) continue;

    // Check if point is definitely to the left of the edge
    if (px < std::min(vi.x, vj.x)) {
      inside = !inside;
      continue;
    }

    // Point is horizontally between edge endpoints - compute exact crossing
    // Compare: px ? (vi.x + (py_test - vi.y) * (vj.x - vi.x) / (vj.y - vi.y))
    // Rearranged to avoid division: (px - vi.x) * (vj.y - vi.y) ? (py_test - vi.y) * (vj.x - vi.x)
    const Wide lhs = (px - vi.x) * (vj.y - vi.y);
    const Wide rhs = (py_test - vi.y) * (vj.x - vi.x);

    if (lhs < rhs) {
      inside = !inside;
    }
  }

  return inside;
}

//------------------------------------------------------------------------------
// Point-on-polygon boundary test
//
// Tests if a point lies exactly on any edge of the polygon.
// Returns true if the point is on the boundary (edge or vertex), false otherwise.
//
// Algorithm:
//   For each edge [vi, vj]:
//   1. Check collinearity: (point - vi) × (vj - vi) == 0
//   2. Check if point is between vi and vj using bounding box test
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointOnPolygon(const Point2<CoordType>& point,
                    const std::vector<Point2<CoordType>>& vertices) {
  const auto N = vertices.size();

  for (size_t i = 0; i < N; ++i) {
    // Get edge vertices
    const auto vi = vertices[i];
    const auto vj = vertices[(i + 1) % N];
    if (collinear(vi, vj, point)) {
      return true;
    }
  }
  return false;  // Point is not on any edge
}

//------------------------------------------------------------------------------
// Point-in-polygon test with indexed vertices
// Convenience wrapper that extracts vertices by index
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolygon(const Point2<CoordType>& point,
                    const std::vector<unsigned>& faceIndices,
                    const std::vector<Point2<CoordType>>& vertices) {
  const size_t N = faceIndices.size();
  std::vector<Point2<CoordType>> faceVerts;
  faceVerts.reserve(N);
  for (const auto idx : faceIndices) {
    faceVerts.push_back(vertices[idx]);
  }
  return pointInPolygon(point, faceVerts);
}

//------------------------------------------------------------------------------
// Point-on-polygon test with indexed vertices
// Convenience wrapper that extracts vertices by index
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointOnPolygon(const Point2<CoordType>& point,
                    const std::vector<unsigned>& faceIndices,
                    const std::vector<Point2<CoordType>>& vertices) {
  const size_t N = faceIndices.size();
  std::vector<Point2<CoordType>> faceVerts;
  faceVerts.reserve(N);
  for (const auto idx : faceIndices) {
    faceVerts.push_back(vertices[idx]);
  }
  return pointOnPolygon(point, faceVerts);
}

//------------------------------------------------------------------------------
// 2D line segment intersection
//
// General routine for finding intersection point
//------------------------------------------------------------------------------
template<typename CoordType>
bool intersection2D(const Point2<CoordType>& a,
                    const Point2<CoordType>& b,
                    const Point2<CoordType>& c,
                    const Point2<CoordType>& n,
                    Point2<CoordType>& result,
                    bool isRay = false) {
  // Integer implementation with overflow protection
  using Wide = WideInt<2>;
  using Big = BigInt<2>;

  const Point2<CoordType> s = b - a;
  const Point2<CoordType> ca = a - c;

  Wide denom = qcross<CoordType>(n, s);
  if (denom == 0) return false;
  Wide t_num = qcross<CoordType>(ca, s);
  Wide u_num = qcross<CoordType>(ca, n);
  if (denom < 0) {
    t_num = -t_num;
    u_num = -u_num;
    denom = -denom;
  }
  bool no_intersection = t_num < 0 || u_num < 0 || u_num > denom;
  if (!isRay) {
    no_intersection = no_intersection || t_num > denom;
  }
  if (no_intersection) return false;
  auto bn = n.template type_cast<Big>()*static_cast<Big>(t_num);
  auto cn = c.template type_cast<Big>();
  auto bd = static_cast<Big>(denom);
  result = ((bd*cn + bn)/bd).template type_cast<CoordType>();
  return true;
}

//------------------------------------------------------------------------------
// 2D line segment intersection
//
// Tests if segment [a, b] intersects segment [c, d].
// If they intersect, returns true and fills result with intersection point.
//
// Integer version: Uses parametric form and cross products to avoid overflow.
// Floating-point version: Uses standard floating-point arithmetic.
//
// Returns false if:
//   - Segments are parallel/collinear
//   - Intersection point is outside either segment's bounds
//------------------------------------------------------------------------------
template<typename CoordType>
bool segmentIntersection2D(const Point2<CoordType>& a,
                           const Point2<CoordType>& b,
                           const Point2<CoordType>& c,
                           const Point2<CoordType>& d,
                           Point2<CoordType>& result) {
  const Point2<CoordType> r = d - c;
  return intersection2D(a, b, c, r, result);
}

 //------------------------------------------------------------------------------
// 2D ray and line segment intersection
//
// Tests if segment [a, b] intersects ray that starts at c and travels in direction n.
// If they intersect, returns true and fills result with intersection point.
//
// Integer version: Uses parametric form and cross products to avoid overflow.
// Floating-point version: Uses standard floating-point arithmetic.
//
// Returns false if:
//   - Segments are parallel/collinear
//   - Intersection point is outside either segment's bounds
//------------------------------------------------------------------------------
template<typename CoordType>
bool segmentRayIntersection2D(const Point2<CoordType>& a,
                              const Point2<CoordType>& b,
                              const Point2<CoordType>& c,
                              const Point2<CoordType>& n,
                              Point2<CoordType>& result) {
  return intersection2D(a, b, c, n, result, true);
}
 
//------------------------------------------------------------------------------
// 3D segment-plane intersection
//------------------------------------------------------------------------------
template<typename CoordType>
int segmentPlaneIntersection3D(const Point3<CoordType>& segStart,
                               const Point3<CoordType>& segEnd,
                               const Point3<CoordType>& v0,
                               const Point3<CoordType>& plane_normal,
                               Point3<WideInt<3>>& result,
                               WideInt<3>& denom) {
  // Implement this
  return 0;
}

//------------------------------------------------------------------------------
// 3D segment-face intersection
//
// Tests if segment [segStart, segEnd] intersects a polygonal face.
//
// Algorithm:
//   1. Use precomputed face plane normal
//   2. Test segment-plane intersection
//   3. Project intersection point and face to 2D
//   4. Test if 2D point is inside 2D polygon
//
// Returns -1 if no intersection exists.
// Returns 0 if ray is coplanar.
// Otherwise 1 and fills result with intersection point.
//------------------------------------------------------------------------------
template<typename CoordType>
int segmentFaceIntersection3D(const Point3<CoordType>& segStart,
                              const Point3<CoordType>& segEnd,
                              const std::vector<unsigned>& faceIndices,
                              const std::vector<Point3<CoordType>>& vertices,
                              const Point3<CoordType>& plane_normal,
                              Point3<CoordType>& result) {
  // Implement this
  return 0;
}

} // namespace polytope

#endif
