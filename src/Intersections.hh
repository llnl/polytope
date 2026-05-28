#ifndef POLYTOPE_INTERSECTIONS_HH
#define POLYTOPE_INTERSECTIONS_HH

//------------------------------------------------------------------------------
// Intersection routines with integer overflow protection
//
// These routines work with quantized (integer) coordinates and prevent overflow
// by casting to larger integer types before performing multiplications.
//
// Type system:
//   - CoordType: Input coordinate type (int for 2D, int64_t for 3D)
//   - WideInt: Larger type for intermediate calculations (int64_t for 2D, __int128 for 3D)
//   - WideInt is automatically deduced from dimension via HashKey<Dim>::CoordHash
//     (named CoordHash in HashKey for historical reasons, but used here for arithmetic)
//
// Strategy:
//   - Cast to WideInt before multiplication to prevent overflow
//   - Delay divisions by carrying denominators through calculations
//   - Return results normalized to CoordType when possible
//------------------------------------------------------------------------------
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>

#include "GeomUtils.hh"

namespace polytope {

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Convex hull intersection methods
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Test if points intersect with a convex hull using half-space intersection, 2D
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolygon_convex(const std::vector<std::vector<int>>& facets,
                           const std::vector<Point2<CoordType>>& vertices,
                           const std::vector<Point2<CoordType>>& points) {
  const auto NF = facets.size();
  const auto NP = points.size();
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
bool convexIntersect(const std::vector<Point2<CoordType>>& pointsA,
                     const std::vector<std::vector<int>>& facetsA,
                     const std::vector<Point2<CoordType>>& pointsB,
                     const std::vector<std::vector<int>>& facetsB) {
  if (pointInPolygon_convex(facetsA, pointsA, pointsB) ||
      pointInPolygon_convex(facetsB, pointsB, pointsA)) {
    return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Point-in-polyhedron test using half-space intersection
//
// Tests if ANY of the test points is inside the convex polyhedron.
// A point is inside if it's on the negative (inside) side of ALL face planes.
// Assumes facets are oriented with outward-pointing normals.
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolyhedron_convex(const std::vector<std::vector<int>>& facets,
                              const std::vector<Point3<CoordType>>& vertices,
                              const std::vector<Point3<CoordType>>& points) {
  const auto NF = facets.size();
  const auto NP = points.size();

  // For each facet, compute normal once and test all points
  for (unsigned ifacet = 0; ifacet < NF; ++ifacet) {
    const auto& v0 = vertices[facets[ifacet][0]];
    const auto& v1 = vertices[facets[ifacet][1]];
    const auto& v2 = vertices[facets[ifacet][2]];

    // Compute the plane normal (outward-pointing)
    const auto normal = computeNormal(v0, v1, v2);

    if (aboveBelow(v0, normal, points)) {
      return true;
    }
  }

  // If we get here, at least one point is inside (not marked outside)
  return false;
}

template<typename CoordType>
bool convexIntersection(const std::vector<Point3<CoordType>>& pointsA,
                        const std::vector<std::vector<int>>& facetsA,
                        const std::vector<Point3<CoordType>>& pointsB,
                        const std::vector<std::vector<int>>& facetsB) {
  // Separating Axis Theorem: Test face normals of A as potential separating axes
  // For each face of A, check if all vertices of B are on the positive (outside) side
  for (const auto& facet : facetsA) {
    const auto& v0 = pointsA[facet[0]];
    const auto& v1 = pointsA[facet[1]];
    const auto& v2 = pointsA[facet[2]];
    const auto normal = computeNormal(v0, v1, v2);
    if (SAT(pointsA, pointsB, normal)) {
      return false;
    }
  }

  // Separating Axis Theorem: Test face normals of B as potential separating axes
  // For each face of B, check if all vertices of A are on the positive (outside) side
  for (const auto& facet : facetsB) {
    const auto& v0 = pointsB[facet[0]];
    const auto& v1 = pointsB[facet[1]];
    const auto& v2 = pointsB[facet[2]];
    const auto normal = computeNormal(v0, v1, v2);
    if (SAT(pointsA, pointsB, normal)) {
      return false;
    }
  }

  // Separating Axis Theorem: Test edge-edge cross products as potential separating axes
  std::vector<std::pair<int, int>> edgesA, edgesB;
  unsigned i, j;
  for (const auto& f : facetsA) {
    const auto n = f.size();
    for (i = 0; i < n; ++i) {
      j = (i+1)%n;
      edgesA.push_back({f[i], f[j]});
    }
  }
  for (const auto& f : facetsB) {
    const auto n = f.size();
    for (i = 0; i < n; ++i) {
      j = (i+1)%n;
      edgesB.push_back({f[i], f[j]});
    }
  }

  unsigned a0, a1, b0, b1;
  for (const auto& edgeAI : edgesA) {
    a0 = edgeAI.first;
    a1 = edgeAI.second;
    const auto& edgeA = pointsA[a1] - pointsA[a0];
    for (const auto& edgeBI : edgesB) {
      b0 = edgeBI.first;
      b1 = edgeBI.second;
      const auto& edgeB = pointsB[b1] - pointsB[b0];
      // Find the cross-product of the edges
      const auto axis = norm_qcross<CoordType>(edgeA, edgeB);
      if (axis.iszero()) {
        continue;
      }
      if (SAT(pointsA, pointsB, axis)) {
        return false;
      }
    }
  }
  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// General (potentially non-convex) intersection methods
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
template<typename CoordType, typename Wide>
bool pointInPolygon_impl(const Point2<CoordType>& point,
                         const std::vector<Point2<CoordType>>& vertices) {
  const auto N = vertices.size();
  bool inside = false;

  // Cast test point once
  const Wide px = static_cast<Wide>(point.x);
  Wide py = static_cast<Wide>(point.y);

  for (size_t i = 0; i < N; ++i) {
    // Get edge vertices (cast to Wide for all operations)
    auto vi = vertices[i].template type_cast<Wide>();
    auto vj = vertices[(i + 1) % N].template type_cast<Wide>();

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

template<typename CoordType>
bool pointInPolygon(const Point2<CoordType>& point,
                    const std::vector<Point2<CoordType>>& vertices) {
  if (typeid(CoordType) == typeid(NarrowInt<2>)) {
    return pointInPolygon_impl<CoordType, WideInt<2>>(point, vertices);
  } else {
    return pointInPolygon_impl<CoordType, WideInt<3>>(point, vertices);
  }
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
template<typename CoordType, typename Wide>
bool pointOnPolygon_impl(const Point2<CoordType>& point,
                         const std::vector<Point2<CoordType>>& vertices) {
  const auto N = vertices.size();

  // Cast test point once
  const auto ph = point.template type_cast<Wide>();

  for (size_t i = 0; i < N; ++i) {
    // Get edge vertices
    const auto vi = vertices[i].template type_cast<Wide>();
    const auto vj = vertices[(i + 1) % N].template type_cast<Wide>();

    // Vector from vi to point and from vi to vj
    const Point2<Wide> vp(ph.x - vi.x, ph.y - vi.y);
    const Point2<Wide> ve(vj.x - vi.x, vj.y - vi.y);

    // Check collinearity using cross product
    const Wide cross = vp.x * ve.y - vp.y * ve.x;
    if (cross != 0) continue;  // Not collinear with this edge

    // Point is collinear with edge - check if it's between endpoints
    // Use bounding box test for each coordinate
    const bool x_between = (ph.x >= std::min(vi.x, vj.x)) &&
                           (ph.x <= std::max(vi.x, vj.x));
    const bool y_between = (ph.y >= std::min(vi.y, vj.y)) &&
                           (ph.y <= std::max(vi.y, vj.y));

    if (x_between && y_between) {
      return true;  // Point is on this edge
    }
  }

  return false;  // Point is not on any edge
}

template<typename CoordType>
bool pointOnPolygon(const Point2<CoordType>& point,
                    const std::vector<Point2<CoordType>>& vertices) {
  if (typeid(CoordType) == typeid(NarrowInt<2>)) {
    return pointOnPolygon_impl<CoordType, WideInt<2>>(point, vertices);
  } else {
    return pointOnPolygon_impl<CoordType, WideInt<3>>(point, vertices);
  }
}

//------------------------------------------------------------------------------
// Point-in-polygon test with indexed vertices
// Convenience wrapper that extracts vertices by index
//------------------------------------------------------------------------------
template<typename CoordType>
bool pointInPolygon(const Point2<CoordType>& point,
                    const std::vector<int>& faceIndices,
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
                    const std::vector<int>& faceIndices,
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
// Tests if segment [a, b] intersects segment [c, d].
// If they intersect, returns true and fills result with intersection point.
//
// Uses parametric form and cross products to avoid floating-point operations.
// Carries normalization denominator to delay division until final result.
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
  using Wide = WideInt<2>;

  // Direction vectors
  const Point2<CoordType> r1 = b - a;  // Direction of segment [a, b]
  const Point2<CoordType> r2 = d - c;  // Direction of segment [c, d]
  const Point2<CoordType> ca = c - a;  // Vector from a to c

  // Compute denominator (r1 × r2)
  Wide denom = qcross<CoordType>(r1, r2);

  // Parallel or collinear segments
  if (denom == 0) return false;

  // Compute parametric coordinates
  // t = (ca × r2) / denom  (parameter along segment [a, b])
  // u = (ca × r1) / denom  (parameter along segment [c, d])
  Wide t_num = qcross<CoordType>(ca, r2);
  Wide u_num = qcross<CoordType>(ca, r1);

  // Normalize signs (make denominator positive)
  if (denom < 0) {
    t_num = -t_num;
    u_num = -u_num;
    denom = -denom;
  }

  // Check if intersection is within both segments [0, 1]
  if (t_num < 0 || t_num > denom || u_num < 0 || u_num > denom) {
    return false;
  }

  // Normalize denom and t_num to prevent overflow in final computation
  // Target: 21 bits to leave room for coordinate multiplication (31 + 21 = 52 < 64)
  auto shift = bitShiftAmount(denom, 21);
  auto denom_norm = denom >> shift;
  Wide t_num_norm = t_num >> shift;

  // Compute intersection point: a + t * r1
  // Keep in high-precision until final conversion
  const auto ah = a.template type_cast<Wide>();
  const auto r1h = r1.template type_cast<Wide>();
  const auto intersect = (ah * denom_norm + r1h * t_num_norm) / denom_norm;

  result = intersect.template type_cast<CoordType>();
  return true;
}

//------------------------------------------------------------------------------
// 3D segment-plane intersection
//
// Tests if segment [segStart, segEnd] intersects the plane defined by:
//   - Point v0 on the plane
//   - plane_normal (not necessarily unit length, should be normalized to fit
//     in CoordType to prevent overflow)
//
// Returns:
//   - -1: segment does not interest plane
//   - 0: starting point lies on the plane
//   - 1: segment intersects plane
//   - result: intersection point (scaled by denom to avoid division)
//   - denom: normalization factor (divide result by denom to get actual point)
//
// This keeps coordinates in integer space by delaying division.
// Caller is responsible for normalizing if needed.
//
// Note: plane_normal should be bit-normalized if it came from a cross product
// to prevent overflow in the dot product calculations.
//------------------------------------------------------------------------------
template<typename CoordType>
int segmentPlaneIntersection3D(const Point3<CoordType>& segStart,
                               const Point3<CoordType>& segEnd,
                               const Point3<CoordType>& v0,
                               const Point3<CoordType>& plane_normal,
                               Point3<WideInt<3>>& result,
                               WideInt<3>& denom) {
  using Wide = WideInt<3>;

  // Compute signed distances from plane
  const Point3<CoordType> PA = v0 - segStart;
  const Point3<CoordType> PB = v0 - segEnd;

  const Wide dist1 = qdot<3, CoordType>(PA, plane_normal);
  const Wide dist2 = qdot<3, CoordType>(PB, plane_normal);

  // Starting point lies on the plane
  if (dist1 == 0) {
    result = segStart.template type_cast<Wide>();
    denom = 1;
    return 0;
  }

  // Segment entirely on one side of plane
  if ((dist1 > 0 && dist2 > 0) || (dist1 < 0 && dist2 < 0)) {
    return -1;
  }

  // Compute parametric coordinate: t = dist1 / (dist1 - dist2)
  Wide t_num = dist1;
  denom = dist1 - dist2;

  // Normalize sign (make denominator positive)
  if (denom < 0) {
    t_num = -t_num;
    denom = -denom;
  }

  // Normalize denom and t_num to prevent overflow in final computation
  // Target: 42 bits to leave room for coordinate multiplication (42 + 42 = 84 < 128)
  auto shift = bitShiftAmount(denom, 42);
  Wide denom_norm = denom >> shift;
  Wide t_num_norm = t_num >> shift;

  // Compute intersection: segStart + t * (segEnd - segStart)
  // Result is scaled by denom_norm (caller must divide)
  const auto startH = segStart.template type_cast<Wide>();
  const auto dirH = (segEnd - segStart).template type_cast<Wide>();
  result = startH * denom_norm + dirH * t_num_norm;
  denom = denom_norm;

  return 1;
}

//------------------------------------------------------------------------------
// 3D segment-face intersection
//
// Tests if segment [segStart, segEnd] intersects a polygonal face.
//
// Algorithm:
//   1. Compute face plane from first three vertices
//   2. Test segment-plane intersection
//   3. Project intersection point and face to 2D
//   4. Test if 2D point is inside 2D polygon
//
// Returns -1 if not intersection exists.
// Returns 0 if ray is coplanar.
// Otherwise 1 and fills result with intersection point.
//------------------------------------------------------------------------------
template<typename CoordType>
int segmentFaceIntersection3D(const Point3<CoordType>& segStart,
                              const Point3<CoordType>& segEnd,
                              const std::vector<int>& faceIndices,
                              const std::vector<Point3<CoordType>>& vertices,
                              Point3<CoordType>& result) {
  using Wide = WideInt<3>;

  const size_t N = faceIndices.size();
  if (N < 3) return -1; // Degenerate face

  // Get three vertices to define the plane
  const Point3<CoordType>& v0 = vertices[faceIndices[0]];
  const Point3<CoordType>& v1 = vertices[faceIndices[1]];
  const Point3<CoordType>& v2 = vertices[faceIndices[2]];

  // Normalize plane normal to fit in CoordType (direction is all that matters)
  const Point3<CoordType> plane_normal = computeNormal(v0, v1, v2);

  // Find segment-plane intersection (scaled by denom)
  Point3<Wide> intersect;
  Wide denom;
  int plane_intersect = segmentPlaneIntersection3D(segStart, segEnd, v0, plane_normal, intersect, denom);
  if (plane_intersect < 0) {
    return -1;
  }
  bool coplanar = (plane_intersect == 0);

  // Determine projection axis (drop largest component of normal)
  Point3<CoordType> abs_normal = plane_normal;
  if (abs_normal.x < 0) abs_normal.x = -abs_normal.x;
  if (abs_normal.y < 0) abs_normal.y = -abs_normal.y;
  if (abs_normal.z < 0) abs_normal.z = -abs_normal.z;
  const int projAxis = abs_normal.maxAxis();

  // Project face vertices to 2D (scaled by denom for consistency)
  std::vector<Point2<Wide>> face2D_full;
  face2D_full.reserve(N);
  Point2<Wide> max_val;
  for (const auto idx : faceIndices) {
    const auto v = vertices[idx].template type_cast<Wide>() * denom;
    face2D_full.push_back(projectTo2D(v, projAxis));
    max_val = face2D_full.back().maxElements(max_val);
  }

  // Project intersection point to 2D
  Point2<Wide> intersect2D_full = projectTo2D(intersect, projAxis);
  max_val = intersect2D_full.maxElements(max_val);

  // Scale down 2D coordinates to prevent overflow in pointInPolygon
  // This is safe because point-in-polygon tests are scale-invariant
  auto max_axis = max_val.maxAxis();
  auto shift = bitShiftAmount(max_val[max_axis], 42);

  // Normalize everything by bit shifting
  Point2<CoordType> intersect2D = intersect2D_full.template bitShift<CoordType>(shift);
  std::vector<Point2<CoordType>> face2D;
  face2D.reserve(N);
  for (const auto& vec : face2D_full) {
    face2D.push_back(vec.template bitShift<CoordType>(shift));
  }

  // Test if 2D point is inside 2D polygon (now with manageable coordinates)
  if (pointInPolygon(intersect2D, face2D) ||
      pointOnPolygon(intersect2D, face2D)) {
    if (coplanar) {
      result = segStart;
      return 0;
    } else {
      result = (intersect / denom).template type_cast<CoordType>();
      return 1;
    }
  }

  return -1;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Half-space clipping primitives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Layer 1: Point classification relative to a plane
//
// Classifies a point as above (+1), on (0), or below (-1) a plane.
// Plane is defined by a point and an outward-pointing normal.
//
// Returns:
//   +1: point is on the positive side (direction of normal)
//    0: point is on the plane
//   -1: point is on the negative side (opposite of normal)
//------------------------------------------------------------------------------
template<typename CoordType>
int classifyPointByPlane(const Point3<CoordType>& point,
                         const Point3<CoordType>& planePoint,
                         const Point3<CoordType>& planeNormal) {
  const auto ap = point - planePoint;
  const auto dist = qdot<3, CoordType>(ap, planeNormal);
  return (dist < 0 ? -1 : dist > 0 ? 1 : 0);
}

//------------------------------------------------------------------------------
// Layer 2: Edge-plane clipping
//
// Clips an edge against a plane, returning the intersection point if it exists.
// This is a convenience wrapper around segmentPlaneIntersection3D that handles
// the normalization by denom.
//
// Returns:
//   -1: edge does not intersect plane (parallel or no crossing)
//    0: edge start point is on the plane
//    1: edge intersects plane, result contains intersection point
//------------------------------------------------------------------------------
template<typename CoordType>
int clipEdgeByPlane(const Point3<CoordType>& edgeStart,
                    const Point3<CoordType>& edgeEnd,
                    const Point3<CoordType>& planePoint,
                    const Point3<CoordType>& planeNormal,
                    Point3<CoordType>& result) {
  using Wide = WideInt<3>;

  Point3<Wide> result_wide;
  Wide denom;

  int status = segmentPlaneIntersection3D(edgeStart, edgeEnd,
                                          planePoint, planeNormal,
                                          result_wide, denom);

  if (status == 0) {
    result = edgeStart;  // Start point on plane
  } else if (status == 1) {
    result = (result_wide / denom).template type_cast<CoordType>();
  }

  return status;
}

//------------------------------------------------------------------------------
// Result structure for face clipping operations
//------------------------------------------------------------------------------
template<typename CoordType>
struct FaceClipResult {
  std::vector<Point3<CoordType>> vertices;  // Clipped vertices
  bool fullyClipped;   // True if entire face was removed
  bool fullyRetained;  // True if face was completely kept (no clipping)

  FaceClipResult() : fullyClipped(false), fullyRetained(false) {}
};

//------------------------------------------------------------------------------
// Layer 3: Face-plane clipping (Sutherland-Hodgman algorithm)
//
// Clips a polygonal face against a plane, keeping vertices on the "keep" side.
//
// Algorithm:
//   - Classify all vertices relative to plane
//   - Walk edges in order
//   - Keep vertices on the "keep" side
//   - Add intersection points where edges cross the plane
//
// Parameters:
//   faceVertices: ordered vertices of the face (counter-clockwise)
//   planePoint: a point on the clipping plane
//   planeNormal: normal vector of plane (should point toward "keep" side)
//   keepPositive: if true, keep positive side; if false, keep negative side
//
// Returns: FaceClipResult with clipped vertices and status flags
//------------------------------------------------------------------------------
template<typename CoordType>
FaceClipResult<CoordType> clipFaceByPlane(
    const std::vector<Point3<CoordType>>& faceVertices,
    const Point3<CoordType>& planePoint,
    const Point3<CoordType>& planeNormal,
    bool keepPositive = true)
{
  FaceClipResult<CoordType> result;
  const size_t N = faceVertices.size();

  if (N < 3) {
    result.fullyClipped = true;
    return result;  // Degenerate face
  }

  // Classify all vertices (reuses Layer 1)
  std::vector<int> sides(N);
  int numPositive = 0, numNegative = 0, numOnPlane = 0;

  for (size_t i = 0; i < N; ++i) {
    sides[i] = classifyPointByPlane(faceVertices[i], planePoint, planeNormal);
    if (sides[i] > 0) ++numPositive;
    else if (sides[i] < 0) ++numNegative;
    else ++numOnPlane;
  }

  // Early exit: entire face on one side
  if (keepPositive) {
    if (numNegative == 0) {
      result.vertices = faceVertices;
      result.fullyRetained = true;
      return result;
    }
    if (numPositive == 0 && numOnPlane == 0) {
      result.fullyClipped = true;
      return result;
    }
  } else {
    if (numPositive == 0) {
      result.vertices = faceVertices;
      result.fullyRetained = true;
      return result;
    }
    if (numNegative == 0 && numOnPlane == 0) {
      result.fullyClipped = true;
      return result;
    }
  }

  // Sutherland-Hodgman: walk edges and build clipped polygon
  result.vertices.reserve(N + 2);  // Clipping can add at most 2 vertices

  for (size_t i = 0; i < N; ++i) {
    const size_t j = (i + 1) % N;
    const auto& vi = faceVertices[i];
    const auto& vj = faceVertices[j];

    // Adjust signs based on which side we're keeping
    const int si = keepPositive ? sides[i] : -sides[i];
    const int sj = keepPositive ? sides[j] : -sides[j];

    // Keep vi if it's on the "keep" side or on the plane
    if (si >= 0) {
      result.vertices.push_back(vi);
    }

    // If edge crosses plane (signs differ and neither is zero), add intersection
    if ((si > 0 && sj < 0) || (si < 0 && sj > 0)) {
      Point3<CoordType> intersection;
      if (clipEdgeByPlane(vi, vj, planePoint, planeNormal, intersection) == 1) {
        result.vertices.push_back(intersection);
      }
    }
  }

  // Check if clipping produced a degenerate result
  if (result.vertices.size() < 3) {
    result.vertices.clear();
    result.fullyClipped = true;
  }

  return result;
}

//------------------------------------------------------------------------------
// Result structure for polyhedron clipping operations
//------------------------------------------------------------------------------
template<typename CoordType>
struct PolyhedronClipResult {
  std::vector<Point3<CoordType>> vertices;  // All unique vertices
  std::vector<std::vector<int>> faces;      // Faces as vertex indices
  bool fullyClipped;   // True if entire polyhedron was removed
  bool fullyRetained;  // True if no clipping occurred

  PolyhedronClipResult() : fullyClipped(false), fullyRetained(false) {}
};

//------------------------------------------------------------------------------
// Layer 4: Polyhedron-plane clipping
//
// Clips a polyhedron against a plane, keeping the portion on the positive side
// of the plane normal.
//
// Algorithm:
//   1. Classify all vertices relative to the plane
//   2. Clip each face against the plane
//   3. Collect new vertices created by edge-plane intersections
//   4. Build a "cap" face from coplanar intersection points
//   5. Rebuild topology with updated vertex indices
//
// Parameters:
//   vertices: input vertex list
//   faces: input face list (each face is a list of vertex indices)
//   planePoint: a point on the clipping plane
//   planeNormal: normal of plane (points toward the side to keep)
//
// Returns: PolyhedronClipResult with clipped geometry
//
// Notes:
//   - Input polyhedron is assumed to be closed and manifold
//   - Face vertex ordering should follow right-hand rule (outward normals)
//   - The cap face will have vertices ordered to maintain manifold property
//------------------------------------------------------------------------------
template<typename CoordType>
PolyhedronClipResult<CoordType> clipPolyhedronByPlane(
    const std::vector<Point3<CoordType>>& vertices,
    const std::vector<std::vector<int>>& faces,
    const Point3<CoordType>& planePoint,
    const Point3<CoordType>& planeNormal)
{
  PolyhedronClipResult<CoordType> result;
  const size_t NV = vertices.size();
  const size_t NF = faces.size();

  if (NV < 4 || NF < 4) {
    result.fullyClipped = true;
    return result;  // Degenerate polyhedron
  }

  // Classify all input vertices (reuses Layer 1)
  std::vector<int> vertexSides(NV);
  int numPositive = 0, numNegative = 0, numOnPlane = 0;

  for (size_t i = 0; i < NV; ++i) {
    vertexSides[i] = classifyPointByPlane(vertices[i], planePoint, planeNormal);
    if (vertexSides[i] > 0) ++numPositive;
    else if (vertexSides[i] < 0) ++numNegative;
    else ++numOnPlane;
  }

  // Early exit: entire polyhedron on one side
  if (numNegative == 0) {
    result.vertices = vertices;
    result.faces = faces;
    result.fullyRetained = true;
    return result;
  }
  if (numPositive == 0 && numOnPlane == 0) {
    result.fullyClipped = true;
    return result;
  }

  // Build vertex map: input vertex or new intersection point -> output index
  std::map<Point3<CoordType>, int> vertexMap;
  auto getOrCreateVertex = [&](const Point3<CoordType>& v) -> int {
    auto it = vertexMap.find(v);
    if (it != vertexMap.end()) {
      return it->second;
    }
    int newIndex = result.vertices.size();
    result.vertices.push_back(v);
    vertexMap[v] = newIndex;
    return newIndex;
  };

  // Track vertices that lie on the clipping plane (for cap face)
  std::vector<int> capVertexIndices;

  // Clip each face (reuses Layer 3)
  for (size_t iFace = 0; iFace < NF; ++iFace) {
    const auto& face = faces[iFace];

    // Extract face vertices
    std::vector<Point3<CoordType>> faceVerts;
    faceVerts.reserve(face.size());
    for (int idx : face) {
      faceVerts.push_back(vertices[idx]);
    }

    // Clip face by plane
    auto clipResult = clipFaceByPlane(faceVerts, planePoint, planeNormal, true);

    if (clipResult.fullyClipped) {
      continue;  // Face was entirely clipped away
    }

    // Build new face with updated indices
    std::vector<int> newFace;
    newFace.reserve(clipResult.vertices.size());

    for (const auto& v : clipResult.vertices) {
      int vIdx = getOrCreateVertex(v);
      newFace.push_back(vIdx);

      // Check if this vertex is on the clipping plane (for cap)
      if (classifyPointByPlane(v, planePoint, planeNormal) == 0) {
        // Only add if not already in cap list
        if (std::find(capVertexIndices.begin(), capVertexIndices.end(), vIdx)
            == capVertexIndices.end()) {
          capVertexIndices.push_back(vIdx);
        }
      }
    }

    if (newFace.size() >= 3) {
      result.faces.push_back(newFace);
    }
  }

  // Build cap face if we have intersection points on the plane
  if (capVertexIndices.size() >= 3) {
    // Need to order cap vertices properly (project to 2D and sort)
    // Determine projection axis (drop largest component of normal)
    Point3<CoordType> abs_normal = planeNormal;
    if (abs_normal.x < 0) abs_normal.x = -abs_normal.x;
    if (abs_normal.y < 0) abs_normal.y = -abs_normal.y;
    if (abs_normal.z < 0) abs_normal.z = -abs_normal.z;
    const int projAxis = abs_normal.maxAxis();

    // Project cap vertices to 2D
    std::vector<Point2<CoordType>> cap2D;
    cap2D.reserve(capVertexIndices.size());
    for (int vIdx : capVertexIndices) {
      cap2D.push_back(projectTo2D(result.vertices[vIdx], projAxis));
    }

    // Order cap vertices by angle around centroid
    Point2<CoordType> centroid(0, 0);
    for (const auto& p : cap2D) {
      centroid.x += p.x;
      centroid.y += p.y;
    }
    centroid.x /= cap2D.size();
    centroid.y /= cap2D.size();

    // Sort by angle (using atan2 equivalent with cross products)
    std::vector<std::pair<int, int>> indexedAngles;  // (index, sort_key)
    for (size_t i = 0; i < cap2D.size(); ++i) {
      auto v = cap2D[i] - centroid;
      // Use a simple angle proxy: atan2-like ordering
      // For quantized coords, we compare using cross products
      indexedAngles.push_back({capVertexIndices[i], i});
    }

    // Simple bubble sort by angle for small cap faces (TODO: optimize)
    auto angleCompare = [&](const std::pair<int,int>& a, const std::pair<int,int>& b) {
      auto va = cap2D[a.second] - centroid;
      auto vb = cap2D[b.second] - centroid;
      // Cross product to determine angle ordering
      auto cross = qcross<CoordType>(va, vb);
      return cross > 0;
    };
    std::sort(indexedAngles.begin(), indexedAngles.end(), angleCompare);

    // Build ordered cap face (reverse if normal points down)
    std::vector<int> capFace;
    for (const auto& pair : indexedAngles) {
      capFace.push_back(pair.first);
    }

    // Reverse cap face orientation if needed to maintain manifold
    // (cap normal should oppose plane normal since we're capping the cut)
    if (planeNormal[projAxis] > 0) {
      std::reverse(capFace.begin(), capFace.end());
    }

    result.faces.push_back(capFace);
  }

  // Check for degenerate result
  if (result.vertices.size() < 4 || result.faces.size() < 4) {
    result.vertices.clear();
    result.faces.clear();
    result.fullyClipped = true;
  }

  return result;
}

} // namespace polytope

#endif
