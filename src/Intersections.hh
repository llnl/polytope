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

#include "Point.hh"
#include "HashKey.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Type alias for wide integer type used in overflow-safe arithmetic
// Borrows from HashKey<Dim>::CoordHash since it defines the appropriate
// wide integer for each dimension (int64_t for 2D, __int128 for 3D)
//------------------------------------------------------------------------------
template<int Dim>
using WideInt = typename HashKey<Dim>::CoordHash;

//------------------------------------------------------------------------------
// Normalize a scalar by bit-shifting to fit within target bit width
//
// Right-shifts the value so it fits within target_bits.
// Preserves sign and is useful for scale-invariant parametric calculations.
//
// Returns: (normalized_value, shift_amount)
//------------------------------------------------------------------------------
template<typename IntType, typename WideInt>
std::pair<IntType, int> normalizeScalarByBitShift(WideInt value,
                                                  const int target_bits = 21) {
  if (value == 0) return {value, 0};

  // Get absolute value for bit counting
  IntType abs_val = value < 0 ? -value : value;

  // Count bits used
  int bits_used = 0;
  IntType temp = abs_val;
  while (temp > 0) {
    ++bits_used;
    temp >>= 1;
  }

  // Determine shift needed
  const int shift_amount = std::max(0, bits_used - target_bits);

  if (shift_amount == 0) {
    return {value, 0};
  }

  return {value >> shift_amount, shift_amount};
}

//------------------------------------------------------------------------------
// Normalize a vector by bit-shifting to fit within target bit width
//
// Finds the largest component (by absolute value) and right-shifts all
// components uniformly so the max fits within target_bits.
//
// This preserves direction while reducing magnitude, making it safe to use
// in scale-invariant calculations (e.g., parametric line-plane intersection).
//
// target_bits should be less than the bit width of CoordType to leave room
// for subsequent operations.
//------------------------------------------------------------------------------
template<typename CoordType>
Point3<CoordType> normalizeByBitShift(const Point3<WideInt<3>>& vec,
                                      const int target_bits = 42) {
  using Wide = WideInt<3>;
  // Find largest absolute value component
  Wide abs_x = vec.x < 0 ? -vec.x : vec.x;
  Wide abs_y = vec.y < 0 ? -vec.y : vec.y;
  Wide abs_z = vec.z < 0 ? -vec.z : vec.z;
  Wide max_val = std::max({abs_x, abs_y, abs_z});

  if (max_val == 0) return Point3<CoordType>(0, 0, 0); // Zero vector

  // Count leading zeros to determine bit width (portable version)
  // This could use __builtin_clzll for int64_t or compiler intrinsics for __int128
  int bits_used = 0;
  Wide temp = max_val;
  while (temp > 0) {
    ++bits_used;
    temp >>= 1;
  }

  // Determine shift needed to fit within target bit width
  const int shift_amount = std::max(0, bits_used - target_bits);

  // Shift all components uniformly and cast to target type
  return Point3<CoordType>(
    static_cast<CoordType>(vec.x >> shift_amount),
    static_cast<CoordType>(vec.y >> shift_amount),
    static_cast<CoordType>(vec.z >> shift_amount)
  );
}

//------------------------------------------------------------------------------
// Quantized dot product with overflow protection
// Casts to WideInt before multiplication to prevent overflow
//------------------------------------------------------------------------------
template<int Dim, typename CoordType>
WideInt<Dim> qdot(const Point<Dim, CoordType>& a,
                  const Point<Dim, CoordType>& b) {
  using Wide = WideInt<Dim>;
  const auto ah = a.template type_cast<Wide>();
  const auto bh = b.template type_cast<Wide>();

  Wide sum = 0;
  for (int i = 0; i < Dim; ++i) {
    sum += ah[i] * bh[i];
  }
  return sum;
}

//------------------------------------------------------------------------------
// Quantized 2D cross product (returns scalar z-component)
// Computes: a.x * b.y - a.y * b.x
//------------------------------------------------------------------------------
template<typename CoordType>
WideInt<2> qcross(const Point2<CoordType>& a,
                  const Point2<CoordType>& b) {
  using Wide = WideInt<2>;
  const auto ah = a.template type_cast<Wide>();
  const auto bh = b.template type_cast<Wide>();
  return (ah.x * bh.y) - (ah.y * bh.x);
}

//------------------------------------------------------------------------------
// Quantized 3D cross product (returns vector)
// Computes: a × b
//------------------------------------------------------------------------------
template<typename CoordType>
Point3<WideInt<3>> qcross(const Point3<CoordType>& a,
                          const Point3<CoordType>& b) {
  using Wide = WideInt<3>;
  const auto ah = a.template type_cast<Wide>();
  const auto bh = b.template type_cast<Wide>();

  const Wide cx = (ah.y * bh.z) - (ah.z * bh.y);
  const Wide cy = (ah.z * bh.x) - (ah.x * bh.z);
  const Wide cz = (ah.x * bh.y) - (ah.y * bh.x);

  return Point3<Wide>(cx, cy, cz);
}

//------------------------------------------------------------------------------
// Project a 3D point to 2D by dropping the specified axis
//   projAxis = 0: drop x, keep (y, z)
//   projAxis = 1: drop y, keep (x, z)
//   projAxis = 2: drop z, keep (x, y)
//------------------------------------------------------------------------------
template<typename IntType>
Point2<IntType> projectTo2D(const Point3<IntType>& point,
                            const int projAxis) {
  switch (projAxis) {
    case 0:  return Point2<IntType>(point.y, point.z);
    case 1:  return Point2<IntType>(point.x, point.z);
    default: return Point2<IntType>(point.x, point.y);
  }
}

//------------------------------------------------------------------------------
// Point-in-polygon test using ray casting algorithm
//
// Casts a horizontal ray from the point to the right and counts edge crossings.
// Odd count = inside, even count = outside.
//
// Handles vertex coincidence by shifting test point up slightly when it lies
// exactly on a horizontal edge endpoint.
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

    if (lhs <= rhs) {
      inside = !inside;
    }
  }

  return inside;
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
  auto [denom_norm, shift] = normalizeScalarByBitShift<CoordType, Wide>(denom, 21);
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
// Returns true if intersection exists and fills:
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
bool segmentPlaneIntersection3D(const Point3<CoordType>& segStart,
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

  // Segment is coplanar with face - any point on segment is valid
  if (dist1 == 0 && dist2 == 0) {
    result = segStart.template type_cast<Wide>();
    denom = 1;
    return true;
  }

  // Segment entirely on one side of plane
  if ((dist1 > 0 && dist2 > 0) || (dist1 < 0 && dist2 < 0)) {
    return false;
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
  auto [denom_norm, shift] = normalizeScalarByBitShift<CoordType, Wide>(denom, 42);
  Wide t_num_norm = t_num >> shift;

  // Compute intersection: segStart + t * (segEnd - segStart)
  // Result is scaled by denom_norm (caller must divide)
  const auto startH = segStart.template type_cast<Wide>();
  const auto dirH = (segEnd - segStart).template type_cast<Wide>();
  result = startH * denom_norm + dirH * t_num_norm;
  denom = denom_norm;

  return true;
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
// Returns true if intersection exists and fills result with intersection point.
//------------------------------------------------------------------------------
template<typename CoordType>
bool segmentFaceIntersection3D(const Point3<CoordType>& segStart,
                               const Point3<CoordType>& segEnd,
                               const std::vector<int>& faceIndices,
                               const std::vector<Point3<CoordType>>& vertices,
                               Point3<CoordType>& result) {
  using Wide = WideInt<3>;

  const size_t N = faceIndices.size();
  if (N < 3) return false; // Degenerate face

  // Get three vertices to define the plane
  const Point3<CoordType> v0 = vertices[faceIndices[0]];
  const Point3<CoordType> v1 = vertices[faceIndices[1]];
  const Point3<CoordType> v2 = vertices[faceIndices[2]];

  // Compute face normal: (v1 - v0) × (v2 - v0)
  const Point3<CoordType> edge1 = v1 - v0;
  const Point3<CoordType> edge2 = v2 - v0;
  const Point3<Wide> plane_normal_full = qcross<CoordType>(edge1, edge2);

  // Normalize plane normal to fit in CoordType (direction is all that matters)
  const Point3<CoordType> plane_normal = normalizeByBitShift<CoordType>(plane_normal_full);

  // Find segment-plane intersection (scaled by denom)
  Point3<Wide> intersect;
  Wide denom;
  if (!segmentPlaneIntersection3D(segStart, segEnd, v0, plane_normal, intersect, denom)) {
    return false;
  }

  // Determine projection axis (drop largest component of normal)
  Point3<CoordType> abs_normal = plane_normal;
  if (abs_normal.x < 0) abs_normal.x = -abs_normal.x;
  if (abs_normal.y < 0) abs_normal.y = -abs_normal.y;
  if (abs_normal.z < 0) abs_normal.z = -abs_normal.z;
  const int projAxis = abs_normal.maxAxis();

  // Project face vertices to 2D (scaled by denom for consistency)
  std::vector<Point2<Wide>> face2D;
  face2D.reserve(N);
  for (const auto idx : faceIndices) {
    const auto v = vertices[idx].template type_cast<Wide>() * denom;
    face2D.push_back(projectTo2D(v, projAxis));
  }

  // Project intersection point to 2D
  const Point2<Wide> intersect2D = projectTo2D(intersect, projAxis);

  // Test if 2D point is inside 2D polygon
  if (pointInPolygon(intersect2D, face2D)) {
    result = (intersect / denom).template type_cast<CoordType>();
    return true;
  }

  return false;
}



} // namespace polytope

#endif
