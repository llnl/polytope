#ifndef __Polytope_GeomUtils__
#define __Polytope_GeomUtils__

#include "polytope.hh"
#include "Point.hh"
#include "HashKey.hh"
#include "EdgeUtils.hh"
#include "Quantizer.hh"

#ifdef POLYTOPE_ENABLE_TRIANGLE
// From predicates.cc
extern double orient2d(double* a, double* b, double* c);
#endif

namespace polytope {

//------------------------------------------------------------------------------
// Type alias for wide integer type used in overflow-safe arithmetic
// Based on CoordType size to ensure safety when 2D methods are called from 3D
// - 32-bit or smaller coords → int64_t (safe for products)
// - 64-bit coords → __int128 (safe for products)
// - Floating-point types → same type (no widening needed)
//------------------------------------------------------------------------------
template<typename CoordType>
struct WideIntHelper {
  using type = typename std::conditional<
    std::is_floating_point<CoordType>::value,
    CoordType,  // For float/double, return same type
    typename std::conditional<
      (sizeof(CoordType) <= 4),
      int64_t,
      __int128
    >::type
  >::type;
};

template<typename CoordType>
using WideInt = typename WideIntHelper<CoordType>::type;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// General helper routines
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Normalize a scalar by bit-shifting to fit within target bit width
//
// Right-shifts the value so it fits within target_bits.
// Preserves sign and is useful for scale-invariant parametric calculations.
//
// Returns: shift_amount
//------------------------------------------------------------------------------
template<typename IntType>
int bitShiftAmount(IntType value,
                   const int target_bits = 21) {
  if (value == 0) return 0;

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
    return 0;
  }

  return shift_amount;
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
template<typename CoordType, typename Wide>
Point2<CoordType> normalizeByBitShift(const Point2<Wide>& vec,
                                      const int target_bits = 42) {
  // Find largest absolute value component
  Wide abs_x = vec.x < 0 ? -vec.x : vec.x;
  Wide abs_y = vec.y < 0 ? -vec.y : vec.y;
  Wide max_val = std::max({abs_x, abs_y});

  if (max_val == 0) return Point2<CoordType>(0, 0);

  // Determine bit shift amount
  auto shift = bitShiftAmount(max_val, target_bits);

  // Shift all components uniformly and cast to target type
  return vec.template bitShift<CoordType>(shift);
}

template<typename CoordType, typename Wide>
Point3<CoordType> normalizeByBitShift(const Point3<Wide>& vec,
                                      const int target_bits = 42) {

  // Find largest absolute value component
  Wide abs_x = vec.x < 0 ? -vec.x : vec.x;
  Wide abs_y = vec.y < 0 ? -vec.y : vec.y;
  Wide abs_z = vec.z < 0 ? -vec.z : vec.z;
  Wide max_val = std::max({abs_x, abs_y, abs_z});

  if (max_val == 0) return Point3<CoordType>(0, 0, 0);

  // Determine bit shift amount
  auto shift = bitShiftAmount(max_val, target_bits);

  // Shift all components uniformly and cast to target type
  return vec.template bitShift<CoordType>(shift);
}

//------------------------------------------------------------------------------
// Quantized dot product with overflow protection
// Casts to WideInt before multiplication to prevent overflow
//------------------------------------------------------------------------------
template<int Dim, typename CoordType>
WideInt<CoordType> qdot(const Point<Dim, CoordType>& a,
                        const Point<Dim, CoordType>& b) {
  using Wide = WideInt<CoordType>;
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
WideInt<CoordType> qcross(const Point2<CoordType>& a,
                          const Point2<CoordType>& b) {
  using Wide = WideInt<CoordType>;
  const auto ah = a.template type_cast<Wide>();
  const auto bh = b.template type_cast<Wide>();
  return (ah.x * bh.y) - (ah.y * bh.x);
}

//------------------------------------------------------------------------------
// Quantized 3D cross product (returns vector)
// Computes: a × b
//------------------------------------------------------------------------------
template<typename CoordType>
Point3<WideInt<CoordType>> qcross(const Point3<CoordType>& a,
                                  const Point3<CoordType>& b) {
  using Wide = WideInt<CoordType>;
  const auto ah = a.template type_cast<Wide>();
  const auto bh = b.template type_cast<Wide>();

  const Wide cx = (ah.y * bh.z) - (ah.z * bh.y);
  const Wide cy = (ah.z * bh.x) - (ah.x * bh.z);
  const Wide cz = (ah.x * bh.y) - (ah.y * bh.x);

  return Point3<Wide>(cx, cy, cz);
}

template<typename CoordType>
Point3<CoordType> norm_qcross(const Point3<CoordType>& v0,
                              const Point3<CoordType>& v1) {
  const auto cross_full = qcross(v0, v1);
  return normalizeByBitShift<CoordType>(cross_full);
}

//------------------------------------------------------------------------------
// Compare the magnitude of 2 points. Return true if p1 > p2
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
bool magComparison(const Point<Dimension, CoordType>& p1,
                   const Point<Dimension, CoordType>& p2) {
  using Wide = WideInt<CoordType>;
  auto p1w = p1.template type_cast<Wide>();
  auto p2w = p2.template type_cast<Wide>();
  Wide mag1 = 0, mag2 = 0;
  for (int d = 0; d < Dimension; ++d) {
    mag1 += p1w[d]*p1w[d];
    mag2 += p2w[d]*p2w[d];
  }
  return (mag1 > mag2) ? true : false;
}

//------------------------------------------------------------------------------
// Get the direction between two double points into an integer type.
//------------------------------------------------------------------------------
template<typename CoordType>
Point2<CoordType> pointDirection(const Point2<double>& p1,
                                 const Point2<double>& p2) {
  auto diff = p2 - p1;
  double len = std::hypot(diff[0], diff[1]);
  if (len == 0.) {
    return Point<2, CoordType>::Zero();
  }
  auto norm = diff/len;
  const double SCALE = std::pow(2.0, HashKey<2>::num1DBits() - 2);
  return (norm*SCALE).template type_cast<CoordType>();
}

//------------------------------------------------------------------------------
// Determine if a ray (an origin and direction) is completely external
// to a bounding box
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
bool isRayExternal(const Point<Dimension, double>& origin,
                   const Point<Dimension, CoordType>& dir) {
  auto& Q = Quantizer<Dimension>::instance();
  Point<Dimension, int> outdirs = Q.externalSides(origin);
  for (int d = 0; d < Dimension; ++d) {
    if (outdirs[d] < 0 && dir[d] < 0) {
      return true;
    } else if (outdirs[d] > 0 && dir[d] > 0) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Given two generator points, return the midpoint
//------------------------------------------------------------------------------
template<typename CoordType>
Point2<CoordType> midPoint(const Point2<CoordType>& gen0,
                           const Point2<CoordType>& gen1) {
  using Wide = WideInt<CoordType>;
  Point2<Wide> sum = gen0.template type_cast<Wide>() + gen1.template type_cast<Wide>();
  return (sum/2).template type_cast<CoordType>();
}

//------------------------------------------------------------------------------
// Given two generator points, return a vector normal to the segment between them
//------------------------------------------------------------------------------
template<typename CoordType>
Point2<CoordType> outwardRay(const Point2<CoordType>& gen0,
                             const Point2<CoordType>& gen1) {
  auto diff = gen1 - gen0;
  return Point2<CoordType>(-diff.y, diff.x);
}

//------------------------------------------------------------------------------
// Given 3 points on a triangle and the circumcenter, determine the outward
// ray direction for the a->b edge.
//------------------------------------------------------------------------------
template<typename CoordType>
Point2<CoordType> outwardRay(const Point2<CoordType>& a,
                             const Point2<CoordType>& b,
                             const Point2<CoordType>& c) {
  //const Point2<double>& circent) {
  auto diff = a - b;
  auto ad = a.template type_cast<double>();
  auto bd = b.template type_cast<double>();
  auto cd = c.template type_cast<double>();
  auto delta = (bd.x - ad.x)*(cd.y - ad.y) - (bd.y - ad.y)*(cd.x - ad.x);
  if (delta > 0.) {
    return Point2<CoordType>(-diff.y, diff.x);
  } else {
    return Point2<CoordType>(diff.y, -diff.x);
  }
  // auto diff = a - b;
  // double test1 = (static_cast<double>(c.x) - circent.x)*static_cast<double>(a.y - b.y);
  // double test2 = (static_cast<double>(c.y) - circent.y)*static_cast<double>(b.x - a.x);
  // if (test1 > -test2) {
  //   return Point2<CoordType>(-diff.y, diff.x);
  // } else {
  //   return Point2<CoordType>(diff.y, -diff.x);
  // }
}

//------------------------------------------------------------------------------
// From 3 points, compute the plane normal and normalize it by bit shifting
//------------------------------------------------------------------------------
template<typename CoordType>
Point3<CoordType> computeNormal(const Point3<CoordType>& v0,
                                const Point3<CoordType>& v1,
                                const Point3<CoordType>& v2) {
  const auto edge1 = v1 - v0;
  const auto edge2 = v2 - v0;
  const auto plane_normal_full = qcross<CoordType>(edge1, edge2);
  return normalizeByBitShift<CoordType>(plane_normal_full);
}

//------------------------------------------------------------------------------
// Determine if a point is collinear with a line
//   1. Check collinearity: (point - vi) × (vj - vi) == 0
//   2. Check if point is between vi and vj using bounding box test
//------------------------------------------------------------------------------
template<typename CoordType>
bool collinear(const Point2<CoordType>& segStart,
               const Point2<CoordType>& segEnd,
               const Point2<CoordType>& point) {
  auto vp = point - segStart;
  auto vj = segEnd - segStart;
  if (qcross(vp, vj) != 0) {
    return false;
  }
  // Point is collinear with line - check if it's between endpoints
  // Use bounding box test for each coordinate
  auto min = segStart.minElements(segEnd);
  auto max = segStart.maxElements(segEnd);
  if (point.allGreaterEqual(min) && point.allLessEqual(max)) {
    return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Project a 3D point to 2D by dropping the specified axis
//   projAxis = 0: drop x, keep (y, z)
//   projAxis = 1: drop y, keep (x, z)
//   projAxis = 2: drop z, keep (x, y)
//------------------------------------------------------------------------------
template<typename CoordType>
Point2<CoordType> projectTo2D(const Point3<CoordType>& point,
                              const int projAxis) {
  switch (projAxis) {
    case 0:  return Point2<CoordType>(point.y, point.z);
    case 1:  return Point2<CoordType>(point.x, point.z);
    default: return Point2<CoordType>(point.x, point.y);
  }
}

//------------------------------------------------------------------------------
// Compare a point to a line or plane to determine if the point is inside the
// interior half-plane, inside the exterior half-plane, or collinear/coplanar.
// Returns -1, 1, or 0, respectively.
// p: Point of interest
// 2D version takes the start and end of the facet.
// 3D version takes a point on the facets plane and the bit-shifted normalized
// plane normal.
//------------------------------------------------------------------------------
template<typename CoordType>
int aboveBelow(const Point2<CoordType>& a,
               const Point2<CoordType>& b,
               const Point2<CoordType>& p) {
  auto ab = b - a;
  auto ap = p - a;
  auto ztest = qcross<CoordType>(ab, ap);
  return -(ztest < 0 ? -1 :
           ztest > 0 ? 1 :
           0);
}

template<typename CoordType>
int aboveBelow(const Point3<CoordType>& a,
               const Point3<CoordType>& n,
               const Point3<CoordType>& p) {
  auto ap = p - a;
  auto ztest = qdot<3>(ap, n);
  return (ztest < 0 ? -1 :
          ztest > 0 ? 1 :
          0);
}

template<int Dimension, typename CoordType>
bool aboveBelow(const Point<Dimension, CoordType>& v0,
                const Point<Dimension, CoordType>& n,
                const std::vector<Point<Dimension, CoordType>>& points) {
  const auto NP = points.size();
  const int result = aboveBelow(v0, n, points[0]);
  if (result == 0) {
    return true;
  }
  for (auto i = 1u; i < NP; ++i) {
    if (result != aboveBelow(v0, n, points[i])) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Determine if a separating axis exists
//------------------------------------------------------------------------------
template<typename CoordType>
bool SAT(const std::vector<Point2<CoordType>>& pointsA,
         const std::vector<Point2<CoordType>>& pointsB,
         const Point2<WideInt<CoordType>>& axis) {
  using AxisType = WideInt<CoordType>;
  using Projection = WideInt<AxisType>;
  POLY_ASSERT(!pointsA.empty());
  POLY_ASSERT(!pointsB.empty());
  if (axis.iszero()) return false;

  auto project = [&axis](const Point2<CoordType>& p) {
    return static_cast<Projection>(p.x)*static_cast<Projection>(axis.x) +
           static_cast<Projection>(p.y)*static_cast<Projection>(axis.y);
  };

  Projection minA = project(pointsA.front()), maxA = minA;
  for (const auto& p : pointsA) {
    const auto ztest = project(p);
    if (ztest < minA) {
      minA = ztest;
    }
    if (ztest > maxA) {
      maxA = ztest;
    }
  }

  Projection minB = project(pointsB.front()), maxB = minB;
  for (const auto& p : pointsB) {
    const auto ztest = project(p);
    if (ztest < minB) {
      minB = ztest;
    }
    if (ztest > maxB) {
      maxB = ztest;
    }
  }
  return maxA < minB || maxB < minA;
}

template<typename CoordType>
bool SAT(const std::vector<Point3<CoordType>>& pointsA,
         const std::vector<Point3<CoordType>>& pointsB,
         const Point3<CoordType>& axis) {
  using Wide = WideInt<CoordType>;
  Wide minA = HashKey<3>::hashMax();
  Wide minB = minA, maxA = -minA, maxB = maxA;
  for (const auto& p : pointsA) {
    auto ztest = qdot<3>(p, axis);
    if (ztest < minA) {
      minA = ztest;
    }
    if (ztest > maxA) {
      maxA = ztest;
    }
  }
  for (const auto& p : pointsB) {
    auto ztest = qdot<3>(p, axis);
    if (ztest < minB) {
      minB = ztest;
    }
    if (ztest > maxB) {
      maxB = ztest;
    }
  }
  return maxA < minB || maxB < minA;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Circumcenter operations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef POLYTOPE_ENABLE_TRIANGLE
// // Returns a point that is the circumcenter
// inline
// Point2<double> circumcenter(const Point2<double>& a,
//                             const Point2<double>& b,
//                             const Point2<double>& c) {
//   // auto ba = b - a;
//   // auto ca = c - a;
//   // auto bc = b - c;
//   // double d = 2*(a.x*bc.y + b.x*ca.y - c.x*ba.y);
//   double a2 = a.x*a.x + a.y*a.y;
//   double b2 = b.x*b.x + b.y*b.y;
//   double c2 = c.x*c.x + c.y*c.y;
//   // Point2<double> out;
//   // out.x = (a2*bc.y + b2*ca.y - c2*ba.y)/d;
//   // out.y = (-a2*bc.x - b2*ca.x + c2*ba.x)/d;
//   // return out;
//   double ap[2] = {a[0], a[1]};
//   double bp[2] = {b[0], b[1]};
//   double cp[2] = {c[0], c[1]};
//   double d = 2*orient2d(ap, bp, cp);
//   double a0[2] = {a2, a[1]}, a1[2] = {a[0], a2};
//   double b0[2] = {b2, b[1]}, b1[2] = {b[0], b2};
//   double c0[2] = {c2, c[1]}, c1[2] = {c[0], c2};
//   return Point2<double>(orient2d(a0,b0,c0)/d, orient2d(a1,b1,c1)/d);
// }

// Returns a point that is the circumcenter
inline
Point2<double> circumcenter(const Point2<double>& a,
                            const Point2<double>& b,
                            const Point2<double>& c) {
  // Differences are formed in long double before subtraction overflow
  const long double ax = static_cast<long double>(a.x);
  const long double ay = static_cast<long double>(a.y);
  const long double bx = static_cast<long double>(b.x);
  const long double by = static_cast<long double>(b.y);
  const long double cx = static_cast<long double>(c.x);
  const long double cy = static_cast<long double>(c.y);

  const long double abx = bx - ax;
  const long double aby = by - ay;
  const long double acx = cx - ax;
  const long double acy = cy - ay;
  const long double cross = abx * acy - aby * acx;
  if (cross == 0.0L) {
    return Point2<double>(0., 0.);
  }
  const long double ab2 = abx * abx + aby * aby;
  const long double ac2 = acx * acx + acy * acy;
  // Circumcenter relative to A:
  //
  // U = A + (ac2 * perp(AB) - ab2 * perp(AC)) / (2 * cross)
  //
  const long double ux =
    ax + (acy * ab2 - aby * ac2) / (2.0L * cross);
  const long double uy =
    ay + (abx * ac2 - acx * ab2) / (2.0L * cross);
  return Point2<double>(static_cast<double>(ux), static_cast<double>(uy));
}
#endif

#ifdef POLYTOPE_ENABLE_TETGEN
inline
Point3<double>
circumcenter(const Point3<double>& p0,
             const Point3<double>& p1,
             const Point3<double>& p2,
             const Point3<double>& p3) {
  const long double x0 = static_cast<long double>(p0.x);
  const long double y0 = static_cast<long double>(p0.y);
  const long double z0 = static_cast<long double>(p0.z);

  const long double ax = static_cast<long double>(p1.x) - x0;
  const long double ay = static_cast<long double>(p1.y) - y0;
  const long double az = static_cast<long double>(p1.z) - z0;

  const long double bx = static_cast<long double>(p2.x) - x0;
  const long double by = static_cast<long double>(p2.y) - y0;
  const long double bz = static_cast<long double>(p2.z) - z0;

  const long double cx = static_cast<long double>(p3.x) - x0;
  const long double cy = static_cast<long double>(p3.y) - y0;
  const long double cz = static_cast<long double>(p3.z) - z0;

  // Solve:
  //
  // 2*a dot u = |a|^2
  // 2*b dot u = |b|^2
  // 2*c dot u = |c|^2
  //
  const long double rhs0 = (ax * ax + ay * ay + az * az) / 2.0L;
  const long double rhs1 = (bx * bx + by * by + bz * bz) / 2.0L;
  const long double rhs2 = (cx * cx + cy * cy + cz * cz) / 2.0L;

  // Matrix:
  //
  // [ ax ay az ]
  // [ bx by bz ] u = rhs
  // [ cx cy cz ]
  //
  const long double det =
    ax * (by * cz - bz * cy)
    - ay * (bx * cz - bz * cx)
    + az * (bx * cy - by * cx);

  // Scale-aware degeneracy test.
  const long double scale =
    std::max({
              std::fabs(ax), std::fabs(ay), std::fabs(az),
              std::fabs(bx), std::fabs(by), std::fabs(bz),
              std::fabs(cx), std::fabs(cy), std::fabs(cz),
              1.0L
      });

  constexpr long double epsilon =
        64.0L * std::numeric_limits<long double>::epsilon();

  if (std::fabs(det) <= epsilon * scale * scale * scale) {
    return Point3<double>(0., 0., 0.); // Coplanar or numerically degenerate
  }

  // Cramer's rule.
  const long double detX =
    rhs0 * (by * cz - bz * cy)
    - ay   * (rhs1 * cz - bz * rhs2)
    + az   * (rhs1 * cy - by * rhs2);

  const long double detY =
    ax   * (rhs1 * cz - bz * rhs2)
    - rhs0 * (bx * cz - bz * cx)
    + az   * (bx * rhs2 - rhs1 * cx);

  const long double detZ =
    ax   * (by * rhs2 - rhs1 * cy)
    - ay   * (bx * rhs2 - rhs1 * cx)
    + rhs0 * (bx * cy - by * cx);

  return Point3<double>(x0 + detX / det,
                        y0 + detY / det,
                        z0 + detZ / det);
}
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Topology utilities
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Orient a 3D facet so its normal points away from a reference point (centroid)
// Version 2: Uses precomputed normal
//------------------------------------------------------------------------------
template<typename CoordType, typename WideType>
void orientFacetOutward(std::vector<unsigned>& facet,
                        const std::vector<Point3<CoordType>>& points,
                        const Point3<CoordType>& precomputedNormal,
                        const Point3<WideType>& centroid,
                        WideType numPoints) {
  if (facet.size() < 3) return;

  // Vector from centroid to facet (use v0 as representative point)
  const auto& v0 = points[facet[0]];
  auto toFace = normalizeByBitShift<CoordType>(numPoints * v0.template type_cast<WideType>() - centroid);
  WideType dot = qdot<3, CoordType>(precomputedNormal, toFace);

  // If dot < 0, normal points inward - reverse vertex order
  if (dot < 0) {
    std::reverse(facet.begin(), facet.end());
  }
}

//------------------------------------------------------------------------------
// Remove coplanar adjacent faces by merging them
// Uses precomputed normals (updates normal array as faces merge)
//------------------------------------------------------------------------------
template<typename CoordType>
void mergeCoplanarFaces(std::vector<std::vector<unsigned>>& faces,
                        std::vector<Point3<CoordType>>& normals,
                        const std::vector<Point3<CoordType>>& points) {
  if (faces.size() < 2) return;

  bool merged = true;
  while (merged) {
    merged = false;
    // Try to find a pair of adjacent coplanar faces to merge
    for (size_t i = 0; i < faces.size()-1 && !merged; ++i) {
      edge::EdgeToFaceMap uniqueEdges;
      edge::addUniqueEdges(faces[i], uniqueEdges);
      for (size_t j = i + 1; j < faces.size() && !merged; ++j) {
        if (normals[i] != normals[j]) continue;
        if (!edge::sharedEdges(faces[j], uniqueEdges)) continue;
        merged = true;
        faces.erase(faces.begin() + j);
        normals.erase(normals.begin() + j);
        // Get new merged face
        faces[i] = edge::traceBoundary(uniqueEdges);
        // Reverse order if normals are different
        if (computeNormal(points[faces[i][0]],
                          points[faces[i][1]],
                          points[faces[i][2]]) == -normals[i]) {
          std::reverse(faces[i].begin(), faces[i].end());
        }
      }
    }
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Precomputation helpers for normals and centroids
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Compute normalized normal for a single face
//------------------------------------------------------------------------------
template<typename CoordType>
Point3<CoordType> computeFaceNormal(const std::vector<unsigned>& face,
                                    const std::vector<Point3<CoordType>>& vertices) {
  if (face.size() < 3) return Point3<CoordType>(0, 0, 0);

  const auto& v0 = vertices[face[0]];
  const auto& v1 = vertices[face[1]];
  const auto& v2 = vertices[face[2]];

  return computeNormal(v0, v1, v2);
}

//------------------------------------------------------------------------------
// Compute normalized normals for all faces
//------------------------------------------------------------------------------
template<typename CoordType>
std::vector<Point3<CoordType>> computeFaceNormals(
    const std::vector<Point3<CoordType>>& vertices,
    const std::vector<std::vector<unsigned>>& faces) {

  std::vector<Point3<CoordType>> normals;
  normals.reserve(faces.size());

  for (const auto& face : faces) {
    normals.push_back(computeFaceNormal(face, vertices));
  }

  return normals;
}

//------------------------------------------------------------------------------
// Compute centroid for a single face (as unnormalized sum)
// Returns: (sum of vertices, count of vertices)
//------------------------------------------------------------------------------
template<typename CoordType>
std::pair<Point3<WideInt<CoordType>>, WideInt<CoordType>>
computeFaceCentroid(const std::vector<unsigned>& face,
                    const std::vector<Point3<CoordType>>& vertices) {
  using Wide = WideInt<CoordType>;
  Point3<Wide> sum(0, 0, 0);

  for (int idx : face) {
    sum = sum + vertices[idx].template type_cast<Wide>();
  }

  return {sum, static_cast<Wide>(face.size())};
}

//------------------------------------------------------------------------------
// Compute centroids for all faces (as unnormalized sums)
// Returns vector of (sum, count) pairs
//------------------------------------------------------------------------------
template<typename CoordType>
std::vector<std::pair<Point3<WideInt<CoordType>>, WideInt<CoordType>>>
computeFaceCentroids(const std::vector<Point3<CoordType>>& vertices,
                     const std::vector<std::vector<unsigned>>& faces) {
  using Wide = WideInt<CoordType>;
  std::vector<std::pair<Point3<Wide>, Wide>> centroids;
  centroids.reserve(faces.size());

  for (const auto& face : faces) {
    centroids.push_back(computeFaceCentroid(face, vertices));
  }

  return centroids;
}

//------------------------------------------------------------------------------
// Compute volume-weighted centroid for a polyhedron (as unnormalized sum)
// Uses tetrahedral decomposition from origin
// Returns: (weighted sum, total weight)
//------------------------------------------------------------------------------
template<typename CoordType>
std::pair<Point3<WideInt<CoordType>>, WideInt<CoordType>>
computePolyhedronCentroid(const std::vector<Point3<CoordType>>& vertices,
                          const std::vector<std::vector<unsigned>>& faces) {
  using Wide = WideInt<CoordType>;
  Point3<Wide> weightedSum(0, 0, 0);
  Wide totalWeight = 0;

  // Decompose polyhedron into tetrahedra from origin
  for (const auto& face : faces) {
    if (face.size() < 3) continue;

    // Triangulate face (simple fan triangulation)
    const auto& v0 = vertices[face[0]];
    for (size_t i = 1; i + 1 < face.size(); ++i) {
      const auto& v1 = vertices[face[i]];
      const auto& v2 = vertices[face[i + 1]];

      // Signed volume of tetrahedron (origin, v0, v1, v2) = dot(v0, cross(v1, v2)) / 6
      // We'll skip the /6 and accumulate the unnormalized weight
      auto cross_full = qcross(v1.template type_cast<Wide>(),
                               v2.template type_cast<Wide>());
      Wide signedVol6 = qdot<3>(v0.template type_cast<Wide>(), cross_full);

      // Centroid of tetrahedron is at (v0 + v1 + v2) / 4 (plus origin / 4)
      // Weighted contribution: signedVol6 * (v0 + v1 + v2) / 4
      // We keep unnormalized: signedVol6 * (v0 + v1 + v2)
      auto tetCenter = v0.template type_cast<Wide>() +
                       v1.template type_cast<Wide>() +
                       v2.template type_cast<Wide>();

      // Normalize both operands before multiplication to prevent overflow
      // signedVol6 can be ~coord^3, tetCenter components can be ~3*coord
      // Product would be ~coord^4 which overflows __int128 for large coords
      // Shift both to fit in ~60 bits so product fits in ~120 bits (< 127 for __int128)
      Wide vol_abs = signedVol6 < 0 ? -signedVol6 : signedVol6;
      Wide max_center = std::max({tetCenter.x < 0 ? -tetCenter.x : tetCenter.x,
                                  tetCenter.y < 0 ? -tetCenter.y : tetCenter.y,
                                  tetCenter.z < 0 ? -tetCenter.z : tetCenter.z});

      int vol_shift = bitShiftAmount(vol_abs, 60);
      int center_shift = bitShiftAmount(max_center, 60);

      Wide vol_norm = signedVol6 >> vol_shift;
      auto center_norm = tetCenter.template bitShift<Wide>(center_shift);

      // Accumulate with normalized values
      weightedSum = weightedSum + center_norm * vol_norm;
      totalWeight += vol_norm;
    }
  }

  // Return unnormalized weighted sum and total weight
  // Actual centroid would be: weightedSum / (4 * totalWeight)
  // But we keep it unnormalized to avoid division
  return {weightedSum, totalWeight * 4};
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Double and pointer operations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int Dimension, typename RealType>
RealType dot(const RealType* a,
             const RealType* b) {
  RealType sum = 0.;
  for (int i = 0; i < Dimension; ++i) {
    sum += a[i]*b[i];
  }
  return sum;
}

template<int Dimension>
double distance(const Point<Dimension, double>& a,
                const Point<Dimension, double>& b) {
  return magnitude(a - b);
}

//------------------------------------------------------------------------------
// Determine if the given points are collinear to some accuracy.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
collinear(const RealType* a, const RealType* b, const RealType* c, const RealType tol) {
  double ab[Dimension], ac[Dimension], abmag = 0.0, acmag = 0.0;
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] = b[j] - a[j];
    ac[j] = c[j] - a[j];
    abmag += ab[j]*ab[j];
    acmag += ac[j]*ac[j];
  }
  if (abmag < tol or acmag < tol) return true;
  abmag = std::sqrt(abmag);
  acmag = std::sqrt(acmag);
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] /= abmag;
    ac[j] /= acmag;
  }
  return std::abs(std::abs(dot<Dimension, double>(ab, ac)) - 1.0) < tol;
}

template<int Dimension>
double magnitude(const Point<Dimension, double>& a) {
  long double dis = 0.;
  for (int d = 0; d < Dimension; ++d) {
    auto dd = static_cast<long double>(a[d]);
    dis += dd*dd;
  }
  dis = std::sqrt(dis);
  return static_cast<double>(dis);
}

template<int Dimension>
long double magnitude(const Point<Dimension, long double>& a) {
  long double dis = 0.;
  for (int d = 0; d < Dimension; ++d) {
    auto dd = a[d];
    dis += dd*dd;
  }
  dis = std::sqrt(dis);
  return dis;
}

template<int Dimension>
Point<Dimension, double> triangleCentroid(const Point<Dimension, double>& a,
                                          const Point<Dimension, double>& b,
                                          const Point<Dimension, double>& c) {
  Point<Dimension, long double> out(0.);
  for (int d = 0; d < Dimension; ++d) {
    long double sum = static_cast<long double>(a[d]) +
      static_cast<long double>(b[d]) +
      static_cast<long double>(c[d]);
    out[d] = sum/3.;
  }
  return out.template type_cast<double>();
}

template<int Dimension, typename RealType>
void UnitVector(RealType* a) {
  if constexpr (Dimension == 2) {
    const RealType mag = std::max(1.E-100, std::sqrt(a[0]*a[0] + a[1]*a[1]));
    a[0] /= mag;
    a[1] /= mag;
  } else if constexpr (Dimension == 3) {
    const RealType mag = std::max(1.0e-100, sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
    a[0] /= mag;
    a[1] /= mag;
    a[2] /= mag;
  }
}

template<int Dimension, typename RealType>
Point<Dimension, RealType> cross(const Point<Dimension, RealType>& a,
                                 const Point<Dimension, RealType>& b) {
  Point<Dimension, RealType> out;
  if constexpr (Dimension == 2) {
    out[0] = a[0]*b[1] - a[1]*b[0];
  } else if constexpr (Dimension == 3) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
  }
  return out;
}

template<int Dimension, typename RealType>
void cross(const RealType* a,
           const RealType* b,
           RealType* c) {
  if constexpr (Dimension == 2) {
    c[2] = a[0]*b[1] - a[1]*b[0];
  } else if constexpr (Dimension == 3) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
}

template<int Dimension>
Point<Dimension, double> normal(const Point<Dimension, double>& a,
                                const Point<Dimension, double>& b,
                                const Point<Dimension, double>& c) {
  auto ad = a.template type_cast<long double>();
  auto bd = b.template type_cast<long double>();
  auto cd = c.template type_cast<long double>();
  auto cprod = cross(bd - ad, cd - ad);
  long double mag = magnitude(cprod);
  return (cprod/mag).template type_cast<double>();
}

}
#endif
