#ifndef POLYTOPE_GEOMUTILS_HH
#define POLYTOPE_GEOMUTILS_HH

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

template<int Dim>
using NarrowInt = typename HashKey<Dim>::IntType;

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

template<typename CoordType>
Point3<CoordType> norm_qcross(const Point3<CoordType>& v0,
                              const Point3<CoordType>& v1) {
  using Wide = WideInt<3>;
  const auto cross_full = qcross<CoordType>(v0, v1);
  return normalizeByBitShift<CoordType>(cross_full);
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
// interior half-plane, inside the exterior half-plane, or colinear/coplanar.
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
  auto ztest = qdot<3, CoordType>(ap, n);
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
  for (auto i = 1; i < NP; ++i) {
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
bool SAT(const std::vector<Point3<CoordType>>& pointsA,
         const std::vector<Point3<CoordType>>& pointsB,
         const Point3<CoordType>& axis) {
  using Wide = WideInt<3>;
  Wide minA = HashKey<3>::hashMax();
  Wide minB = minA, maxA = -minA, maxB = maxA;
  for (const auto& p : pointsA) {
    auto ztest = qdot<3, CoordType>(p, axis);
    if (ztest < minA) {
      minA = ztest;
    }
    if (ztest > maxA) {
      maxA = ztest;
    }
  }
  for (const auto& p : pointsB) {
    auto ztest = qdot<3, CoordType>(p, axis);
    if (ztest < minB) {
      minB = ztest;
    }
    if (ztest > maxB) {
      maxB = ztest;
    }
  }
  return maxA < minB || maxB < minA;
}
}
#endif
