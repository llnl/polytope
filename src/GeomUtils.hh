#ifndef __Polytope_GeomUtils__
#define __Polytope_GeomUtils__

#include "polytope.hh"
#include "Point.hh"
#include "MortonKeyTraits.hh"
#include "EdgeUtils.hh"
#include "Quantizer.hh"

namespace polytope {

template<int Dimension>
using WideInt = typename MortonKeyTraits<Dimension>::Wide;

template<int Dimension>
using BigInt = typename MortonKeyTraits<Dimension>::Big;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// General helper routines
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//------------------------------------------------------------------------------
// Quantized dot product with overflow protection
// Casts to WideInt before multiplication to prevent overflow
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
WideInt<Dimension> qdot(const Point<Dimension, CoordType>& a,
                        const Point<Dimension, CoordType>& b) {
  using Wide = WideInt<Dimension>;
  const auto ah = a.template type_cast<Wide>();
  const auto bh = b.template type_cast<Wide>();

  Wide sum = 0;
  for (int i = 0; i < Dimension; ++i) {
    sum += ah[i] * bh[i];
  }
  return sum;
}

//------------------------------------------------------------------------------
// Larger quantized dot product with overflow protection
// Casts to BigInt before multiplication to prevent overflow
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType, typename Wide>
BigInt<Dimension> qqdot(const Point<Dimension, Wide>& a,
                        const Point<Dimension, CoordType>& b) {
  using Big = BigInt<Dimension>;
  const auto ah = a.template type_cast<Big>();
  const auto bh = b.template type_cast<Big>();

  Big sum = 0;
  for (int i = 0; i < Dimension; ++i) {
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
// Compare the magnitude of 2 points. Return true if p1 > p2
//------------------------------------------------------------------------------
template<int Dimension, typename CoordType>
bool magComparison(const Point<Dimension, CoordType>& p1,
                   const Point<Dimension, CoordType>& p2) {
  using Wide = WideInt<Dimension>;
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
  const double SCALE = std::pow(2.0, MortonKeyTraits<2>::bitsPerCoordinate - 2);
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
  using Wide = WideInt<2>;
  Point2<Wide> sum = gen0.template type_cast<Wide>()
    + gen1.template type_cast<Wide>();
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
         const Point2<WideInt<2>>& axis) {
  using Wide = WideInt<2>;
  using Big = BigInt<2>;
  POLY_ASSERT(!pointsA.empty());
  POLY_ASSERT(!pointsB.empty());
  if (axis.iszero()) return false;

  Big minA = qqdot(axis, pointsA.front());
  Big maxA = minA;
  for (const auto& p : pointsA) {
    const auto ztest = qqdot(axis, p);
    if (ztest < minA) {
      minA = ztest;
    }
    if (ztest > maxA) {
      maxA = ztest;
    }
  }

  Big minB = qqdot(axis, pointsB.front());
  Big maxB = minB;
  for (const auto& p : pointsB) {
    const auto ztest = qqdot(axis, p);
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
  using Wide = WideInt<3>;
  Wide minA = MortonKeyTraits<3>::maxKey();
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
