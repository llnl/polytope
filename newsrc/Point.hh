//----------------------------------------------------------------------------//
// 2D and 3D integral Point types used internally in polytope.  Not really
// for external consumption!.
//----------------------------------------------------------------------------//
#ifndef __polytope_Point__
#define __polytope_Point__

#include "polytope_serialize.hh"
#include "polytope_internal.hh"

#include <iostream>
#include <iterator>
#include <cmath>

namespace polytope {

namespace {
inline std::ostream& operator<<(std::ostream& os, __int128 n) {
  if (n == 0) return os << "0";
  if (n < 0) {
    os << "-";
    n = -n;
  }
  std::string s;
  while (n > 0) {
    s += (char)('0' + (n % 10));
    n /= 10;
  }
  std::reverse(s.begin(), s.end());
  return os << s;
}
}

//------------------------------------------------------------------------------
// A integer version of the simple 2D point.
//------------------------------------------------------------------------------

template<int Dimension, typename CoordType> class Point {};

template<typename CoordType>
using Point2 = Point<2, CoordType>;

template<typename CoordType>
using Point3 = Point<3, CoordType>;

template<typename CoordType>
struct Point<2, CoordType> {
  CoordType x, y;
  unsigned index;
  // Constructors
  Point(): x(0), y(0), index(0) {}

  Point(const CoordType xy) :
    x(xy), y(xy), index(0) {}

  Point(const CoordType xi,
        const CoordType yi,
        const unsigned i = 0) :
    x(xi), y(yi), index(i) {}

  Point(const CoordType* ri,
        const unsigned i = 0) :
    x(ri[0]), y(ri[1]), index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType& dx,
        const unsigned i = 0):
    x(static_cast<CoordType>(ri[0]/dx + 0.5)),
    y(static_cast<CoordType>(ri[1]/dx + 0.5)),
    index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType* dx,
        const unsigned i = 0):
    x(static_cast<CoordType>(ri[0]/dx[0] + 0.5)),
    y(static_cast<CoordType>(ri[1]/dx[1] + 0.5)),
    index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType* rlow,
        const RealType* dx, const unsigned i = 0):
    x(static_cast<CoordType>((ri[0] - rlow[0])/dx[0] + 0.5)),
    y(static_cast<CoordType>((ri[1] - rlow[1])/dx[1] + 0.5)),
    index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType* rlow,
        const RealType& dx, const unsigned i = 0):
    x(static_cast<CoordType>((ri[0] - rlow[0])/dx + 0.5)),
    y(static_cast<CoordType>((ri[1] - rlow[1])/dx + 0.5)),
    index(i) {}

  // Operators
  bool operator==(const Point<2, CoordType>& rhs) const { return (x == rhs.x and y == rhs.y); }
  bool operator!=(const Point<2, CoordType>& rhs) const { return !(*this == rhs); }
  // Element-wise comparison operators
  bool allLess(const Point& rhs) const {
    return (x < rhs.x && y < rhs.y);
  }
  bool allLessEqual(const Point& rhs) const {
    return (x <= rhs.x && y <= rhs.y);
  }
  bool allGreater(const Point& rhs) const {
    return (x > rhs.x && y > rhs.y);
  }
  bool allGreaterEqual(const Point& rhs) const {
    return (x >= rhs.x && y >= rhs.y);
  }
  // NOTE: Comparison operators are lexicographic, not element-wise. Use above for element-wise.
  bool operator<(const Point& rhs) const {
    return (x < rhs.x || (x == rhs.x && y < rhs.y));
  }
  bool operator<=(const Point& rhs) const {
    return (x < rhs.x || (x == rhs.x && y <= rhs.y));
  }
  bool operator>(const Point& rhs) const {
    return (x > rhs.x || (x == rhs.x && y > rhs.y));
  }
  bool operator>=(const Point& rhs) const {
    return (x > rhs.x || (x == rhs.x && y >= rhs.y));
  }

  Point& operator+=(const Point& rhs) { x += rhs.x; y += rhs.y; return *this; }
  Point& operator-=(const Point& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
  Point& operator*=(const CoordType& rhs) { x *= rhs; y *= rhs; return *this; }
  Point& operator/=(const CoordType& rhs) { x /= rhs; y /= rhs; return *this; }
  Point& operator/=(const Point& rhs) { x /= rhs.x; y /= rhs.y; return *this; }
  Point operator+(const Point& rhs) const { Point result(*this); result += rhs; return result; }
  Point operator-(const Point& rhs) const { Point result(*this); result -= rhs; return result; }
  Point operator*(const Point& rhs) const { return Point(x*rhs.x, y*rhs.y); }
  Point operator*(const CoordType& rhs) const { Point result(*this); result *= rhs; return result; }
  Point operator/(const CoordType& rhs) const { Point result(*this); result /= rhs; return result; }
  Point operator/(const Point& rhs) const { return Point(x/rhs.x, y/rhs.y); }
  Point operator-() const { return Point(-x, -y); }
  CoordType  operator[](const size_t i) const { POLY_ASSERT(i < 2); return *(&x + i); }
  CoordType& operator[](const size_t i)       { POLY_ASSERT(i < 2); return *(&x + i); }

  void clipPoint(const Point& lorhs, const Point& hirhs) {
    x = std::min(hirhs.x, std::max(lorhs.x, x));
    y = std::min(hirhs.y, std::max(lorhs.y, y));
  }

  bool iszero() const { return (x == 0 && y == 0) ? true : false; }
  void zero() { x = 0; y = 0; }
  void one() { x = 1; y = 1; }

  template<typename IntType, typename RealType>
  Point<2, IntType> convertXi(const Point<2, RealType>& blo,
                              const Point<2, RealType>& dx) const {
    // Quantize: RealType -> IntType
    POLY_ASSERT(typeid(CoordType) == typeid(RealType));
    IntType xOut, yOut;
    xOut = static_cast<IntType>((this->x - blo.x)/dx.x + 0.5);
    yOut = static_cast<IntType>((this->y - blo.y)/dx.y + 0.5);
    return Point<2, IntType>(xOut, yOut);
  }

  template<typename RealType>
  Point<2, RealType> convertx(const Point<2, RealType>& blo,
                              const Point<2, RealType>& dx) const {
    //POLY_ASSERT(typeid(CoordType) != typeid(RealType));
    RealType xOut, yOut;
    // Dequantize: IntType -> RealType
    xOut = dx.x*(static_cast<RealType>(this->x) - 0.5) + blo.x;
    yOut = dx.y*(static_cast<RealType>(this->y) - 0.5) + blo.y;
    return Point<2, RealType>(xOut, yOut);
  }

  template<typename RealType>
  Point<2, RealType> type_cast() const {
    return Point<2, RealType>(static_cast<RealType>(x), static_cast<RealType>(y));
  }

  // Return the min and max elements in each direction
  Point minElements(const Point& in) const {
    CoordType xOut = (in.x < this->x) ? in.x : this->x;
    CoordType yOut = (in.y < this->y) ? in.y : this->y;
    return Point(xOut, yOut);
  }

  Point maxElements(const Point& in) const {
    CoordType xOut = (in.x > this->x) ? in.x : this->x;
    CoordType yOut = (in.y > this->y) ? in.y : this->y;
    return Point(xOut, yOut);
  }

  int maxAxis() const {
    return (x >= y) ? 0 : 1;
  }

  template<typename IntType>
  Point<2, IntType> bitShift(const int shift) const {
    return Point<2, IntType>(
      static_cast<IntType>(x >> shift),
      static_cast<IntType>(y >> shift)
    );
  }
};

// It's nice being able to print these things.
template<typename CoordType>
std::ostream&
operator<<(std::ostream& os, const Point<2, CoordType>& p) {
  os << "[" << p.x << ", " << p.y << "]";//(" << p.index << ")";
  return os;
}

template<typename CoordType>
std::ostream&
operator<<(std::ostream& os, const std::vector<Point<2, CoordType>>& pv) {
  os << "v = [";
  for (const auto& p : pv) {
    os << p << "," << std::endl;
  }
  os << "]";
  return os;
}

template<typename CoordType, typename CoordHash>
CoordHash dot(const Point2<CoordType>& a,
              const Point2<CoordType>& b) {
  return (a.x*b.x) + (a.y*b.y);
}

// Serialization.
template<typename CoordType>
struct Serializer<Point<2, CoordType> > {

  static void serializeImpl(const Point<2, CoordType>& value,
                            std::vector<char>& buffer) {
    serialize(value.x, buffer);
    serialize(value.y, buffer);
    serialize(value.index, buffer);
  }

  static void deserializeImpl(Point<2, CoordType>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.x, bufItr, endItr);
    deserialize(value.y, bufItr, endItr);
    deserialize(value.index, bufItr, endItr);
  }
};

//------------------------------------------------------------------------------
// A integer version of the simple 3D point.
//------------------------------------------------------------------------------
template<typename CoordType>
struct Point<3, CoordType> {
  CoordType x, y, z;
  unsigned index;
  // Constructors
  Point(): x(0), y(0), z(0), index(0) {}
  Point(const CoordType xyz) :
    x(xyz), y(xyz), z(xyz), index(0) {}
  Point(const CoordType xi,
        const CoordType yi,
        const CoordType zi,
        const unsigned i = 0) :
    x(xi), y(yi), z(zi), index(i) {}

  Point(const CoordType* ri,
        const unsigned i = 0) :
    x(ri[0]), y(ri[1]), z(ri[2]), index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType* dx,
        const unsigned i = 0):
    x(static_cast<CoordType>(ri[0]/dx[0] + 0.5)),
    y(static_cast<CoordType>(ri[1]/dx[1] + 0.5)),
    z(static_cast<CoordType>(ri[2]/dx[2] + 0.5)),
    index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType& dx,
        const unsigned i = 0):
    x(static_cast<CoordType>(ri[0]/dx + 0.5)),
    y(static_cast<CoordType>(ri[1]/dx + 0.5)),
    z(static_cast<CoordType>(ri[2]/dx + 0.5)),
    index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType* rlow,
        const RealType* dx, const unsigned i = 0):
    x(static_cast<CoordType>((ri[0] - rlow[0])/dx[0] + 0.5)),
    y(static_cast<CoordType>((ri[1] - rlow[1])/dx[1] + 0.5)),
    z(static_cast<CoordType>((ri[2] - rlow[2])/dx[2] + 0.5)),
    index(i) {}

  template<typename RealType>
  Point(const RealType* ri, const RealType* rlow,
        const RealType& dx, const unsigned i = 0):
    x(static_cast<CoordType>((ri[0] - rlow[0])/dx + 0.5)),
    y(static_cast<CoordType>((ri[1] - rlow[1])/dx + 0.5)),
    z(static_cast<CoordType>((ri[2] - rlow[2])/dx + 0.5)),
    index(i) {}

  // Operators
  bool operator==(const Point& rhs) const { return (x == rhs.x and y == rhs.y and z == rhs.z); }
  bool operator!=(const Point& rhs) const { return !(*this == rhs); }
  // Element-wise comparison operators
  bool allLess(const Point& rhs) const {
    return (x < rhs.x && y < rhs.y && z < rhs.z);
  }
  bool allLessEqual(const Point& rhs) const {
    return (x <= rhs.x && y <= rhs.y && z <= rhs.z);
  }
  bool allGreater(const Point& rhs) const {
    return (x > rhs.x && y > rhs.y && z > rhs.z);
  }
  bool allGreaterEqual(const Point& rhs) const {
    return (x >= rhs.x && y >= rhs.y && z >= rhs.z);
  }
  // NOTE: Comparison operators are lexicographic, not element-wise. Use above for element-wise.
  bool operator<(const Point& rhs) const {
    return (x < rhs.x || (x == rhs.x && (y < rhs.y || (y == rhs.y && z < rhs.z))));
  }
  bool operator<=(const Point& rhs) const {
    return (x < rhs.x || (x == rhs.x && (y < rhs.y || (y == rhs.y && z <= rhs.z))));
  }
  bool operator>(const Point& rhs) const {
    return (x > rhs.x || (x == rhs.x && (y > rhs.y || (y == rhs.y && z > rhs.z))));
  }
  bool operator>=(const Point& rhs) const {
    return (x >= rhs.x || (x == rhs.x && (y > rhs.y || (y == rhs.y && z >= rhs.z))));
  }

  template<typename RealType>
  Point(const RealType& xi, const RealType& yi, const RealType& zi,
        const RealType& xlow, const RealType& ylow, const RealType& zlow,
        const RealType& dx,
        const unsigned i = 0):
    x(static_cast<CoordType>((xi - xlow)/dx + 0.5)),
    y(static_cast<CoordType>((yi - ylow)/dx + 0.5)),
    z(static_cast<CoordType>((zi - zlow)/dx + 0.5)),
    index(i) {}

  Point& operator+=(const Point& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
  Point& operator-=(const Point& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
  Point& operator*=(const CoordType& rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
  Point& operator/=(const CoordType& rhs) { x /= rhs; y /= rhs; z /= rhs; return *this; }
  Point& operator/=(const Point& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }
  Point operator+(const Point& rhs) const { Point result(*this); result += rhs; return result; }
  Point operator-(const Point& rhs) const { Point result(*this); result -= rhs; return result; }
  Point operator*(const Point& rhs) const { return Point(x*rhs.x, y*rhs.y, z*rhs.z); }
  Point operator*(const CoordType& rhs) const { Point result(*this); result *= rhs; return result; }
  Point operator/(const CoordType& rhs) const { Point result(*this); result /= rhs; return result; }
  Point operator/(const Point& rhs) const { return Point(x/rhs.x, y/rhs.y, z/rhs.z); }
  Point operator-() const { return Point(-x, -y, -z); }
  CoordType  operator[](const size_t i) const { POLY_ASSERT(i < 3); return *(&x + i); }
  CoordType& operator[](const size_t i)       { POLY_ASSERT(i < 3); return *(&x + i); }

  void clipPoint(const Point& lorhs, const Point& hirhs) {
    x = std::min(hirhs.x, std::max(lorhs.x, x));
    y = std::min(hirhs.y, std::max(lorhs.y, y));
    z = std::min(hirhs.z, std::max(lorhs.z, z));
  }

  bool iszero() const { return (x == 0 && y == 0 && z == 0) ? true : false; }
  void zero() { x = 0; y = 0; z = 0; }
  void one() { x = 1; y = 1; z = 1; }

  template<typename IntType, typename RealType>
  Point<3, IntType> convertXi(const Point<3, RealType>& blo,
                              const Point<3, RealType>& dx) const {
    // Quantize: RealType -> IntType
    POLY_ASSERT(typeid(CoordType) == typeid(RealType));
    IntType xOut, yOut, zOut;
    xOut = static_cast<IntType>((this->x - blo.x)/dx.x + 0.5);
    yOut = static_cast<IntType>((this->y - blo.y)/dx.y + 0.5);
    zOut = static_cast<IntType>((this->z - blo.z)/dx.z + 0.5);
    return Point<3, IntType>(xOut, yOut, zOut, index);
  }

  template<typename RealType>
  Point<3, RealType> convertx(const Point<3, RealType>& blo,
                              const Point<3, RealType>& dx) const {
    POLY_ASSERT(typeid(CoordType) != typeid(RealType));
    RealType xOut, yOut, zOut;
    // Dequantize: IntType -> RealType
    xOut = dx.x*(static_cast<RealType>(this->x) - 0.5) + blo.x;
    yOut = dx.y*(static_cast<RealType>(this->y) - 0.5) + blo.y;
    zOut = dx.z*(static_cast<RealType>(this->z) - 0.5) + blo.z;
    return Point<3, RealType>(xOut, yOut, zOut, index);
  }

  template<typename RealType>
  Point<3, RealType> type_cast() const {
    return Point<3, RealType>(static_cast<RealType>(x),
                              static_cast<RealType>(y),
                              static_cast<RealType>(z));                            
  }

  // Return the min and max elements in each direction
  Point minElements(const Point& in) const {
    CoordType xOut = (in.x < this->x) ? in.x : this->x;
    CoordType yOut = (in.y < this->y) ? in.y : this->y;
    CoordType zOut = (in.z < this->z) ? in.z : this->z;
    return Point(xOut, yOut, zOut);
  }

  Point maxElements(const Point& in) const {
    CoordType xOut = (in.x > this->x) ? in.x : this->x;
    CoordType yOut = (in.y > this->y) ? in.y : this->y;
    CoordType zOut = (in.z > this->z) ? in.z : this->z;
    return Point(xOut, yOut, zOut);
  }

  int maxAxis() const {
    if (x >= y && x >= z) {
      return 0;
    } else if (y >= z) {
      return 1;
    }
    return 2;
  }

  template<typename IntType>
  Point<3, IntType> bitShift(const int shift) const {
    return Point<3, IntType>(
      static_cast<IntType>(x >> shift),
      static_cast<IntType>(y >> shift),
      static_cast<IntType>(z >> shift)
    );
  }
};

template<typename CoordType, typename CoordHash>
CoordHash dot(const Point3<CoordType>& a, const Point3<CoordType>& b) {
  return (a.x*b.x) + (a.y*b.y);
}

// It's nice being able to print these things.
template<typename CoordType>
std::ostream&
operator<<(std::ostream& os, const Point<3, CoordType>& p) {
  os << "(" << p.x << " " << p.y << " " << p.z <<  ")(" << p.index << ")";
  return os;
}

// Serialization.
template<typename CoordType>
struct Serializer<Point<3, CoordType> > {

  static void serializeImpl(const Point<3, CoordType>& value,
                            std::vector<char>& buffer) {
    serialize(value.x, buffer);
    serialize(value.y, buffer);
    serialize(value.z, buffer);
    serialize(value.index, buffer);
  }

  static void deserializeImpl(Point<3, CoordType>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.x, bufItr, endItr);
    deserialize(value.y, bufItr, endItr);
    deserialize(value.z, bufItr, endItr);
    deserialize(value.index, bufItr, endItr);
  }
};

//------------------------------------------------------------------------------
// Provide a special comparator for point types with some fuzz.
//------------------------------------------------------------------------------
template<typename CoordType>
struct PointComparator {
  CoordType mfuzz;
  PointComparator(const CoordType fuzz): mfuzz(fuzz) {}
  bool operator()(const Point<2, CoordType>& lhs, const Point<2, CoordType>& rhs) const {
    return (rhs.x - lhs.x > mfuzz                                      ? true :
            std::abs(rhs.x - lhs.x) <= mfuzz and rhs.y - lhs.y > mfuzz ? true :
            false);
  }
  bool operator()(const Point<3, CoordType>& lhs, const Point<3, CoordType>& rhs) const {
    return (rhs.x - lhs.x > mfuzz                                                                           ? true :
            std::abs(rhs.x - lhs.x) <= mfuzz and rhs.y - lhs.y > mfuzz                                      ? true :
            std::abs(rhs.x - lhs.x) <= mfuzz and std::abs(rhs.y - lhs.y) <= mfuzz and rhs.z - lhs.z > mfuzz ? true :
            false);
  }
};

// General functions
template<int Dimension, typename CoordType>
inline
Point<Dimension, CoordType>
operator*(const CoordType val, const Point<Dimension, CoordType>& vec) {
  return vec*val;
}

// Roll flattened coordinates into Points
template<int Dimension, typename CoordType>
std::vector<Point<Dimension, CoordType>>
extractCoords(const std::vector<CoordType>& allpoints) {
  auto n = allpoints.size()/Dimension;
  std::vector<Point<Dimension, CoordType>> result(n);
  for(auto i = 0; i < n; ++i) {
    result[i] = Point<Dimension, CoordType>(&(allpoints[Dimension*i]), i);
  }
  return result;
}

template<int Dimension, typename CoordType>
std::vector<CoordType>
flattenCoords(const std::vector<Point<Dimension, CoordType>>& allpoints) {
  auto n = allpoints.size();
  auto n2 = Dimension*n;
  std::vector<CoordType> result(n2);
  for(auto i = 0; i < n; ++i) {
    for(auto d = 0; d < Dimension; ++d) {
      result[Dimension*i+d] = allpoints[i][d];
    }
  }
  return result;
}

template<int Dimension, typename CoordType>
void
findBoundingElements(const std::vector<Point<Dimension, CoordType>>& allpoints,
                     Point<Dimension, CoordType>& minPoint,
                     Point<Dimension, CoordType>& maxPoint) {
  for (const auto& p : allpoints) {
    minPoint = minPoint.minElements(p);
    maxPoint = maxPoint.maxElements(p);
  }
}

template<int Dimension, typename CoordType>
Point<Dimension, CoordType> round(const Point<Dimension, double>& point) {
  Point<Dimension, CoordType> out;
  for (int d = 0; d < Dimension; ++d) {
    out[d] = static_cast<CoordType>(std::floor(point[d]));
  }
  return out;
}

} // namespace polytope
#endif
