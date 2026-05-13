//----------------------------------------------------------------------------//
// 2D and 3D integral Point types used internally in polytope.  Not really
// for external consumption!.
//----------------------------------------------------------------------------//
#ifndef __polytope_Point__
#define __polytope_Point__

#include <iostream>
#include <iterator>

#include "polytope_serialize.hh"
#include "polytope_internal.hh"

namespace polytope {

//------------------------------------------------------------------------------
// A integer version of the simple 2D point.
//------------------------------------------------------------------------------
template<typename CoordType>
struct Point2 {
  CoordType x, y;
  unsigned index;
  Point2(): x(0), y(0), index(0) {}

  Point2(const CoordType& xi,
         const CoordType& yi,
         const unsigned i = 0) :
    x(xi), y(yi), index(i) {}

  Point2& operator=(const Point2& rhs) {
    x = rhs.x;
    y = rhs.y;
    index = rhs.index;
    return *this;
  }
  bool operator==(const Point2& rhs) const { return (x == rhs.x and y == rhs.y); }
  bool operator!=(const Point2& rhs) const { return !(*this == rhs); }
  bool operator<(const Point2& rhs) const {
    return (x < rhs.x                ? true :
            x == rhs.x and y < rhs.y ? true :
            false);
  }
  bool operator>(const Point2& rhs) const {
    return (x > rhs.x                ? true :
            x == rhs.x and y > rhs.y ? true :
            false);
  }

  template<typename RealType>
  Point2(const Point2<RealType>& real_point,
         const Point2<RealType>& dx,
         const unsigned i = 0):
    x(static_cast<CoordType>(real_point.x/dx.x + 0.5)),
    y(static_cast<CoordType>(real_point.y/dx.y + 0.5)),
    index(i) {}

  template<typename RealType>
  Point2(const Point2<RealType>& real_point,
         const Point2<RealType>& blo,
         const Point2<RealType>& dx,
         const unsigned i = 0):
    x(static_cast<CoordType>((real_point.x - blo.x)/dx.x + 0.5)),
    y(static_cast<CoordType>((real_point.y - blo.y)/dx.y + 0.5)),
    index(i) {}

  template<typename RealType>
  Point2(const RealType& xi, const RealType& yi, const RealType& dx, const unsigned i = 0): 
    x(static_cast<CoordType>(xi/dx + 0.5)),
    y(static_cast<CoordType>(yi/dx + 0.5)),
    index(i) {}

  template<typename RealType>
  Point2(const RealType& xi, const RealType& yi, 
         const RealType& xlow, const RealType& ylow,
         const RealType& dx,
         const unsigned i = 0): 
    x(static_cast<CoordType>((xi - xlow)/dx + 0.5)),
    y(static_cast<CoordType>((yi - ylow)/dx + 0.5)),
    index(i) {}

  template<typename RealType> RealType realx(const RealType& xmin, const RealType& dx) const { return static_cast<RealType>(x*dx) + xmin; }
  template<typename RealType> RealType realy(const RealType& ymin, const RealType& dy) const { return static_cast<RealType>(y*dy) + ymin; }

  template<typename RealType> Point2<RealType>
  realPoint(const Point2<RealType>& dx) const {
    Point2<RealType> real_point;
    real_point.x = static_cast<RealType>(dx.x*(x - 0.5));
    real_point.y = static_cast<RealType>(dx.y*(y - 0.5));
    return real_point;
  }

  template<typename RealType> Point2<RealType>
  realPoint(const Point2<RealType>& blo,
            const Point2<RealType>& dx) const {
    Point2<RealType> real_point;
    real_point.x = static_cast<RealType>(dx.x*(x - 0.5) + blo.x);
    real_point.y = static_cast<RealType>(dx.y*(y - 0.5) + blo.y);
    return real_point;
  }

  Point2& operator+=(const Point2& rhs) { x += rhs.x; y += rhs.y; return *this; }
  Point2& operator-=(const Point2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
  Point2& operator*=(const CoordType& rhs) { x *= rhs; y *= rhs; return *this; }
  Point2& operator/=(const CoordType& rhs) { x /= rhs; y /= rhs; return *this; }
  Point2 operator+(const Point2& rhs) const { Point2 result(*this); result += rhs; return result; }
  Point2 operator-(const Point2& rhs) const { Point2 result(*this); result -= rhs; return result; }
  Point2 operator*(const CoordType& rhs) const { Point2 result(*this); result *= rhs; return result; }
  Point2 operator/(const CoordType& rhs) const { Point2 result(*this); result /= rhs; return result; }
  Point2 operator-() const { return Point2(-x, -y); }
  CoordType  operator[](const size_t i) const { POLY_ASSERT(i < 2); return *(&x + i); }
  CoordType& operator[](const size_t i)       { POLY_ASSERT(i < 2); return *(&x + i); }
};

// It's nice being able to print these things.
template<typename CoordType>
std::ostream&
operator<<(std::ostream& os, const Point2<CoordType>& p) {
  os << "(" << p.x << " " << p.y << ")(" << p.index << ")";
  return os;
}

// Serialization.
template<typename CoordType>
struct Serializer<Point2<CoordType> > {

  static void serializeImpl(const Point2<CoordType>& value,
                            std::vector<char>& buffer) {
    serialize(value.x, buffer);
    serialize(value.y, buffer);
    serialize(value.index, buffer);
  }

  static void deserializeImpl(Point2<CoordType>& value,
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
struct Point3 {
  CoordType x, y, z;
  unsigned index;
  Point3(): x(0), y(0), z(0), index(0) {}
  Point3(const CoordType& xi, const CoordType& yi, const CoordType& zi, const unsigned i = 0): x(xi), y(yi), z(zi), index(i) {}
  Point3& operator=(const Point3& rhs) {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    index = rhs.index;
    return *this;
  }
  bool operator==(const Point3& rhs) const { return (x == rhs.x and y == rhs.y and z == rhs.z); }
  bool operator!=(const Point3& rhs) const { return !(*this == rhs); }
  bool operator<(const Point3& rhs) const {
    return (x < rhs.x                               ? true :
            x == rhs.x and y < rhs.y                ? true :
            x == rhs.x and y == rhs.y and z < rhs.z ? true :
            false);
  }
  bool operator>(const Point3& rhs) const {
    return (x > rhs.x                               ? true :
            x == rhs.x and y > rhs.y                ? true :
            x == rhs.x and y == rhs.y and z > rhs.z ? true :
            false);
  }

  template<typename RealType>
  Point3(const Point3<RealType>& real_point,
         const Point3<RealType>& dx,
         const unsigned i = 0):
    x(static_cast<CoordType>(real_point.x/dx.x + 0.5)),
    y(static_cast<CoordType>(real_point.y/dx.y + 0.5)),
    z(static_cast<CoordType>(real_point.z/dx.z + 0.5)),
    index(i) {}

  template<typename RealType>
  Point3(const Point3<RealType>& real_point,
         const Point3<RealType>& blo,
         const Point3<RealType>& dx,
         const unsigned i = 0):
    x(static_cast<CoordType>((real_point.x - blo.x)/dx.x + 0.5)),
    y(static_cast<CoordType>((real_point.y - blo.y)/dx.y + 0.5)),
    z(static_cast<CoordType>((real_point.z - blo.z)/dx.z + 0.5)),
    index(i) {}

  template<typename RealType>
  Point3(const RealType& xi, const RealType& yi, const RealType& zi, const RealType& dx, const unsigned i = 0): 
    x(static_cast<CoordType>(xi/dx + 0.5)),
    y(static_cast<CoordType>(yi/dx + 0.5)),
    z(static_cast<CoordType>(zi/dx + 0.5)),
    index(i) {}
  template<typename RealType>
  Point3(const RealType& xi, const RealType& yi, const RealType& zi, 
         const RealType& xlow, const RealType& ylow, const RealType& zlow, 
         const RealType& dx, 
         const unsigned i = 0): 
    x(static_cast<CoordType>((xi - xlow)/dx + 0.5)),
    y(static_cast<CoordType>((yi - ylow)/dx + 0.5)),
    z(static_cast<CoordType>((zi - zlow)/dx + 0.5)),
    index(i) {}
  template<typename RealType> RealType realx(const RealType& xmin, const RealType& dx) const { return static_cast<RealType>(x*dx) + xmin; }
  template<typename RealType> RealType realy(const RealType& ymin, const RealType& dy) const { return static_cast<RealType>(y*dy) + ymin; }
  template<typename RealType> RealType realz(const RealType& zmin, const RealType& dz) const { return static_cast<RealType>(z*dz) + zmin; }

  template<typename RealType> Point3<RealType>
  realPoint(const Point3<RealType>& dx) const {
    Point3<RealType> real_point;
    real_point.x = static_cast<RealType>(dx.x*(x - 0.5));
    real_point.y = static_cast<RealType>(dx.y*(y - 0.5));
    real_point.z = static_cast<RealType>(dx.z*(z - 0.5));
    return real_point;
  }

  template<typename RealType> Point3<RealType>
  realPoint(const Point3<RealType>& blo,
            const Point3<RealType>& dx) const {
    Point3<RealType> real_point;
    real_point.x = static_cast<RealType>(dx.x*(x - 0.5) + blo.x);
    real_point.y = static_cast<RealType>(dx.y*(y - 0.5) + blo.y);
    real_point.z = static_cast<RealType>(dx.z*(z - 0.5) + blo.z);
    return real_point;
  }

  Point3& operator+=(const Point3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
  Point3& operator-=(const Point3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
  Point3& operator*=(const CoordType& rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
  Point3& operator/=(const CoordType& rhs) { x /= rhs; y /= rhs; z /= rhs; return *this; }
  Point3 operator+(const Point3& rhs) const { Point3 result(*this); result += rhs; return result; }
  Point3 operator-(const Point3& rhs) const { Point3 result(*this); result -= rhs; return result; }
  Point3 operator*(const CoordType& rhs) const { Point3 result(*this); result *= rhs; return result; }
  Point3 operator/(const CoordType& rhs) const { Point3 result(*this); result /= rhs; return result; }
  Point3 operator-() const { return Point3(-x, -y, -z); }
  CoordType  operator[](const size_t i) const { POLY_ASSERT(i < 3); return *(&x + i); }
  CoordType& operator[](const size_t i)       { POLY_ASSERT(i < 3); return *(&x + i); }
};

// It's nice being able to print these things.
template<typename CoordType>
std::ostream&
operator<<(std::ostream& os, const Point3<CoordType>& p) {
  os << "(" << p.x << " " << p.y << " " << p.z <<  ")(" << p.index << ")";
  return os;
}

// Serialization.
template<typename CoordType>
struct Serializer<Point3<CoordType> > {

  static void serializeImpl(const Point3<CoordType>& value,
                            std::vector<char>& buffer) {
    serialize(value.x, buffer);
    serialize(value.y, buffer);
    serialize(value.z, buffer);
    serialize(value.index, buffer);
  }

  static void deserializeImpl(Point3<CoordType>& value,
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
  bool operator()(const Point2<CoordType>& lhs, const Point2<CoordType>& rhs) const {
    return (rhs.x - lhs.x > mfuzz                                      ? true :
            std::abs(rhs.x - lhs.x) <= mfuzz and rhs.y - lhs.y > mfuzz ? true :
            false);
  }
  bool operator()(const Point3<CoordType>& lhs, const Point3<CoordType>& rhs) const {
    return (rhs.x - lhs.x > mfuzz                                                                           ? true :
            std::abs(rhs.x - lhs.x) <= mfuzz and rhs.y - lhs.y > mfuzz                                      ? true :
            std::abs(rhs.x - lhs.x) <= mfuzz and std::abs(rhs.y - lhs.y) <= mfuzz and rhs.z - lhs.z > mfuzz ? true :
            false);
  }
};

}

#endif
