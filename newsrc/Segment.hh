//----------------------------------------------------------------------------//
// 2D and 3D integral Segment types used internally in polytope.  Not really
// for external consumption!.
//----------------------------------------------------------------------------//
#ifndef __polytope_Segment__
#define __polytope_Segment__

#include <iostream>
#include <iterator>

#include "polytope_serialize.hh"
#include "polytope_internal.hh"
#include "Point.hh"

namespace polytope {

//template<int Dimension, typename CoordType>
struct Segment {
  using CoordType = HashKey<2>::IntType;
  using PointType = Point<2, CoordType>;
  PointType a, b;
  Segment(): a(), b() {}
  Segment(CoordType x0, CoordType y0, CoordType x1, CoordType y1): a(x0,y0), b(x1,y1) {}
  bool operator==(const Segment& rhs) const { return a == rhs.a and b == rhs.b; }
  bool operator!=(const Segment& rhs) const { return !(*this == rhs); }
  bool operator<(const Segment& rhs) const {
    return (a < rhs.a                ? true :
            a == rhs.a and b < rhs.b ? true :
            false);
  }
};

// It's nice being able to print these things.
std::ostream&
operator<<(std::ostream& os, const Segment& s) {
  os << "Segment(" << s.a << " " << s.b << ")";
  return os;
}

// Serialization.
template<>
struct Serializer<Segment> {

  static void serializeImpl(const Segment& value,
                            std::vector<char>& buffer) {
    serialize(value.a, buffer);
    serialize(value.b, buffer);
  }

  static void deserializeImpl(Segment& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.a, bufItr, endItr);
    deserialize(value.b, bufItr, endItr);
  }
};

}

#endif
