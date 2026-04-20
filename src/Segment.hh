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
#include "DimensionTraits.hh"

namespace polytope {

template<int nDim>
struct Segment {
  using Point = DimensionTraits<2, double>::IntPoint;
  typedef DimensionTraits<2, double>::CoordHash CoordHash;
  Point a, b;
  Segment(): a(), b() {}
  Segment(CoordHash x0, CoordHash y0, CoordHash x1, CoordHash y1): a(x0,y0), b(x1,y1) {}
  bool operator==(const Segment& rhs) const { return a == rhs.a and b == rhs.b; }
  bool operator!=(const Segment& rhs) const { return !(*this == rhs); }
  bool operator<(const Segment& rhs) const {
    return (a < rhs.a                ? true :
            a == rhs.a and b < rhs.b ? true :
            false);
  }
};

// It's nice being able to print these things.
template<int nDim>
std::ostream&
operator<<(std::ostream& os, const Segment<nDim>& s) {
  os << "Segment(" << s.a << " " << s.b << ")";
  return os;
}

// Serialization.
template<int nDim>
struct Serializer<Segment<nDim> > {

  static void serializeImpl(const Segment<nDim>& value,
                            std::vector<char>& buffer) {
    serialize(value.a, buffer);
    serialize(value.b, buffer);
  }

  static void deserializeImpl(Segment<nDim>& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.a, bufItr, endItr);
    deserialize(value.b, bufItr, endItr);
  }
};

}

#endif
