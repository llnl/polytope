//------------------------------------------------------------------------------
// This file is used to register our types with the Boost.Polygon library.
//------------------------------------------------------------------------------
#ifndef __polytope_RegisterBoostPolygonType__
#define __polytope_RegisterBoostPolygonType__

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>

// Include Boost multiprecision types
#include <boost/multiprecision/cpp_int.hpp>

#include "Point.hh"
//#include "Segment.hh"
#include "HashKey.hh"

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost {
namespace polygon {

using IntType = typename polytope::HashKey<2>::IntType;
using IntPoint = polytope::Point2<IntType>;

template <>
struct geometry_concept<IntPoint> { typedef point_concept type; };
  
template <>
struct point_traits<IntPoint> {
  typedef IntPoint point_type;
  typedef IntType coordinate_type;

  static inline coordinate_type get(const IntPoint& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};

template <>
struct point_mutable_traits<IntPoint> {
  typedef IntPoint point_type;
  typedef IntType coordinate_type;

  static inline void set(point_type& point, orientation_2d orient, coordinate_type value) {
    if (orient == HORIZONTAL)
      point.x = value;
    else
      point.y = value;
  }
  static inline point_type construct(coordinate_type x, coordinate_type y) {
    return point_type(x,y);
  }
};


// template <>
// struct geometry_concept<IntSegment> { typedef segment_concept type; };

// template <>
// struct point_traits<IntSegment> {
//   typedef IntType coordinate_type;
//   typedef IntPoint point_type;

//   static inline point_type get(const IntSegment& segment, direction_1d dir) {
//     return dir.to_int() ? segment.b : segment.a;
//   }
// };

inline IntPoint BoostToPolytope(const point_data<IntType>& point, const int index = 0) {
  return IntPoint(point.x(), point.y(), index);
}

inline std::vector<IntPoint> BoostToPolytope(const polygon_data<IntType>& polygon) {
  const auto N = polygon.size();
  std::vector<IntPoint> points;
  points.reserve(N);
  unsigned i = 0;
  for (const auto& p : polygon) {
    points.push_back(BoostToPolytope(p, i++));
  }
  if (N > 1 && points.front() == points.back()) {
    points.pop_back();
  }
  return points;
}

inline std::vector<IntPoint> BoostToPolytope(const polygon_with_holes_data<IntType>& polygon) {
  const auto N = polygon.size();
  std::vector<IntPoint> points;
  points.reserve(N);
  unsigned i = 0;
  for (const auto& p : polygon) {
    points.push_back(BoostToPolytope(p, i++));
  }
  if (N > 1 && points.front() == points.back()) {
    points.pop_back();
  }
  return points;
}

inline std::vector<IntPoint> outerPoints(const polygon_with_holes_data<IntType>& polygon) {
  std::vector<IntPoint> out;
  for(auto it = begin_points(polygon); it != end_points(polygon); ++it) {
    out.push_back(BoostToPolytope(*it));
  }
  return out;
}

inline std::vector<std::vector<IntPoint>> innerPoints(const polygon_with_holes_data<IntType>& polygon) {
  std::vector<std::vector<IntPoint>> out;
  auto hole_it = begin_holes(polygon);
  auto hole_end = end_holes(polygon);
  for(; hole_it != hole_end; ++hole_it) {
    const auto& hole = *hole_it;
    out.push_back(std::vector<IntPoint>());
    for (auto& point : hole) {
      out.back().push_back(BoostToPolytope(point));
    }
  }
  return out;
}

} //end boost namespace
} //end polygon namespace

#endif
