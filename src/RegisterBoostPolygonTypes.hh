//------------------------------------------------------------------------------
// This file is used to register our types with the Boost.Polygon library.
//------------------------------------------------------------------------------
#ifndef __polytope_RegisterBoostPolygonType__
#define __polytope_RegisterBoostPolygonType__

#include <cstdint>

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/voronoi.hpp>

// Include Boost multiprecision types
#include <boost/multiprecision/cpp_int.hpp>

#include "Point.hh"
#include "HashKey.hh"

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost {
namespace polygon {

#ifdef POLYTOPE_ENABLE_HIBIT2D
namespace detail {
template <>
struct voronoi_ctype_traits<std::int64_t> {
  typedef std::int64_t int_type;

  // For 32-bit inputs, int_x2_type is a 64-bit int.
  // For 64-bit inputs, we must use the library's multiprecision type scaled up.
  typedef extended_int<128> int_x2_type;
  typedef extended_int<128> uint_x2_type; // Maps signed logic to extended block

  // The core multiprecision int needs to scale to 512 bits (or 256 depending on depth,
  // but 512 safely covers maximum predicate multiplication chain for 64-bit inputs)
  typedef extended_int<512> big_int_type;

  // Output vertex type (64-bit double standard IEEE-754 floating point)
  typedef double fpt_type;

  // Extended exponent floating-point wrapper to handle high dynamic range calculations
  typedef extended_exponent_fpt<fpt_type> efpt_type;
  typedef ulp_comparison<fpt_type> ulp_cmp_type;

  // Required type converters for the internal predicates
  typedef type_converter_fpt to_fpt_converter_type;
  typedef type_converter_efpt to_efpt_converter_type;
};
} // namespace detail
#endif

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
