//------------------------------------------------------------------------------
// This file is used to register our types with the Boost.Polygon library.
// It also provides routines for converting between Boost and Polytope types.
//------------------------------------------------------------------------------
#ifndef __polytope_RegisterBoostPolygonType__
#define __polytope_RegisterBoostPolygonType__

#include <cstdint>

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/voronoi.hpp>

// Include Boost multiprecision types
#include <boost/multiprecision/cpp_int.hpp>

#include "QuantizedKeyTraits.hh"
#include "Cell.hh"

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

using QuantizedCoordinate2D = polytope::QuantizedCoordinate<2>;
using QuantizedPoint2D = polytope::QuantizedPoint<2>;
using QuantizedCell = polytope::Cell<2, QuantizedCoordinate2D>;

template <>
struct geometry_concept<QuantizedPoint2D> { typedef point_concept type; };

template <>
struct point_traits<QuantizedPoint2D> {
  typedef QuantizedPoint2D point_type;
  typedef QuantizedCoordinate2D coordinate_type;

  static inline coordinate_type get(const QuantizedPoint2D& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};

template <>
struct point_mutable_traits<QuantizedPoint2D> {
  typedef QuantizedPoint2D point_type;
  typedef QuantizedCoordinate2D coordinate_type;

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

inline QuantizedPoint2D BoostToPolytope(const point_data<QuantizedCoordinate2D>& point, const int index = 0) {
  return QuantizedPoint2D(point.x(), point.y(), index);
}

inline std::vector<QuantizedPoint2D> BoostToPolytope(const polygon_data<QuantizedCoordinate2D>& polygon) {
  const auto N = polygon.size();
  std::vector<QuantizedPoint2D> points;
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

inline std::vector<QuantizedPoint2D> BoostToPolytope(const polygon_with_holes_data<QuantizedCoordinate2D>& polygon) {
  const auto N = polygon.size();
  std::vector<QuantizedPoint2D> points;
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

inline std::vector<QuantizedPoint2D> outerPoints(const polygon_with_holes_data<QuantizedCoordinate2D>& polygon) {
  std::vector<QuantizedPoint2D> out;
  for(auto it = begin_points(polygon); it != end_points(polygon); ++it) {
    out.push_back(BoostToPolytope(*it));
  }
  return out;
}

inline std::vector<std::vector<QuantizedPoint2D>> innerPoints(const polygon_with_holes_data<QuantizedCoordinate2D>& polygon) {
  std::vector<std::vector<QuantizedPoint2D>> out;
  auto hole_it = begin_holes(polygon);
  auto hole_end = end_holes(polygon);
  for(; hole_it != hole_end; ++hole_it) {
    const auto& hole = *hole_it;
    out.push_back(std::vector<QuantizedPoint2D>());
    for (auto& point : hole) {
      out.back().push_back(BoostToPolytope(point));
    }
  }
  return out;
}

inline polygon_with_holes_data<QuantizedCoordinate2D>
polytopeToBoost(const QuantizedCell& cell) {
  polygon_with_holes_data<QuantizedCoordinate2D> polygon;
  set_points(polygon, cell.points().begin(), cell.points().end());
  return polygon;
}

} //end boost namespace
} //end polygon namespace

#endif
