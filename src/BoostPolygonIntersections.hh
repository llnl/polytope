//-----------------------------------------------------------------------------//
// BoostPolygonIntersections
//
// Boost.Polygon-backed operations on quantized 2D polygons.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_BoostPolygonIntersections__
#define __Polytope_BoostPolygonIntersections__

#include <vector>

#include "RegisterBoostPolygonTypes.hh"

namespace polytope {

namespace bp = boost::polygon;
using namespace boost::polygon::operators;

using Polygon = bp::polygon_data<QuantizedCoordinate<2>>;
using PolygonWithHoles = bp::polygon_with_holes_data<QuantizedCoordinate<2>>;
using PolygonSet = bp::polygon_set_data<QuantizedCoordinate<2>>;

inline std::vector<PolygonWithHoles>
boostUnion(const PolygonWithHoles& first,
           const PolygonWithHoles& second) {
  std::vector<PolygonWithHoles> result;
  bp::assign(result, first | second);
  return result;
}

inline std::vector<PolygonWithHoles>
boostIntersect(const PolygonWithHoles& first,
               const PolygonWithHoles& second) {
  std::vector<PolygonWithHoles> result;
  bp::assign(result, first & second);
  return result;
}

// Clip a Polytope cell against a Boost polygon.
inline std::vector<PolygonWithHoles>
boostIntersect(const Cell<2, QuantizedCoordinate<2>>::CellType& cell,
               const PolygonWithHoles& boundary) {
  const auto boostCell = bp::polytopeToBoost(cell);
  std::vector<PolygonWithHoles> result;
  bp::assign(result, boostCell & boundary);
  return result;
}

// Merge a polygon into outPolygon when their union is a single polygon.
inline bool
validUnion(const PolygonWithHoles& polygon,
           PolygonWithHoles& outPolygon) {
  const auto unionResult = boostUnion(outPolygon, polygon);
  if (unionResult.size() != 1) return false;
  outPolygon = unionResult.front();
  return true;
}

} // namespace polytope

#endif
