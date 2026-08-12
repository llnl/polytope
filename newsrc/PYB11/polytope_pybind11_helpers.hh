//------------------------------------------------------------------------------
// Wrapper-only helpers for the renovated Polytope Python bindings.
//------------------------------------------------------------------------------
#ifndef __Polytope_new_PYB11_helpers__
#define __Polytope_new_PYB11_helpers__

#include "Point.hh"
#include "polytope.hh"

#ifdef POLYTOPE_ENABLE_BOOST
#include "BoostTessellator.hh"
#endif
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <vector>

namespace py = pybind11;

namespace polytope {
namespace pybind11_helpers {

// Accept either [(x, y), ...] / [(x, y, z), ...] or a flat coordinate list.
template<int Dimension, typename RealType>
std::vector<RealType>
copyCoords(const py::object& coords) {
  std::vector<RealType> result;
  if (py::isinstance<py::list>(coords) || py::isinstance<py::tuple>(coords)) {
    const auto seq = coords.cast<py::sequence>();
    if (seq.size() == 0) return result;

    const auto first = seq[0];
    if (py::isinstance<py::list>(first) || py::isinstance<py::tuple>(first)) {
      for (const auto item: seq) {
        const auto point = py::reinterpret_borrow<py::sequence>(item);
        if (point.size() != Dimension) {
          throw py::value_error("Coordinate tuple has the wrong dimension");
        }
        for (const auto value: point) result.push_back(value.cast<RealType>());
      }
    } else {
      for (const auto value: seq) result.push_back(value.cast<RealType>());
    }
  } else {
    result = coords.cast<std::vector<RealType>>();
  }
  if (result.size() % Dimension != 0) {
    throw py::value_error("Coordinate list length is not divisible by dimension");
  }
  return result;
}

template<int Dimension, typename CoordType>
py::list
pointsAsTuples(const std::vector<Point<Dimension, CoordType>>& points) {
  py::list result;
  for (const auto& point: points) {
    py::tuple tup(Dimension);
    for (auto i = 0; i < Dimension; ++i) tup[i] = point[i];
    result.append(tup);
  }
  return result;
}

} // namespace pybind11_helpers
} // namespace polytope

#endif
