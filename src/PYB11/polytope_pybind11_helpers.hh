//------------------------------------------------------------------------------
// Wrapper-only helpers for the renovated Polytope Python bindings.
//------------------------------------------------------------------------------
#ifndef __Polytope_PYB11_helpers__
#define __Polytope_PYB11_helpers__

#include "QuantizedKeyTraits.hh"
#include "Point.hh"
#include "polytope.hh"

#ifdef POLYTOPE_ENABLE_SILO
#include "SiloWriter.hh"
#endif

#ifdef POLYTOPE_ENABLE_BOOST
#include "BoostTessellator.hh"
#endif
#ifdef POLYTOPE_ENABLE_TRIANGLE
#include "TriangleTessellator.hh"
#endif
#ifdef POLYTOPE_ENABLE_MPI
#include "DistributedTessellator.hh"
#endif

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <array>

namespace py = pybind11;

namespace polytope {
namespace pybind11_helpers {

inline
std::string
unsignedInt128ToString(unsigned __int128 value) {
  if (value == 0) return "0";

  std::string result;
  while (value > 0) {
    const auto digit = static_cast<unsigned>(value % 10);
    result.push_back(static_cast<char>('0' + digit));
    value /= 10;
  }
  std::reverse(result.begin(), result.end());
  return result;
}

inline
py::object
int128ToPy(const __int128 value) {
  const bool negative = value < 0;
  const auto magnitude = negative ?
    static_cast<unsigned __int128>(-(value + 1)) + 1 :
    static_cast<unsigned __int128>(value);
  auto text = unsignedInt128ToString(magnitude);
  if (negative) text.insert(text.begin(), '-');

  char* end = nullptr;
  PyObject* result = PyLong_FromString(text.c_str(), &end, 10);
  if (result == nullptr or end == nullptr or *end != '\0') {
    throw py::error_already_set();
  }
  return py::reinterpret_steal<py::object>(result);
}

inline
py::object
uint128ToPy(const unsigned __int128 value) {
  const auto text = unsignedInt128ToString(value);
  char* end = nullptr;
  PyObject* result = PyLong_FromString(text.c_str(), &end, 10);
  if (result == nullptr or end == nullptr or *end != '\0') {
    throw py::error_already_set();
  }
  return py::reinterpret_steal<py::object>(result);
}

inline
unsigned __int128
parseUnsignedInt128(const std::string& text,
                    std::size_t pos,
                    const unsigned __int128 maxValue) {
  unsigned __int128 result = 0;
  bool foundDigit = false;
  for (; pos < text.size(); ++pos) {
    const auto ch = static_cast<unsigned char>(text[pos]);
    if (!std::isdigit(ch)) {
      throw py::value_error("Expected a Python integer");
    }
    foundDigit = true;
    const auto digit = static_cast<unsigned>(text[pos] - '0');
    if (result > (maxValue - digit)/10) {
      throw py::value_error("Python integer is outside the supported __int128 range");
    }
    result = 10*result + digit;
  }
  if (!foundDigit) {
    throw py::value_error("Expected a Python integer");
  }
  return result;
}

inline
__int128
pyToInt128(const py::object& value) {
  const auto pyint = py::module_::import("builtins").attr("int")(value);
  const std::string text = py::str(pyint);
  std::size_t pos = 0;
  bool negative = false;
  if (text[pos] == '-' or text[pos] == '+') {
    negative = text[pos] == '-';
    ++pos;
  }

  const auto maxPositive = (static_cast<unsigned __int128>(1) << 127) - 1;
  const auto maxMagnitude = negative ? (static_cast<unsigned __int128>(1) << 127) : maxPositive;
  const auto magnitude = parseUnsignedInt128(text, pos, maxMagnitude);
  if (magnitude == 0) return 0;
  return negative ?
    (-static_cast<__int128>(magnitude - 1) - 1) :
    static_cast<__int128>(magnitude);
}

inline
unsigned __int128
pyToUInt128(const py::object& value) {
  const auto pyint = py::module_::import("builtins").attr("int")(value);
  const std::string text = py::str(pyint);
  std::size_t pos = 0;
  if (text[pos] == '+') {
    ++pos;
  } else if (text[pos] == '-') {
    throw py::value_error("Negative Python integer cannot be converted to unsigned __int128");
  }
  return parseUnsignedInt128(text, pos, ~static_cast<unsigned __int128>(0));
}

template<int Dimension>
py::object
keyToPy(const QuantizedKey<Dimension>& value) {
#ifdef POLYTOPE_ENABLE_HIBIT2D
  return int128ToPy(value);
#else
  if constexpr (Dimension == 3) {
    return int128ToPy(value);
  } else {
    return py::int_(value);
  }
#endif
}

template<int Dimension>
QuantizedKey<Dimension>
pyToKey(const py::object& value) {
#ifdef POLYTOPE_ENABLE_HIBIT2D
  return pyToInt128(value);
#else
  if constexpr (Dimension == 3) {
    return pyToInt128(value);
  } else {
    return value.cast<QuantizedKey<Dimension>>();
  }
#endif
}

template<int Dimension, typename CoordType>
Point<Dimension, CoordType>
pyToPoint(const py::object& value,
          const unsigned index = 0) {
  using PointType = Point<Dimension, CoordType>;
  if (py::isinstance<PointType>(value)) {
    auto result = value.cast<PointType>();
    return result;
  }
  Point<Dimension, CoordType> out;
  const py::sequence seq = py::reinterpret_borrow<py::sequence>(value);
  const py::ssize_t size = py::len(seq);
  if (size != Dimension) {
    throw py::value_error("Point is wrong dimension");
  }
  for (py::ssize_t d = 0; d < size; ++d) {
    out[d] = py::cast<CoordType>(seq[d]);
  }
  out.index = index;
  return out;
}

// Helper routine to convert flattened or nested python sequences to
// vector of points
inline
bool
isPythonSequence(const py::handle& value) {
  return PySequence_Check(value.ptr()) and
         not py::isinstance<py::str>(value) and
         not py::isinstance<py::bytes>(value);
}

// Accept either [(x, y), ...] / [(x, y, z), ...] or a flat coordinate list.
template<int Dimension, typename CoordType>
std::vector<Point<Dimension, CoordType>>
copyCoords(const py::object& coords) {
  using PointType = Point<Dimension, CoordType>;
  using PointVector = std::vector<PointType>;
  if (py::isinstance<PointVector>(coords)) {
    auto result = coords.cast<PointVector>();
    return result;
  } else if (py::isinstance<PointType>(coords)) {
    return {pyToPoint<Dimension, CoordType>(coords, 0)};
  }
  PointVector result;
  POLY_CHECK2(isPythonSequence(coords), "Must pass a sequence to copyCoords");
  const py::sequence seq = py::reinterpret_borrow<py::sequence>(coords);
  const py::ssize_t size = py::len(seq);
  if (size == 0) return result;
  const auto first = seq[0];
  // If it is a nested list
  if (py::isinstance<PointType>(first) || isPythonSequence(first)) {
    result.reserve(size);
    for (py::ssize_t i = 0; i < size; ++i) {
      result.push_back(pyToPoint<Dimension, CoordType>(seq[i], i));
    }
    return result;
  }
  std::vector<CoordType> svec;
  svec.reserve(size);
  for (const py::handle value : seq) {
    svec.push_back(value.cast<CoordType>());
  }
  result = extractCoords<Dimension, CoordType>(svec);
  return result;
}

// Helper routines for converting triangle and neighbor lists
// Assumes input object is something like [[0, 1, 2], [2, 3, 4],...
inline
std::vector<std::array<int, 3>>
copyPyToTriList(const py::object& values,
                const std::string& func_name) {
  POLY_CHECK2(isPythonSequence(values), "Must pass a sequence to " << func_name);
  const py::sequence seq = py::reinterpret_borrow<py::sequence>(values);
  const py::ssize_t size = py::len(seq);
  std::vector<std::array<int, 3>> result;
  if (size == 0) return result;
  result.reserve(size);
  POLY_CHECK2(isPythonSequence(seq[0]), "Must nested sequence to " << func_name);
  for (py::ssize_t i = 0; i < size; ++i) {
    const py::sequence seqi = py::reinterpret_borrow<py::sequence>(seq[i]);
    POLY_ASSERT2(py::len(seqi) == 3, "Inner sequence in " << func_name << " must be length 3");
    result.push_back({seqi[0].cast<int>(), seqi[1].cast<int>(), seqi[2].cast<int>()});
  }
  return result;
}

// Helper routines for converting facet and hole lists
template<typename ValueType>
std::vector<ValueType>
copyPyToVector(const py::handle& values,
               const std::string& name) {
  if (not isPythonSequence(values)) {
    throw py::type_error(name + " must be a nested sequence");
  }

  std::vector<ValueType> result;
  const auto seq = py::reinterpret_borrow<py::sequence>(values);
  const py::ssize_t size = py::len(seq);
  result.reserve(size);
  for (const auto value: seq) result.push_back(value.cast<ValueType>());
  return result;
}

template<typename ValueType>
std::vector<std::vector<ValueType>>
copyNestedPyToVector(const py::object& values,
                     const std::string& name) {
  if (not isPythonSequence(values)) {
    throw py::type_error(name + " must be a nested sequence");
  }

  std::vector<std::vector<ValueType>> result;
  const auto seq = values.cast<py::sequence>();
  const py::ssize_t size = py::len(seq);
  result.reserve(size);
  for (const auto values: seq) {
    result.push_back(copyPyToVector<ValueType>(values, name));
  }
  return result;
}

inline
std::vector<std::vector<unsigned>>
copyFacetList(const py::object& facets,
              const std::string& name) {
  return copyNestedPyToVector<unsigned>(facets, name);
}

inline
std::vector<std::vector<std::vector<unsigned>>>
copyHoleList(const py::object& holes,
             const std::string& name) {
  if (not isPythonSequence(holes)) {
    throw py::type_error(name + " must be a nested sequence");
  }

  std::vector<std::vector<std::vector<unsigned>>> result;
  const auto seq = holes.cast<py::sequence>();
  result.reserve(seq.size());
  for (const auto hole: seq) {
    result.push_back(copyFacetList(py::reinterpret_borrow<py::object>(hole), name));
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

template<int Dimension, typename CoordType>
py::list
nestedPointsAsTuples(const std::vector<std::vector<Point<Dimension, CoordType>>>& pointvec) {
  py::list result;
  for (const auto& points : pointvec) {
    py::list r2;
    for (const auto& point: points) {
      py::tuple tup(Dimension);
      for (auto i = 0; i < Dimension; ++i) tup[i] = point[i];
      r2.append(tup);
    }
    result.append(r2);
  }
  return result;
}

} // namespace pybind11_helpers
} // namespace polytope

namespace pybind11 {
namespace detail {

template<>
struct type_caster<__int128> {
public:
  PYBIND11_TYPE_CASTER(__int128, _("int"));

  bool load(handle src, bool) {
    value = polytope::pybind11_helpers::pyToInt128(py::reinterpret_borrow<py::object>(src));
    return true;
  }

  static handle cast(const __int128 src, return_value_policy, handle) {
    return polytope::pybind11_helpers::int128ToPy(src).release();
  }
};

template<>
struct type_caster<unsigned __int128> {
public:
  PYBIND11_TYPE_CASTER(unsigned __int128, _("int"));

  bool load(handle src, bool) {
    value = polytope::pybind11_helpers::pyToUInt128(py::reinterpret_borrow<py::object>(src));
    return true;
  }

  static handle cast(const unsigned __int128 src, return_value_policy, handle) {
    return polytope::pybind11_helpers::uint128ToPy(src).release();
  }
};

} // namespace detail
} // namespace pybind11

#endif
