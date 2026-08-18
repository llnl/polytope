//------------------------------------------------------------------------------
// Wrapper-only helpers for the renovated Polytope Python bindings.
//------------------------------------------------------------------------------
#ifndef __Polytope_PYB11_helpers__
#define __Polytope_PYB11_helpers__

#include "HashKey.hh"
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
coordHashToPy(const typename HashKey<Dimension>::CoordHash& value) {
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
typename HashKey<Dimension>::CoordHash
pyToCoordHash(const py::object& value) {
#ifdef POLYTOPE_ENABLE_HIBIT2D
  return pyToInt128(value);
#else
  if constexpr (Dimension == 3) {
    return pyToInt128(value);
  } else {
    return value.cast<typename HashKey<Dimension>::CoordHash>();
  }
#endif
}

#ifdef POLYTOPE_ENABLE_MPI
template<int Dimension>
class PyDistributedTessellator: public DistributedTessellator<Dimension> {
public:
  using Base = Tessellator<Dimension, double>;

  explicit PyDistributedTessellator(py::object serialTessellator):
    DistributedTessellator<Dimension>(serialTessellator.cast<Base&>()),
    m_serialTessellator(std::move(serialTessellator)) {
  }

private:
  py::object m_serialTessellator;
};
#endif

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

inline
bool
isPythonSequence(const py::handle& value) {
  return PySequence_Check(value.ptr()) and
         not py::isinstance<py::str>(value) and
         not py::isinstance<py::bytes>(value);
}

template<typename ValueType>
std::vector<ValueType>
copyPyToVector(const py::handle& values,
               const std::string& name) {
  if (not isPythonSequence(values)) {
    throw py::type_error(name + " must be a nested sequence");
  }

  std::vector<ValueType> result;
  const auto seq = py::reinterpret_borrow<py::sequence>(values);
  result.reserve(seq.size());
  for (const auto value: seq) result.push_back(value.cast<ValueType>());
  return result;
}

inline
std::vector<std::vector<unsigned>>
copyFacetList(const py::object& facets,
              const std::string& name) {
  if (not isPythonSequence(facets)) {
    throw py::type_error(name + " must be a nested sequence");
  }

  std::vector<std::vector<unsigned>> result;
  const auto seq = facets.cast<py::sequence>();
  result.reserve(seq.size());
  for (const auto facet: seq) {
    result.push_back(copyPyToVector<unsigned>(facet, name));
  }
  return result;
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
