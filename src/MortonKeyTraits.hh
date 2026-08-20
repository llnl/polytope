//-----------------------------------------------------------------------------//
// MortonKeyTraits
//
// Generalized class for handling Morton keys.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_MortonKeyTraits__
#define __Polytope_MortonKeyTraits__

#include <cstdint>
#include <cstddef>
#include <functional>

#include "Point.hh"

namespace polytope {

struct KeyHasher128 {
  std::size_t operator()(const __int128 key) const noexcept {
    constexpr std::uint64_t goldenRatio = 0x9e3779b97f4a7c15ULL;

    const auto lower = static_cast<std::uint64_t>(key);
    const auto upper = static_cast<std::uint64_t>(key >> 64);

    auto result = static_cast<std::size_t>(lower);
    result ^= static_cast<std::size_t>(upper) + goldenRatio +
              (result << 6) + (result >> 2);
    return result;
  }
};

template<int Dimension> struct MortonKeyTraits;

// Reversible Morton (Z-order) encoding for non-negative 2D coordinates.
// Key bits are interleaved in x0, y0, x1, y1, ... order.
template<> struct MortonKeyTraits<2> {
#ifdef POLYTOPE_ENABLE_HIBIT2D
  using Key = __int128;
  using Coordinate = std::int64_t;
  using KeyHasher = KeyHasher128;
  static constexpr int keyBitWidth = 128;
  // Exceeding 52 bits prevents consistent conversion to double.
  static constexpr int bitsPerCoordinate = 52;
#else
  using Key = std::int64_t;
  using Coordinate = int;
  using KeyHasher = std::hash<Key>;
  static constexpr int keyBitWidth = 64;
  static constexpr int bitsPerCoordinate = 31;
#endif

  using PointType = Point<2, Coordinate>;

  static constexpr int dimension = 2;
  static constexpr int wideBitWidth = 2 * bitsPerCoordinate;
  static constexpr int bigBitWidth = 3 * bitsPerCoordinate;

  static constexpr int numDim() { return dimension; }
  static constexpr int coordinateBits() { return bitsPerCoordinate; }

  static_assert(dimension * bitsPerCoordinate <= keyBitWidth,
                "The key must hold all interleaved coordinate bits.");

  // When multiplying two coordinates together
  static constexpr int wideWidth = wideBitWidth;
  using Wide = ac_int<wideWidth, true>;
  // When multiplying a Wide and a coordinate int
  static constexpr int bigWidth = bigBitWidth;
  using Big = ac_int<bigWidth, true>;

  static constexpr Coordinate maxCoordinate() {
    return (1ULL << (bitsPerCoordinate - 1)) - 1ULL;
  }
  static constexpr Key maxKey() {
    return ((unsigned Key)1 << (keyBitWidth - 1ULL)) - 1ULL;
  }

  static Key encode(const PointType& point) {
    Key key = 0;
    for (int bit = 0; bit < bitsPerCoordinate; ++bit) {
      key |= static_cast<Key>((point.x >> bit) & 1ULL) << (dimension * bit);
      key |= static_cast<Key>((point.y >> bit) & 1ULL) << (dimension * bit + 1);
    }
    return key;
  }

  static PointType decode(const Key key) {
    Coordinate x = 0;
    Coordinate y = 0;
    for (int bit = 0; bit < bitsPerCoordinate; ++bit) {
      x |= ((key >> (dimension * bit)) & 1ULL) << bit;
      y |= ((key >> (dimension * bit + 1)) & 1ULL) << bit;
    }
    return PointType(x, y);
  }
};

// Reversible Morton (Z-order) encoding for non-negative 3D coordinates.
// Key bits are interleaved in x0, y0, z0, x1, y1, z1, ... order.
template<> struct MortonKeyTraits<3> {
  using Key = __int128;
  using Coordinate = std::int64_t;
  using KeyHasher = KeyHasher128;
  using PointType = Point<3, Coordinate>;

  static constexpr int dimension = 3;
  static constexpr int keyBitWidth = 128;
  static constexpr int bitsPerCoordinate = 42;
  static constexpr int wideBitWidth = 2 * bitsPerCoordinate;
  static constexpr int bigBitWidth = 3 * bitsPerCoordinate;

  static constexpr int numDim() { return dimension; }
  static constexpr int coordinateBits() { return bitsPerCoordinate; }

  static_assert(dimension * bitsPerCoordinate <= keyBitWidth,
                "The key must hold all interleaved coordinate bits.");

  // When multiplying two coordinates together
  static constexpr int wideWidth = wideBitWidth;
  using Wide = ac_int<wideWidth, true>;
  // When multiplying a Wide and a coordinate int
  static constexpr int bigWidth = bigBitWidth;
  using Big = ac_int<bigWidth, true>;

  static constexpr Coordinate maxCoordinate() {
    return (1ULL << (bitsPerCoordinate - 1)) - 1ULL;
  }
  static constexpr Key maxKey() {
    return ((unsigned Key)1 << (keyBitWidth - 1ULL)) - 1ULL;
  }

  static Key encode(const PointType& point) {
    Key key = 0;
    for (int bit = 0; bit < bitsPerCoordinate; ++bit) {
      key |= static_cast<Key>((point.x >> bit) & 1ULL) << (dimension * bit);
      key |= static_cast<Key>((point.y >> bit) & 1ULL) << (dimension * bit + 1);
      key |= static_cast<Key>((point.z >> bit) & 1ULL) << (dimension * bit + 2);
    }
    return key;
  }

  static PointType decode(const Key key) {
    Coordinate x = 0;
    Coordinate y = 0;
    Coordinate z = 0;
    for (int bit = 0; bit < bitsPerCoordinate; ++bit) {
      x |= ((key >> (dimension * bit)) & 1ULL) << bit;
      y |= ((key >> (dimension * bit + 1)) & 1ULL) << bit;
      z |= ((key >> (dimension * bit + 2)) & 1ULL) << bit;
    }
    return PointType(x, y, z);
  }
};

template<int Dimension>
using MortonKey = typename MortonKeyTraits<Dimension>::Key;

template<int Dimension>
using QuantizedCoordinate = typename MortonKeyTraits<Dimension>::Coordinate;

template<int Dimension>
using QuantizedPoint = typename MortonKeyTraits<Dimension>::PointType;

template<int Dimension>
using MortonKeyHasher = typename MortonKeyTraits<Dimension>::KeyHasher;

}
#endif
