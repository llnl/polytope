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

template<int Dimension,
         typename KeyT,
         typename UnsignedKeyT,
         typename CoordinateT,
         typename KeyHasherT,
         int KeyBits,
         int CoordinateBits>
struct MortonKeyTraitsBase {
  using Key = KeyT;
  using UnsignedKey = UnsignedKeyT;
  using Coordinate = CoordinateT;
  using KeyHasher = KeyHasherT;
  using PointType = Point<Dimension, Coordinate>;

  static constexpr int dimension = Dimension;
  static constexpr int keyBitWidth = KeyBits;
  static constexpr int bitsPerCoordinate = CoordinateBits;
  static constexpr int wideBitWidth = 2 * bitsPerCoordinate;
  static constexpr int bigBitWidth = 3 * bitsPerCoordinate;

  static_assert(dimension * bitsPerCoordinate <= keyBitWidth,
                "The key must hold all interleaved coordinate bits.");

  static constexpr int numDim() { return dimension; }
  static constexpr int coordinateBits() { return bitsPerCoordinate; }

  // When multiplying two coordinates together.
  static constexpr int wideWidth = wideBitWidth;
  using Wide = ac_int<wideWidth, true>;
  // When multiplying a Wide and a coordinate integer.
  static constexpr int bigWidth = bigBitWidth;
  using Big = ac_int<bigWidth, true>;

  static constexpr Coordinate maxCoordinate() {
    return (1ULL << (bitsPerCoordinate - 1)) - 1ULL;
  }

  static constexpr Key maxKey() {
    return (static_cast<UnsignedKey>(1) << (keyBitWidth - 1ULL)) - 1ULL;
  }

  static Key encode(const PointType& point) {
    Key key = 0;
    for (int bit = 0; bit < bitsPerCoordinate; ++bit) {
      for (int axis = 0; axis < dimension; ++axis) {
        const int keyBit = dimension * bit + axis;
        key |= static_cast<Key>((point[axis] >> bit) & 1ULL) << keyBit;
      }
    }
    return key;
  }

  static PointType decode(const Key key) {
    PointType point;
    for (int bit = 0; bit < bitsPerCoordinate; ++bit) {
      for (int axis = 0; axis < dimension; ++axis) {
        const int keyBit = dimension * bit + axis;
        point[axis] |= ((key >> keyBit) & 1ULL) << bit;
      }
    }
    return point;
  }
};

// Reversible Morton (Z-order) encoding for non-negative 2D coordinates.
// Key bits are interleaved in x0, y0, x1, y1, ... order.
#ifdef POLYTOPE_ENABLE_HIBIT2D
// Exceeding 52 coordinate bits prevents consistent conversion to double.
template<> struct MortonKeyTraits<2> :
  MortonKeyTraitsBase<2, __int128, unsigned __int128, std::int64_t, KeyHasher128, 128, 52> {};
#else
template<> struct MortonKeyTraits<2> :
  MortonKeyTraitsBase<2, std::int64_t, std::uint64_t, int, std::hash<std::int64_t>, 64, 31> {};
#endif

// Reversible Morton (Z-order) encoding for non-negative 3D coordinates.
// Key bits are interleaved in x0, y0, z0, x1, y1, z1, ... order.
template<> struct MortonKeyTraits<3> :
  MortonKeyTraitsBase<3, __int128, unsigned __int128, std::int64_t, KeyHasher128, 128, 42> {};

// Specialized aliases to use throughout code
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
