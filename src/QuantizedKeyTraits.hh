//-----------------------------------------------------------------------------//
// QuantizedKeyTraits
//
// Common representation for keys built from quantized coordinates.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_QuantizedKeyTraits__
#define __Polytope_QuantizedKeyTraits__

#include <cstddef>
#include <cstdint>
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

template<int Dimension> struct QuantizedKeyTraits;

template<int Dimension,
         typename KeyT,
         typename UnsignedKeyT,
         typename CoordinateT,
         typename KeyHasherT,
         int KeyBits,
         int CoordinateBits>
struct QuantizedKeyTraitsBase {
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
                "The key must hold all coordinate bits.");

  static constexpr int numDim() { return dimension; }
  static constexpr int coordinateBits() { return bitsPerCoordinate; }

  static constexpr int wideWidth = wideBitWidth;
  using Wide = ac_int<wideWidth, true>;

  static constexpr int bigWidth = bigBitWidth;
  using Big = ac_int<bigWidth, true>;

  static constexpr Coordinate maxCoordinate() {
    return (1ULL << (bitsPerCoordinate - 1)) - 1ULL;
  }

  // The maximum value representable by the signed key storage type.
  static constexpr Key maxKey() {
    return (static_cast<UnsignedKey>(1) << (keyBitWidth - 1ULL)) - 1ULL;
  }
};

#ifdef POLYTOPE_ENABLE_HIBIT2D
// 52 bits per dimension. Exceeding this prevents consistent conversion to double.
template<> struct QuantizedKeyTraits<2> :
  QuantizedKeyTraitsBase<2, __int128, unsigned __int128, std::int64_t,
                         KeyHasher128, 128, 52> {};
#else
// 31 bits per dimension.
template<> struct QuantizedKeyTraits<2> :
  QuantizedKeyTraitsBase<2, std::int64_t, std::uint64_t, int,
                         std::hash<std::int64_t>, 64, 31> {};
#endif

template<> struct QuantizedKeyTraits<3> :
  QuantizedKeyTraitsBase<3, __int128, unsigned __int128, std::int64_t,
                         KeyHasher128, 128, 42> {};

template<int Dimension>
using QuantizedKey = typename QuantizedKeyTraits<Dimension>::Key;

template<int Dimension>
using QuantizedCoordinate = typename QuantizedKeyTraits<Dimension>::Coordinate;

template<int Dimension>
using QuantizedPoint = typename QuantizedKeyTraits<Dimension>::PointType;

template<int Dimension>
using QuantizedKeyHasher = typename QuantizedKeyTraits<Dimension>::KeyHasher;

} // namespace polytope

#endif
