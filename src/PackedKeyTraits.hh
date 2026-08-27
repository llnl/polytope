//-----------------------------------------------------------------------------//
// PackedKeyTraits
//
// Fixed-width consecutive-coordinate key encoding.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_PackedKeyTraits__
#define __Polytope_PackedKeyTraits__

#include "QuantizedKeyTraits.hh"

namespace polytope {

template<int Dimension> struct PackedKeyTraits;

template<int Dimension,
         typename KeyT,
         typename UnsignedKeyT,
         typename CoordinateT,
         typename KeyHasherT,
         int KeyBits,
         int CoordinateBits>
using PackedKeyTraitsBase = QuantizedKeyTraitsBase<Dimension, KeyT, UnsignedKeyT,
                                                   CoordinateT, KeyHasherT, KeyBits,
                                                   CoordinateBits>;

template<int Dimension>
struct PackedKeyTraits : QuantizedKeyTraits<Dimension> {
  using Base = QuantizedKeyTraits<Dimension>;
  using Key = typename Base::Key;
  using UnsignedKey = typename Base::UnsignedKey;
  using PointType = typename Base::PointType;

  // x occupies the least-significant field, followed by y and z.
  static Key encode(const PointType& point) {
    UnsignedKey key = 0;
    constexpr UnsignedKey coordinateMask =
      (static_cast<UnsignedKey>(1) << Base::bitsPerCoordinate) - 1;

    for (int axis = 0; axis < Dimension; ++axis) {
      key |= (static_cast<UnsignedKey>(point[axis]) & coordinateMask)
             << (axis * Base::bitsPerCoordinate);
    }
    return static_cast<Key>(key);
  }

  static PointType decode(const Key key) {
    PointType point;
    constexpr UnsignedKey coordinateMask =
      (static_cast<UnsignedKey>(1) << Base::bitsPerCoordinate) - 1;
    const auto unsignedKey = static_cast<UnsignedKey>(key);

    for (int axis = 0; axis < Dimension; ++axis) {
      point[axis] = static_cast<typename Base::Coordinate>(
        (unsignedKey >> (axis * Base::bitsPerCoordinate)) & coordinateMask);
    }
    return point;
  }
};

template<int Dimension>
using PackedKey = QuantizedKey<Dimension>;

template<int Dimension>
using PackedKeyHasher = QuantizedKeyHasher<Dimension>;

} // namespace polytope

#endif
