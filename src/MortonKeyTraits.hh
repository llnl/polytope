//-----------------------------------------------------------------------------//
// MortonKeyTraits
//
// Morton (Z-order) key encoding for quantized coordinates.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_MortonKeyTraits__
#define __Polytope_MortonKeyTraits__

#include "QuantizedKeyTraits.hh"

namespace polytope {

template<int Dimension> struct MortonKeyTraits;

template<int Dimension,
         typename KeyT,
         typename UnsignedKeyT,
         typename CoordinateT,
         typename KeyHasherT,
         int KeyBits,
         int CoordinateBits>
using MortonKeyTraitsBase = QuantizedKeyTraitsBase<Dimension, KeyT, UnsignedKeyT,
                                                   CoordinateT, KeyHasherT, KeyBits,
                                                   CoordinateBits>;

template<int Dimension>
struct MortonKeyTraits : QuantizedKeyTraits<Dimension> {
  using Base = QuantizedKeyTraits<Dimension>;
  using Key = typename Base::Key;
  using PointType = typename Base::PointType;

  // Key bits are interleaved in x0, y0, z0, x1, y1, z1, ... order.
  static Key encode(const PointType& point) {
    Key key = 0;
    for (int bit = 0; bit < Base::bitsPerCoordinate; ++bit) {
      for (int axis = 0; axis < Dimension; ++axis) {
        const int keyBit = Dimension * bit + axis;
        key |= static_cast<Key>((point[axis] >> bit) & 1ULL) << keyBit;
      }
    }
    return key;
  }

  static PointType decode(const Key key) {
    PointType point;
    for (int bit = 0; bit < Base::bitsPerCoordinate; ++bit) {
      for (int axis = 0; axis < Dimension; ++axis) {
        const int keyBit = Dimension * bit + axis;
        point[axis] |= ((key >> keyBit) & 1ULL) << bit;
      }
    }
    return point;
  }
};

// Compatibility aliases for existing Morton-facing APIs.
template<int Dimension>
using MortonKey = QuantizedKey<Dimension>;

template<int Dimension>
using MortonKeyHasher = QuantizedKeyHasher<Dimension>;

} // namespace polytope

#endif
