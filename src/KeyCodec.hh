//-----------------------------------------------------------------------------//
// KeyCodec
//
// Runtime selection between supported quantized-point key encodings.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_KeyCodec__
#define __Polytope_KeyCodec__

#include "MortonKeyTraits.hh"
#include "PackedKeyTraits.hh"
#include "QuantizedKeyTraits.hh"

namespace polytope {

enum class KeyEncoding {
  Morton,
  Packed
};

template<int Dimension>
class KeyCodec {
public:
  using Traits = QuantizedKeyTraits<Dimension>;
  using MortonTraits = MortonKeyTraits<Dimension>;
  using PackedTraits = PackedKeyTraits<Dimension>;
  using Key = typename Traits::Key;
  using Coordinate = typename Traits::Coordinate;
  using KeyHasher = typename Traits::KeyHasher;
  using PointType = typename Traits::PointType;

  explicit KeyCodec(const KeyEncoding encoding = KeyEncoding::Morton):
    m_encoding(encoding) {}

  KeyEncoding encoding() const {
    return m_encoding;
  }

  void setEncoding(const KeyEncoding encoding) {
    m_encoding = encoding;
  }

  // A display-ready name for the currently selected key encoding.
  const std::string keyName() const noexcept {
    switch (m_encoding) {
      case KeyEncoding::Morton:
        return "Morton";

      case KeyEncoding::Packed:
        return "Packed";
    }

    return "Unknown";
  }

  static constexpr int numDim() {
    return Traits::numDim();
  }

  static constexpr int coordinateBits() {
    return Traits::coordinateBits();
  }

  static constexpr Coordinate maxCoordinate() {
    return Traits::maxCoordinate();
  }

  Key encode(const PointType& point) const {
    switch (m_encoding) {
      case KeyEncoding::Morton:
        return MortonTraits::encode(point);

      case KeyEncoding::Packed:
        return PackedTraits::encode(point);
    }

    return MortonTraits::encode(point);
  }

  PointType decode(const Key key) const {
    switch (m_encoding) {
      case KeyEncoding::Morton:
        return MortonTraits::decode(key);

      case KeyEncoding::Packed:
        return PackedTraits::decode(key);
    }

    return MortonTraits::decode(key);
  }

private:
  KeyEncoding m_encoding;
};

} // namespace polytope

#endif
