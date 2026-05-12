//-----------------------------------------------------------------------------//
// KeyTraits
//
// Generalized class for handling hash keys.
//----------------------------------------------------------------------------//
#ifndef __Polytope_HashKey__
#define __Polytope_HashKey__

#include <cstdint>
#include <cstddef>
#include <memory>
#include <functional>

#include "polytope.hh"
#include "polytope_serialize.hh"

namespace polytope {

class HashKey2D {
public:
  HashKey2D() noexcept = default;
  explicit HashKey2D(uint64_t key) noexcept : m_key(key) {}

  uint64_t value() const noexcept { return m_key; }

  bool operator==(const HashKey2D& other) const noexcept {
    return m_key == other.m_key;
  }

  bool operator!=(const HashKey2D& other) const noexcept {
    return !(*this == other);
  }

  static unsigned bitsPerDim() { return 32; }

  // Interleave using Morton indexing
  template<typename CoordType> inline
  void interleave(const Point2<CoordType>& point) {
    m_key = 0;
    auto x = static_cast<uint32_t>(point.x);
    auto y = static_cast<uint32_t>(point.y);
    for (auto i = 0; i < bitsPerDim(); ++i) {
      m_key |= ((x >> i) & 1UL) << (2 * i);
      m_key |= ((y >> i) & 1UL) << (2 * i + 1);
    }
  }

  // Deinterleave Morton's indexing
  template<typename CoordType> inline
  Point2<CoordType> deinterleave() {
    uint32_t x = 0;
    uint32_t y = 0;
    for (unsigned i = 0; i < bitsPerDim(); ++i) {
      x |= ((m_key >> (2*i))   & 1UL) << i;
      y |= ((m_key >> (2*i+1)) & 1UL) << i;
    }
    return Point2<CoordType>(static_cast<CoordType>(x),
                             static_cast<CoordType>(y));
  }

  uint64_t m_key = 0;
};

class HashKey3D {
public:
  HashKey3D() noexcept = default;
  HashKey3D(uint64_t lo, uint64_t hi) noexcept : m_lo(lo), m_hi(hi) {}

  uint64_t lo() const noexcept { return m_lo; }

  uint64_t hi() const noexcept { return m_hi; }

  bool operator==(const HashKey3D& other) const noexcept {
    return m_lo == other.m_lo && m_hi == other.m_hi;
  }

  bool operator!=(const HashKey3D& other) const noexcept {
    return !(*this == other);
  }

  static unsigned bitsPerDim() { return 42; }

  // Interleave using Morton indexing
  template<typename CoordType> inline
  void interleave(const Point3<CoordType>& point) {
    auto x = static_cast<uint64_t>(point.x);
    auto y = static_cast<uint64_t>(point.y);
    auto z = static_cast<uint64_t>(point.z);
    m_lo = 0;
    m_hi = 0;
    for (auto i = 0; i < bitsPerDim(); ++i) {
      setBit128(3*i,   (x >> i) & 1ULL);
      setBit128(3*i+1, (y >> i) & 1ULL);
      setBit128(3*i+2, (z >> i) & 1ULL);
    }
  }

  // Deinterleave Morton's indexing
  template<typename CoordType> inline
  Point3<CoordType> deinterleave() {
    uint64_t x = 0, y = 0, z = 0;
    for (unsigned i = 0; i < bitsPerDim(); ++i) {
      x |= getBit128(3*i)   << i;
      y |= getBit128(3*i+1) << i;
      z |= getBit128(3*i+2) << i;
    }
    return Point3<CoordType>(static_cast<uint64_t>(x),
                             static_cast<uint64_t>(y),
                             static_cast<uint64_t>(z));
  }

  // Bit operations
  inline void setBit128(unsigned bitIndex, uint64_t bit) noexcept {
    if (bitIndex < 64) {
      m_lo |= (bit & 1ULL) << bitIndex;
    } else {
      m_hi |= (bit & 1ULL) << (bitIndex - 64);
    }
  }

  inline uint64_t getBit128(unsigned bitIndex) const noexcept {
    if (bitIndex < 64) {
      return (m_lo >> bitIndex) & 1ULL;
    }
    return (m_hi >> (bitIndex - 64)) & 1ULL;
  }

  uint64_t m_lo = 0;
  uint64_t m_hi = 0;
};

// Serialization
template<>
struct Serializer<HashKey2D> {
  static void serializeImpl(const HashKey2D& value,
                            std::vector<char>& buffer) {
    serialize(value.m_key, buffer);
  }

  static void deserializeImpl(HashKey2D& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    deserialize(value.m_key, bufItr, endItr);
  }
};

template<>
struct Serializer<HashKey3D> {
  static void serializeImpl(const HashKey3D& value,
                            std::vector<char>& buffer) {
    serialize(value.m_lo, buffer);
    serialize(value.m_hi, buffer);
  }

  static void deserializeImpl(HashKey3D& value,
                              std::vector<char>::const_iterator& bufItr,
                              const std::vector<char>::const_iterator& endItr) {
    deserialize(value.m_lo, bufItr, endItr);
    deserialize(value.m_hi, bufItr, endItr);
  }
};
}
#endif
