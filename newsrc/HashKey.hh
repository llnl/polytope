//-----------------------------------------------------------------------------//
// HashKey
//
// Generalized class for handling hash keys.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_HashKey__
#define __Polytope_HashKey__

#include <cstdint>
#include <cstddef>
#include <memory>
#include <functional>

#include "Point.hh"

namespace polytope {

struct Int128Hash {
  std::size_t operator()(__int128 v) const noexcept {
    std::uint64_t lo = static_cast<std::uint64_t>(v);
    std::uint64_t hi = static_cast<std::uint64_t>(v >> 64);

    std::size_t seed = static_cast<std::size_t>(lo);
    seed ^= static_cast<std::size_t>(hi) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

template<int Dimension> struct HashKey;

template<> struct HashKey<2> {
  using CoordHash = int64_t;
  using IntType = int; // Number of bits must exceed num1DBits
  using HashType = std::hash<CoordHash>;
  using IntPoint = Point<2, IntType>;

  static constexpr int       numBits()   { return 64; }
  static constexpr int       num1DBits() { return 31; }
  static constexpr IntType   coordMax()  { return (1ULL << (num1DBits() - 1ULL)) - 1ULL; }
  static constexpr CoordHash hashMax()   { return ((unsigned CoordHash)1 << (numBits() - 1ULL)) - 1ULL; }

  static CoordHash hash(const IntPoint& point) {
    CoordHash key = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      key |= (static_cast<CoordHash>((point.x >> i) & 1UL) << (2*i));
      key |= (static_cast<CoordHash>((point.y >> i) & 1UL) << (2*i + 1));
    }
    return key;
  }

  static IntPoint unhash(const CoordHash& key) {
    IntType x = 0, y = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      x |= ((key >> (2*i))   & 1UL) << i;
      y |= ((key >> (2*i+1)) & 1UL) << i;
    }
    return IntPoint(x, y);
  }
};

template<> struct HashKey<3> {
  using CoordHash = __int128;
  using HashType = Int128Hash;
  using IntType = int64_t; // Number of bits must exceed num1DBits
  using IntPoint = Point<3, IntType>;

  static constexpr int       numBits()   { return 128; }
  static constexpr int       num1DBits() { return 42; }
  static constexpr IntType   coordMax()  { return (1ULL << (num1DBits() - 1ULL)) - 1ULL; }
  static constexpr CoordHash hashMax()   { return ((unsigned CoordHash)1 << (numBits() - 1ULL)) - 1ULL; }

  static CoordHash hash(const IntPoint& point) {
    CoordHash key = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      key |= (static_cast<CoordHash>((point.x >> i) & 1ULL) << (3*i));
      key |= (static_cast<CoordHash>((point.y >> i) & 1ULL) << (3*i + 1));
      key |= (static_cast<CoordHash>((point.z >> i) & 1ULL) << (3*i + 2));
    }
    return key;
  }

  static IntPoint unhash(const CoordHash& key) {
    IntType x = 0, y = 0, z = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      x |= ((key >> (3*i))   & 1ULL) << i;
      y |= ((key >> (3*i+1)) & 1ULL) << i;
      z |= ((key >> (3*i+2)) & 1ULL) << i;
    }
    return IntPoint(x, y, z);
  }
};

}
#endif
