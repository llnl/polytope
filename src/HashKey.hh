//-----------------------------------------------------------------------------//
// HashKey
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

template<int Dimension, typename RealType> struct HashKey;

template<typename RealType> struct HashKey<2, RealType> {
  using CoordHash = uint64_t;
  using CoordPoint = Point2<CoordHash>;
  using RealPoint = Point2<RealType>;

  static constexpr unsigned  flagBit()   { return 63; }
  static constexpr CoordHash FlagMask() {
    return static_cast<CoordHash>(1) << flagBit();
  }

  // Check if hash corresponds to the inner or outer box
  static bool getOuterFlag(const CoordHash& hash) {
    return (hash & FlagMask()) != 0;
  }

  // Set hash to correspond to outer box
  static void enableOuterFlag(CoordHash& hash) {
    hash |= FlagMask();
  }

  // Set hash to correspond to inner box
  static void disableOuterFlag(CoordHash& hash) {
    hash &= ~FlagMask();
  }

  static constexpr unsigned  num1DBits() { return 31; }
  static constexpr CoordHash coordMax()  { return (1ULL << (num1DBits() - 1ULL)) - 1ULL; }

  static CoordHash hash(const CoordPoint& point) {
    CoordHash key = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      key |= ((point.x >> i) & 1UL) << (2*i);
      key |= ((point.y >> i) & 1UL) << (2*i + 1);
    }
    return key;
  }

  static CoordPoint unhash(const CoordHash& key) {
    CoordPoint point(0, 0);
    for (auto i = 0; i < num1DBits(); ++i) {
      point.x |= ((key >> (2*i))   & 1UL) << i;
      point.y |= ((key >> (2*i+1)) & 1UL) << i;
    }
    return point;
  }

  static CoordHash hashPosition(const RealPoint& pos,
                                const RealPoint& length) {
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return hash(CoordPoint(pos, dx));
  }

  static CoordHash hashPosition(const RealPoint& pos,
                                const RealPoint& bhi,
                                const RealPoint& blo) {
    RealPoint length = bhi - blo;
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return hash(CoordPoint(pos, blo, dx));
  }

  static RealPoint unhashPosition(const CoordHash& hashed_pos,
                                  const RealPoint& length) {
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return unhash(hashed_pos).realPoint(dx);
  }

  static CoordHash unhashPosition(const CoordHash& hashed_pos,
                                  const RealPoint& bhi,
                                  const RealPoint& blo) {
    RealPoint length = bhi - blo;
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return unhash(hashed_pos).realPoint(blo, dx);
  }

};

template<typename RealType> struct HashKey<3, RealType> {
  using CoordHash = unsigned __int128;
  using CoordPoint = Point3<uint64_t>;
  using RealPoint = Point3<RealType>;

  static constexpr unsigned  flagBit()   { return 127; }
  static constexpr CoordHash FlagMask() {
    return static_cast<CoordHash>(1) << flagBit();
  }

  static bool getOuterFlag(const CoordHash& hash) {
    return (hash & FlagMask()) != 0;
  }

  static void enableOuterFlag(CoordHash& hash) {
    hash |= FlagMask();
  }

  static void disableOuterFlag(CoordHash& hash) {
    hash &= ~FlagMask();
  }

  static constexpr unsigned  num1DBits() { return 42; }
  static constexpr CoordHash coordMax()  { return (1ULL << (num1DBits() - 1ULL)) - 1ULL; }

  static CoordHash hash(const CoordPoint& point) {
    CoordHash key = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      key |= (static_cast<CoordHash>((point.x >> i) & 1ULL) << (3*i));
      key |= (static_cast<CoordHash>((point.y >> i) & 1ULL) << (3*i + 1));
      key |= (static_cast<CoordHash>((point.z >> i) & 1ULL) << (3*i + 2));
    }
    return key;
  }

  static CoordPoint unhash(const CoordHash& key) {
    uint64_t x = 0, y = 0, z = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      x |= ((key >> (3*i))   & 1ULL) << i;
      y |= ((key >> (3*i+1)) & 1ULL) << i;
      z |= ((key >> (3*i+2)) & 1ULL) << i;
    }
    return CoordPoint(x, y, z);
  }

  static CoordHash hashPosition(const RealPoint& pos,
                                const RealPoint& length) {
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return hash(CoordPoint(pos, dx));
  }

  static CoordHash hashPosition(const RealPoint& pos,
                                const RealPoint& bhi,
                                const RealPoint& blo) {
    RealPoint length = bhi - blo;
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return hash(CoordPoint(pos, blo, dx));
  }

  static RealPoint unhashPosition(const CoordHash& hashed_pos,
                                  const RealPoint& length) {
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return unhash(hashed_pos).realPoint(dx);
  }

  static RealPoint unhashPosition(const CoordHash& hashed_pos,
                                  const RealPoint& bhi,
                                  const RealPoint& blo) {
    RealPoint length = bhi - blo;
    RealPoint dx = length / static_cast<RealType>(coordMax());
    return unhash(hashed_pos).realPoint(blo, dx);
  }

};
}
#endif
