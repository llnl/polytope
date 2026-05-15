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

template<int Dimension> struct HashKey;

template<> struct HashKey<2> {
  using CoordHash = uint64_t;
  using IntType = unsigned int; // Number of bits must exceed num1DBits
  using IntPoint = typename PointType<2, IntType>::type;

  static constexpr unsigned flagBit()   { return 63; }
  static constexpr unsigned num1DBits() { return 31; }
  static constexpr IntType  coordMax()  { return (1ULL << (num1DBits() - 1ULL)) - 1ULL; }

  // Bit mask for the flag bit
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

  static CoordHash hash(const IntPoint& point) {
    CoordHash key = 0;
    for (auto i = 0; i < num1DBits(); ++i) {
      key |= ((point.x >> i) & 1UL) << (2*i);
      key |= ((point.y >> i) & 1UL) << (2*i + 1);
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
  using CoordHash = unsigned __int128;
  using IntType = uint64_t; // Number of bits must exceed num1DBits
  using IntPoint = typename PointType<3, IntType>::type;

  static constexpr unsigned flagBit()   { return 127; }
  static constexpr unsigned num1DBits() { return 42; }
  static constexpr IntType coordMax()  { return (1ULL << (num1DBits() - 1ULL)) - 1ULL; }

  // Bit mask for the flag bit
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

template<int Dimension, typename CoordType>
struct NewDimensionTraits {
  using PointType = typename PointType<Dimension, CoordType>::type;

  // Convert flattened array into vector of points
  static std::vector<PointType>
  extractCoords(const std::vector<CoordType>& allCoords,
                const std::vector<unsigned>& indices) {
    std::vector<PointType> result;
    result.reserve(indices.size());
    for (std::vector<unsigned>::const_iterator itr = indices.begin();
         itr != indices.end(); ++itr) {
      const unsigned i = *itr;
      result.push_back(PointType(&(allCoords[Dimension*i]), i));
    }
    return result;
  }

  static std::vector<PointType>
  extractCoords(const std::vector<CoordType>& allCoords) {
    std::vector<PointType> result;
    auto n = allCoords.size()/Dimension;
    result.reserve(n);
    for (auto i = 0; i < n; ++i) {
      result.push_back(PointType(&(allCoords[Dimension*i]), i));
    }
    return result;
  }
};

}
#endif
