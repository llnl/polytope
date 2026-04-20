//-----------------------------------------------------------------------------//
// KeyTraits
//
// Encapsulate how we think about keys for the space filling curves.
// This is specialized assuming we're working with 64 bit uint64_t.
//----------------------------------------------------------------------------//
#ifndef __Polytope_Utilities_KeyTraits__
#define __Polytope_Utilities_KeyTraits__

#include <stdint.h>

namespace polytope {

struct KeyTraits {
  using Key = int;
  static const uint32_t numbits;
  static const uint32_t numbits1d;
  static const Key zero;
  static const Key one;
  static const Key two;
  static const Key minKey1d;
  static const Key maxKey1d;
  static const Key maxKey;
};

}

#endif
