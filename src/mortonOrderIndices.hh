//-----------------------------------------------------------------------------//
// mortonOrderIndicies
//
// Compute the Morton ordered hashed indicies for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Ported from the Spheral++ algorithm.
//----------------------------------------------------------------------------//
#ifndef __Polytope_mortonOrderIndices__
#define __Polytope_mortonOrderIndices__

#include <vector>

#include "polytope.hh"
#include "KeyTraits.hh"

namespace {

//------------------------------------------------------------------------------
// Stand alone functions to do the interleaved bit hashing of positions.
//------------------------------------------------------------------------------
// 2D
inline
polytope::KeyTraits::Key
hashPosition(const polytope::Point2<polytope::KeyTraits::Key>& xoff) {

  typedef polytope::KeyTraits::Key Key;

  // Get the bits for each dimension.
  POLY_ASSERT(xoff.x <= polytope::KeyTraits::maxKey1d);
  POLY_ASSERT(xoff.y <= polytope::KeyTraits::maxKey1d);

  // Interleave the bits.
  Key result = polytope::KeyTraits::zero;
  Key mask = polytope::KeyTraits::one;
  for (size_t i = 0; i != polytope::KeyTraits::numbits1d; ++i, mask <<= 1) {
    result += (((xoff.x & mask) << i) +
               ((xoff.y & mask) << (i + 1)));
  }

  // That's it.
  POLY_ASSERT(result <= polytope::KeyTraits::maxKey);
  return result;
}

// 3-D
inline
polytope::KeyTraits::Key
hashPosition(const polytope::Point3<polytope::KeyTraits::Key>& xoff) {

  typedef polytope::KeyTraits::Key Key;

  // Get the bits for each dimension.
  POLY_ASSERT(xoff.x <= polytope::KeyTraits::maxKey1d);
  POLY_ASSERT(xoff.y <= polytope::KeyTraits::maxKey1d);
  POLY_ASSERT(xoff.z <= polytope::KeyTraits::maxKey1d);

  // Interleave the bits.
  Key result = polytope::KeyTraits::zero;
  Key mask = polytope::KeyTraits::one;
  for (size_t i = 0; i != polytope::KeyTraits::numbits1d; ++i, mask <<= 1) {
    result += (((xoff.x & mask) << i) +
               ((xoff.y & mask) << (i + 1)) +
               ((xoff.z & mask) << (i + 2)));
  }

  // That's it.
  POLY_ASSERT(result <= polytope::KeyTraits::maxKey);
  return result;
}

}

//------------------------------------------------------------------------------
// Hash the positions into their tree ordered indicies.
//------------------------------------------------------------------------------
namespace polytope {

template<typename PointType>
std::vector<KeyTraits::Key>
mortonOrderIndices(const std::vector<PointType>& points) {

  typedef typename KeyTraits::Key Key;

  // Prepare the result.
  const unsigned npoints = points.size();
  std::vector<KeyTraits::Key> result;
  result.reserve(npoints);

  // Go over all points.
  for (typename std::vector<PointType>::const_iterator itr = points.begin();
       itr != points.end();
       ++itr) result.push_back(hashPosition(*itr));

  POLY_ASSERT(result.size() == npoints);
  return result;
}

}

#endif
