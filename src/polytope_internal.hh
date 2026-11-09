// polytope_internal.hh
//
// Put common includes for polytope here that you don't necessarily 
// want exposed in the public interface.
#ifndef __Polytope_Internal__
#define __Polytope_Internal__

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include "polytope.hh"

// An POLY_ASSERT macro, if one isn't already defined.

namespace polytope {

// Forward declare our helper abort method.
[[noreturn]] void internal_abort();

#ifdef POLYTOPE_ENABLE_DEBUG
#define POLY_ASSERT(x) \
  do { \
    if (!(x)) { \
      std::cerr << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
      internal_abort(); \
    } \
  } while (false)
#define POLY_ASSERT2(x, msg) \
  do { \
    if (!(x)) { \
      std::cerr << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
      internal_abort(); \
    } \
  } while (false)
#define POLY_BEGIN_CONTRACT_SCOPE { 
#define POLY_END_CONTRACT_SCOPE }
#else
#define POLY_ASSERT(x) do { } while (false)
#define POLY_ASSERT2(x, msg) do { } while (false)
#define POLY_BEGIN_CONTRACT_SCOPE if (false) {
#define POLY_END_CONTRACT_SCOPE }
#endif

#define POLY_CONTRACT_VAR(x) \
  do { \
    if (false) (void)&(x); \
  } while (false)

// A requirement contract that is always on to check user input.
#define POLY_VERIFY(x) \
  do { \
    if (!(x)) { \
      std::cerr << "Requirement " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
      internal_abort(); \
    } \
  } while (false)
#define POLY_VERIFY2(x, msg) \
  do { \
    if (!(x)) { \
      std::cerr << "Requirement " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
      internal_abort(); \
    } \
  } while (false)

#define POLY_CHECK(x) \
  do { \
    if (!(x)) { \
      std::cerr << "Check " << #x << " failed at " << __FILE__ << ":" << __LINE__ << std::endl; \
      internal_abort(); \
    } \
  } while (false)

#define POLY_CHECK2(x, msg) \
  do { \
    if (!(x)) { \
      std::cerr << "Check " << #x << " failed at " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
      internal_abort(); \
    } \
  } while (false)

} //end namespace polytope

#endif
