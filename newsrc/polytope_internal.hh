// polytope_internal.hh
//
// Put common includes for polytope here that you don't necessarily 
// want exposed in the public interface.
#ifndef POLYTOPE_INTERNAL_HH
#define POLYTOPE_INTERNAL_HH

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include "polytope.hh"

// An POLY_ASSERT macro, if one isn't already defined.

namespace polytope {

// Forward declare our helper abort method.
void internal_abort();

#ifdef POLYTOPE_ENABLE_DEBUG
#define POLY_ASSERT(x) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
    internal_abort(); \
  }
#define POLY_ASSERT2(x, msg) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    internal_abort(); \
  }
#define POLY_BEGIN_CONTRACT_SCOPE { 
#define POLY_END_CONTRACT_SCOPE }
#else
#define POLY_ASSERT(x)
#define POLY_ASSERT2(x, msg)
#define POLY_BEGIN_CONTRACT_SCOPE if (false) {
#define POLY_END_CONTRACT_SCOPE }
#endif

#define POLY_CONTRACT_VAR(x) if (0 && &x == &x){}

// A requirement contract that is always on to check user input.
#define POLY_VERIFY(x)                                     \
  if (!(x)) {                                              \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
    internal_abort();                                         \
  }
#define POLY_VERIFY2(x, msg)                                            \
  if (!(x)) {                                                           \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    internal_abort();                                                   \
  }

#define POLY_CHECK(x)                                                   \
  if (!(x)) {                                                           \
    std::cout << "Check " << #x << " failed at " << __FILE__ << ":" << __LINE__ << std::endl; \
    internal_abort();                                                   \
  }

#define POLY_CHECK2(x, msg)                                             \
  if (!(x)) {                                                           \
    std::cout << "Check " << #x << " failed at " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    internal_abort();                                                   \
  }

} //end namespace polytope

#endif
