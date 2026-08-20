//----------------------------------------------------------------------------//
// Wrappers for the ac_types library and some helper functions.
//----------------------------------------------------------------------------//
#ifndef __Polytope_ac_wrapper__
#define __Polytope_ac_wrapper__

#include "ac_int.h"
#include <type_traits>
#include <cstdint>
#include <cstddef>

// Routines for testing if a type is an ac_int
template<typename T>
struct is_ac_int : std::false_type {};

template<int W, bool Signed>
struct is_ac_int<ac_int<W, Signed>> : std::true_type {};

template<typename T>
inline constexpr bool is_ac_int_v =
  is_ac_int<std::remove_cv_t<std::remove_reference_t<T>>>::value;

// ac types does will only convert back to generic types
// with explicit calls
template<typename From, typename To>
To ac_converter(const From& in) {
  if constexpr (is_ac_int_v<From>) {
    if constexpr (std::is_same_v<To, double>) {
      return in.to_double();
    } else if constexpr (std::is_same_v<To, int>) {
      return in.to_int();
    } else if constexpr (std::is_same_v<To, long int>) {
      return in.to_long();
    } else if constexpr (std::is_same_v<To, int64_t>) {
      return in.to_int64();
    } else if constexpr (std::is_same_v<To, unsigned int>) {
      return in.to_uint();
    } else if constexpr (std::is_same_v<To, unsigned long int>) {
      return in.to_ulong();
    } else if constexpr (std::is_same_v<To, uint64_t>) {
      return in.to_uint64();
    } else {
      return To(in);
    }
  } else {
    return static_cast<To>(in);
  }
}

#endif
