set(LANG_STR "CXX")
if(POLYTOPE_ENABLE_HIP)
  set(LANG_STR "HIP")
endif()

# General compiler flags.
set(POLYTOPE_CXX_COMPILE_FLAGS
  -Wall
  -Wno-sign-compare
  -Wno-unused-local-typedefs
  -Wno-misleading-indentation)

set(POLYTOPE_CXX_LINK_FLAGS)

if (NOT MSVC)
  list(APPEND POLYTOPE_CXX_COMPILE_FLAGS -fPIC)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
  list(APPEND POLYTOPE_CXX_COMPILE_FLAGS -qPIC)
endif()

if (POLYTOPE_ENABLE_ASAN)
  list(APPEND POLYTOPE_CXX_COMPILE_FLAGS -fsanitize=address)
  list(APPEND POLYTOPE_CXX_LINK_FLAGS -fsanitize=address)
endif()

if (POLYTOPE_ENABLE_UBSAN)
  list(APPEND POLYTOPE_CXX_COMPILE_FLAGS -fsanitize=undefined)
  list(APPEND POLYTOPE_CXX_LINK_FLAGS -fsanitize=undefined)
endif()

set_property(GLOBAL PROPERTY POLYTOPE_CXX_COMPILE_FLAGS
  "$<$<COMPILE_LANGUAGE:${LANG_STR}>:${POLYTOPE_CXX_COMPILE_FLAGS}>")

set_property(GLOBAL PROPERTY POLYTOPE_CXX_LINK_FLAGS "${POLYTOPE_CXX_LINK_FLAGS}")
