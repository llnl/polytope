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

if (NOT MSVC)
  list(APPEND POLYTOPE_CXX_COMPILE_FLAGS -fPIC)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
  list(APPEND POLYTOPE_CXX_COMPILE_FLAGS -qPIC)
endif()
