# Minimum CMake version.
cmake_minimum_required (VERSION 3.21)
include(CMakeDependentOption)

# Set Polytope's module path.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(cmake/PolytopeVersion.cmake)
project(Polytope LANGUAGES C CXX VERSION ${POLYTOPE_VERSION})

# We require at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(POLYTOPE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(POLYTOPE_VERSION_MINOR ${PROJECT_VERSION_MINOR})

set(POLYTOPE_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "Path to top source directory")

#-------------------------------------------------------------------------------
# Set full rpath information by default
#-------------------------------------------------------------------------------
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# List containing target dependencies
set(POLYTOPE_TPL_DEPENDS )

#option(POLYTOPE_ENABLE_HIP "Build with HIP enabled" OFF) # Currently unsupported
#set(ENABLE_MPI OFF CACHE BOOL "")
set(POLYTOPE_ENABLE_MPI ${ENABLE_MPI})
option(POLYTOPE_ENABLE_BOOST "Enable Boost Voronoi tessellator" OFF)
option(POLYTOPE_ENABLE_SILO "Enable Silo read/write features" OFF)
option(POLYTOPE_ENABLE_PYTHON "Enable Python interface" OFF)
option(POLYTOPE_ENABLE_TRIANGLE "Enable Triangle Voronoi tessellator" OFF)
option(POLYTOPE_ENABLE_TETGEN "Enable TetGen Voronoi tessellator" OFF)
set(TESTING_PREREQS (POLYTOPE_ENABLE_BOOST AND POLYTOPE_ENABLE_SILO))
option(POLYTOPE_ENABLE_TESTING "Enable testing" OFF)
option(POLYTOPE_ENABLE_ASAN "Enable ASAN" OFF)
option(POLYTOPE_ENABLE_UBSAN "Enable UBSAN" OFF)
if (POLYTOPE_ENABLE_TEST AND NOT TESTING_PREREQS)
  message(FATAL_ERROR "Testing relies on Boost and Silo. Must enable Boost and Silo to enable testing.")
endif()
# Voro is not supported
#option(POLYTOPE_ENABLE_VORO "Enable Voro Voronoi tessellator" OFF)
option(POLYTOPE_ENABLE_APPLE "Enable Apple installation support" OFF)
option(ENABLE_STATIC_TPLS "Used for libraries that require older HDF5 versions" OFF)
option(BUILD_SHARED_LIBS "Build Polytope as shared library" OFF)

include(cmake/Compiler.cmake)

# Configure and include BLT
if (POLYTOPE_ENABLE_PYTHON AND DEFINED python_DIR)
  set(Python_EXECUTABLE ${python_DIR}/bin/python3)
  set(Python3_EXECUTABLE ${python_DIR}/bin/python3)
endif()
set(BLT_CXX_STD "c++17" CACHE STRING "")

if(NOT EXISTS "${POLYTOPE_ROOT_DIR}/cmake/blt/SetupBLT.cmake")
  message(FATAL_ERROR "No BLT instance found. Update submodules.")
endif()

include(cmake/blt/SetupBLT.cmake)

#-------------------------------------------------------------------------------#
# Set a default build type if none was specified
#-------------------------------------------------------------------------------#
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build (debug, release, etc)." FORCE)

  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(POLYTOPE_ENABLE_DEBUG ON)
endif()

include(cmake/InstallTPLs.cmake)

# old installs do not work with MPI enabled
if (NOT ENABLE_MPI)
  add_subdirectory(src)
endif()
add_subdirectory(newsrc)

if (POLYTOPE_ENABLE_TESTING)
  enable_testing()
  #add_subdirectory(tests)
  add_subdirectory(newtests)
endif()

# Export package to allow it to be installed with find_package
install(EXPORT polytope-targets DESTINATION lib/cmake)
include(CMakePackageConfigHelpers)
set(POLYTOPE_INSTALL_INCLUDE_DIR "include/polytope" CACHE STRING "")
set(POLYTOPE_INSTALL_LIB_DIR "lib" CACHE STRING "")
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/polytope-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/polytope-config.cmake
  INSTALL_DESTINATION "lib"
  PATH_VARS
  POLYTOPE_INSTALL_INCLUDE_DIR
  POLYTOPE_INSTALL_LIB_DIR)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polytope-config.cmake
  DESTINATION "lib/cmake")
# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE
  "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
include(CPack)
