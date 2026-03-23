# Minimum CMake version.
cmake_minimum_required (VERSION 3.21)

# Set Polytope's module path.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(cmake/PolytopeVersion.cmake)
project(Polytope LANGUAGES C CXX VERSION ${POLYTOPE_VERSION})

# We require at least C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(POLYTOPE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(POLYTOPE_VERSION_MINOR ${PROJECT_VERSION_MINOR})

set(POLYTOPE_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "Path to top source directory")

# List containing target dependencies
set(POLYTOPE_TPL_DEPENDS )

#option(POLYTOPE_ENABLE_HIP "Build with HIP enabled" OFF) # Currently unsupported
option(POLYTOPE_ENABLE_MPI "Build with MPI enabled" OFF)
option(POLYTOPE_ENABLE_BOOST "Enable Boost Voronoi tessellator" ON)
option(POLYTOPE_ENABLE_SILO "Enable Silo read/write features" OFF)
option(POLYTOPE_ENABLE_PYTHON "Enable Python interface" OFF)
option(POLYTOPE_ENABLE_TRIANGLE "Enable Triangle Voronoi tesselator" OFF)
option(POLYTOPE_ENABLE_TETGEN "Enable TetGen Voronoi tesselator" OFF)
option(POLYTOPE_ENABLE_TESTING "Enable testing" ON)
option(POLYTOPE_ENABLE_C "Enable C interface" OFF)
option(ENABLE_STATIC_TPLS "Used for libraries that require older HDF5 versions" OFF)
option(BUILD_SHARED_LIBS "Build Polytope as shared library" OFF)

include(cmake/Compiler.cmake)

# Configure and include BLT
if (POLYTOPE_ENABLE_PYTHON AND DEFINED python_DIR)
  set(Python_EXECUTABLE ${python_DIR}/bin/python3)
  set(Python3_EXECUTABLE ${python_DIR}/bin/python3)
endif()
set(BLT_CXX_STD "c++14" CACHE STRING "")

if(NOT EXISTS "${POLYTOPE_ROOT_DIR}/cmake/blt/SetupBLT.cmake")
  message(FATAL_ERROR "No BLT instance found. Update submodules.")
endif()

include(cmake/blt/SetupBLT.cmake)

include(cmake/InstallTPLs.cmake)

add_subdirectory(src)

if (POLYTOPE_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if (POLYTOPE_ENABLE_C)
  add_subdirectory(polytope_c)
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
