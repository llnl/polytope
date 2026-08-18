#-----------------------------------------------------------------------------------
# Polytope TPLs
#-----------------------------------------------------------------------------------

# MPI
#-----------------------------------------------------------------------------------
if (POLYTOPE_ENABLE_MPI)
  find_package(MPI REQUIRED)
  list(APPEND POLYTOPE_TPL_DEPENDS MPI::MPI_C MPI::MPI_CXX)

  if (MPIEXEC_EXECUTABLE)
    set(POLYTOPE_MPIEXEC "${MPIEXEC_EXECUTABLE}" CACHE FILEPATH "MPI launcher for Polytope tests")
  elseif (MPIEXEC)
    set(POLYTOPE_MPIEXEC "${MPIEXEC}" CACHE FILEPATH "MPI launcher for Polytope tests")
  endif()

  if (MPIEXEC_NUMPROC_FLAG)
    set(POLYTOPE_MPIEXEC_NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING "MPI launcher process-count flag")
  else()
    set(POLYTOPE_MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "MPI launcher process-count flag")
  endif()
endif()

# Manually imported libs
set(IMPORTED_LIBS )

# Libs brought in through find_package
set(FP_TPLS )
set(FP_DIRS )

# Python
#-----------------------------------------------------------------------------------
if (POLYTOPE_ENABLE_PYTHON)
  # Find the appropriate Python
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  set(POLYTOPE_SITE_PACKAGES_PATH "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages" )
  list(APPEND POLYTOPE_TPL_DEPENDS Python3::Python)

  # Set the PYB11Generator path
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${POLYTOPE_ROOT_DIR}/extern/PYB11Generator" CACHE PATH "")
  endif()
  # Set the pybind11 path
  if (NOT PYBIND11_ROOT_DIR)
    set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
    set(PYBIND11_NOPYTHON TRUE)
  endif()
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)
  list(APPEND POLYTOPE_TPL_DEPENDS pybind11_headers)
  list(APPEND IMPORTED_LIBS pybind11_headers)
endif()

# BOOST
#-----------------------------------------------------------------------------------
if(POLYTOPE_ENABLE_BOOST AND NOT BOOST_FOUND)
  if (DEFINED boost_DIR)
    find_package(Boost 1.50 REQUIRED NO_DEFAULT_PATH PATHS ${boost_DIR})
  else()
    message(FATAL_ERROR "Must provide boost_DIR to enable Boost")
  endif()

  if (TARGET Boost::headers)
    list(APPEND POLYTOPE_TPL_DEPENDS Boost::headers)
  elseif (TARGET Boost::boost)
    list(APPEND POLYTOPE_TPL_DEPENDS Boost::boost)
  else()
    add_library(polytope_boost_headers INTERFACE)
    target_include_directories(polytope_boost_headers SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
    list(APPEND POLYTOPE_TPL_DEPENDS polytope_boost_headers)
  endif()
endif()

# Silo/HDF5
#-----------------------------------------------------------------------------------
if(POLYTOPE_ENABLE_SILO)
  if(NOT ENABLE_STATIC_TPL)
    list(APPEND FP_TPLS hdf5)
    list(APPEND FP_DIRS ${hdf5_DIR})
    find_package(hdf5 REQUIRED NO_DEFAULT_PATH PATHS ${hdf5_DIR})
    message("Found HDF5 External Package.")
    if(ENABLE_STATIC_TPL)
      list(APPEND POLYTOPE_TPL_DEPENDS hdf5-static hdf5_hl-static)
    else()
      list(APPEND POLYTOPE_TPL_DEPENDS hdf5-shared hdf5_hl-shared)
    endif()
  else()
    set(HDF5_LIBS ${hdf5_DIR}/lib/libhdf5.a ${hdf5_DIR}/lib/libhdf5_hl.a)
    blt_import_library(NAME hdf5
      LIBRARIES ${HDF5_LIBS}
      INCLUDES ${hdf5_DIR}/include
      TREAT_INCLUDES_AS_SYSTEM ON
      EXPORTABLE ON)
    list(APPEND POLYTOPE_TPL_DEPENDS hdf5)
    list(APPEND IMPORTED_LIBS hdf5)
  endif()
  # find_package(Silo... is currently broken but should be working by the next release
  set(SILO_LIB_NAME libsiloh5.a)
  if (POLYTOPE_ENABLE_APPLE)
    set(SILO_LIB_NAME libsiloh5.dylib)
  endif()
  file(GLOB_RECURSE SILO_LIB "${silo_DIR}/*${SILO_LIB_NAME}")
  blt_import_library(NAME silo
    LIBRARIES ${SILO_LIB}
    TREAT_INCLUDES_AS_SYSTEM ON
    INCLUDES ${silo_DIR}/include
    EXPORTABLE ON)
  list(APPEND POLYTOPE_TPL_DEPENDS silo)
  list(APPEND IMPORTED_LIBS silo)
endif()

# Qhull
#-----------------------------------------------------------------------------------
if(POLYTOPE_ENABLE_QHULL AND NOT Qhull_FOUND)
  list(APPEND FP_TPLS Qhull)
  list(APPEND FP_DIRS ${qhull_DIR})
  find_package(Qhull REQUIRED NO_DEFAULT_PATH PATHS ${qhull_DIR})
  list(APPEND POLYTOPE_TPL_DEPENDS Qhull::qhull_r Qhull::qhullcpp)
endif()

# Triangle
#-----------------------------------------------------------------------------------
# Spack does not install Triangle in any useful way so we have to install it ourselves
# To use it, make sure the triangle source code is in ${triangle_SRC_DIR}.
if(POLYTOPE_ENABLE_TRIANGLE)
  if(NOT DEFINED triangle_SRC_DIR)
    message(FATAL_ERROR "Must provide triangle_SRC_DIR to enable Triangle")
  endif()
  set(triangle_sources ${triangle_SRC_DIR}/triangle.c)
  set(triangle_headers ${triangle_SRC_DIR}/triangle.h)
  blt_add_library(NAME triangle
    HEADERS ${triangle_headers}
    SOURCES ${triangle_sources}
    SHARED TRUE)
  target_compile_definitions(triangle PRIVATE TRILIBRARY ANSI_DECLARATORS CDT_ONLY)
  target_compile_definitions(triangle PRIVATE REAL=double VOID=void)
  install(TARGETS triangle EXPORT polytope-targets DESTINATION lib)
  install(FILES ${triangle_headers} DESTINATION include/triangle)
  list(APPEND POLYTOPE_TPL_DEPENDS triangle)
  include_directories(${triangle_SRC_DIR})
endif()

# Tetgen
#-----------------------------------------------------------------------------------
if(POLYTOPE_ENABLE_TETGEN)
  blt_import_library(NAME tetgen
    LIBRARIES ${tetgen_DIR}/lib/libtet.a
    INCLUDES ${tetgen_DIR}/include
    TREAT_INCLUDES_AS_SYSTEM ON
    EXPORTABLE ON)
  list(APPEND POLYTOPE_TPL_DEPENDS tetgen)
  list(APPEND IMPORTED_LIBS tetgen)
endif()

# Caliper
#-----------------------------------------------------------------------------------
# TODO: add Caliper support
# if(POLYTOPE_ENABLE_TIMERS)
#   if(NOT DEFINED caliper_DIR)
#     message(FATAL_ERROR "Must provide caliper_DIR if enabling timers")
#   endif()
#   list(APPEND FP_TPLS caliper)
#   list(APPEND FP_DIRS ${caliper_DIR})
#   find_package(caliper REQUIRED NO_DEFAULT_PATHS PATHS ${caliper_DIR}/share/cmake/caliper)
#   list(APPEND POLYTOPE_TPL_DEPENDS caliper)
# endif()

foreach(lib ${IMPORTED_LIBS})
  get_target_property(_is_imported ${lib} IMPORTED)
  if(NOT ${_is_imported})
    install(TARGETS ${lib}
      EXPORT polytope-targets
      DESTINATION lib/cmake)
    set_target_properties(${lib} PROPERTIES EXPORT_NAME polytope::${lib})
  endif()
endforeach()

set_property(GLOBAL APPEND PROPERTY POLYTOPE_TPL_DEPENDS "${POLYTOPE_TPL_DEPENDS}")
