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
if(POLYTOPE_ENABLE_BOOST)
  if (DEFINED boost_DIR)
    list(APPEND FP_DIRS ${boost_DIR})
    list(APPEND FP_TPLS Boost)
    find_package(Boost 1.50 REQUIRED NO_DEFAULT_PATH PATHS ${boost_DIR})
  else()
    message(FATAL_ERROR "Must provide boost_DIR if enabling Boost")
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
  # find_package(Silo is currently broken but should be working by the next release
  set(SILO_LIB_NAME libsiloh5.a)
  if (POLYTOPE_ENABLE_APPLE)
    set(SILO_LIB_NAME libsiloh5.dylib)
  endif()
  file(GLOB SILO_LIB "${silo_DIR}/*/${SILO_LIB_NAME}")
  blt_import_library(NAME silo
    LIBRARIES ${SILO_LIB}
    TREAT_INCLUDES_AS_SYSTEM ON
    INCLUDES ${silo_DIR}/include
    EXPORTABLE ON)
  list(APPEND POLYTOPE_TPL_DEPENDS silo)
  list(APPEND IMPORTED_LIBS silo)
endif()

if(POLYTOPE_ENABLE_TETGEN)
  blt_import_library(NAME tetgen
    LIBRARIES ${tetgen_DIR}/lib/libtet.a
    INCLUDES ${tetgen_DIR}/include
    TREAT_INCLUDES_AS_SYSTEM ON
    EXPORTABLE ON)
  list(APPEND POLYTOPE_TPL_DEPENDS tetgen)
  list(APPEND IMPORTED_LIBS tetgen)
endif()

if(POLYTOPE_ENABLE_VORO)
  blt_import_library(NAME voro
    LIBRARIES ${voro_DIR}/lib/libvoro++.so
    INCLUDES ${voro_DIR}/include
    TREAT_INCLUDES_AS_SYSTEM ON
    EXPORTABLE ON)
  list(APPEND POLYTOPE_TPL_DEPENDS voro)
  list(APPEND IMPORTED_LIBS voro)
endif()

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
