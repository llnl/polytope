
# Build a configuration header file from our options.
include_directories(.)
include_directories(${PROJECT_BINARY_DIR}/src)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/src/polytope.hh"
)

set(polytope_headers
  ac_wrapper.hh
  Cell.hh
  Clipping2D.hh
  Communicator.hh
  Cube.hh
  EdgeUtils.hh
  EscapePod.hh
  GeomUtils.hh
  MortonKeyTraits.hh
  Intersections.hh
  Partitioner.hh
  PLC.hh
  Point.hh
  QuantTessellation.hh
  Quantizer.hh
  QuantPLC.hh
  Shapes.hh
  Tessellation.hh
  Tessellator.hh
  TessellatorInline.hh
  TetgenTessellator.hh
  Serializer.hh
  polytope_internal.hh
  findBoundaryElements.hh
)

set(polytope_sources
  Communicator.cc
  QuantTessellation2D.cc
  QuantTessellation3D.cc
  Tessellation.cc
  QuantPLC.cc
  polytope_internal_abort.cc
)

list(APPEND polytope_headers "${PROJECT_BINARY_DIR}/src/polytope.hh")

# Include Silo stuff if we've found it.
if (POLYTOPE_ENABLE_SILO)
  list(APPEND polytope_headers
    SiloReader.hh
    SiloUtils.hh
    SiloWriter.hh)
  list(APPEND polytope_sources
    SiloWriter2D.cc
    SiloWriter3D.cc
    SiloReader2D.cc
    SiloReader3D.cc)
endif()

# Check if we can include each of the serial tessellators:
# Boost tessellator
if (POLYTOPE_ENABLE_BOOST)
  list(APPEND polytope_headers
    RegisterBoostPolygonTypes.hh
    BoostPolygonIntersections.hh
    BoostTessellator.hh)
  list(APPEND polytope_sources
    BoostTessellator.cc)
endif()

# Tetgen tessellator
if (POLYTOPE_ENABLE_TETGEN)
  list(APPEND polytope_sources TetgenTessellator.cc)
endif ()

# Triangle tessellator
if (POLYTOPE_ENABLE_TRIANGLE)
  list(APPEND polytope_headers TriangleTessellator.hh)
  list(APPEND polytope_sources TriangleTessellator.cc)
  if (NOT POLYTOPE_ENABLE_TETGEN)
    list(APPEND polytope_source predicates.cc)
  endif()
endif()

# Check if we can incude the distributed tessellators:
if (POLYTOPE_ENABLE_MPI)
  list(APPEND polytope_headers
    DistributedTessellator.hh
    ParallelUtils.hh)
  list(APPEND polytope_sources
    DistributedTessellator.cc)
endif()

get_property(POLYTOPE_TPL_DEPENDS GLOBAL PROPERTY POLYTOPE_TPL_DEPENDS)
get_property(POLYTOPE_CXX_COMPILE_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_COMPILE_FLAGS)
get_property(POLYTOPE_CXX_LINK_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_LINK_FLAGS)
blt_add_library(NAME polytopeC
  HEADERS    ${polytope_headers}
  SOURCES    ${polytope_sources}
  DEPENDS_ON ${POLYTOPE_TPL_DEPENDS}
  SHARED     ${BUILD_SHARED_LIBS}
)

target_compile_options(polytopeC PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})
target_link_options(polytopeC PUBLIC ${POLYTOPE_CXX_LINK_FLAGS})

# We must set the polytope python target to be "polytope"
# so the module is named polytope. So our C++ target is called
# polytopeC. But we will export it as polytope
set_target_properties(polytopeC PROPERTIES OUTPUT_NAME "polytope")
set_target_properties(polytopeC PROPERTIES EXPORT_NAME "polytope")
# Library install targets
install(TARGETS polytopeC EXPORT polytope-targets DESTINATION lib)
install(FILES ${polytope_headers} DESTINATION include/polytope)

blt_add_executable(NAME loadEscapePod
  HEADERS ${polytope_headers}
  SOURCES loadEscapePod.cc ${polytope_sources}
  DEPENDS_ON polytopeC ${POLYTOPE_TPL_DEPENDS}
  OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin
  OUTPUT_NAME loadEscapePod
  INCLUDES ${POLYTOPE_ROOT_DIR}/src ${PROJECT_BINARY_DIR}/src
)

target_compile_options(loadEscapePod PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})
target_link_options(loadEscapePod PUBLIC ${POLYTOPE_CXX_LINK_FLAGS})

if(POLYTOPE_ENABLE_PYTHON)
  add_subdirectory(PYB11)
endif()
