
# Build a configuration header file from our options.
include_directories(.)
include_directories(${PROJECT_BINARY_DIR}/newsrc)
configure_file(
  "${PROJECT_SOURCE_DIR}/newsrc/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/newsrc/polytope.hh"
)

set(polytope_headers
  Cell.hh
  Clipping2D.hh
  Communicator.hh
  EdgeUtils.hh
  GeomUtils.hh
  HashKey.hh
  Intersections.hh
  PLC.hh
  ReducedPLC.hh
  RegisterBoostPolygonTypes.hh
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

list(APPEND polytope_headers "${PROJECT_BINARY_DIR}/newsrc/polytope.hh")

# Include Silo stuff if we've found it.
if (POLYTOPE_ENABLE_SILO)
  list(APPEND polytope_headers
    SiloReader.hh
    SiloUtils.hh
    SiloWriter.hh)
  list(APPEND polytope_sources
    SiloWriter_2d.cc
    SiloWriter_3d.cc
    SiloReader_2d.cc
    SiloReader_3d.cc
    SiloUtils.cc)
endif()

# Check if we can include each of the serial tessellators:
# Boost tessellator
if (POLYTOPE_ENABLE_BOOST)
  list(APPEND polytope_headers
    BoostTessellator.hh)
  list(APPEND polytope_sources
    BoostTessellator.cc)
endif()

# Tetgen tessellator
# if (POLYTOPE_ENABLE_TETGEN)
#   list(APPEND polytope_sources TetgenTessellator.cc)
# endif ()

# Triangle tessellator
if (POLYTOPE_ENABLE_TRIANGLE)
  list(APPEND polytope_headers TriangleTessellator.hh)
  list(APPEND polytope_sources TriangleTessellator.cc)
  if (NOT POLYTOPE_ENABLE_TETGEN)
    list(APPEND polytope_source predicates.cc)
  endif()
endif()

# VORO++ tessellator
if (POLYTOPE_ENABLE_VORO)
  # list(APPEND polytope_headers VoroPP_2d.hh VoroPP_3d.hh)
  # list(APPEND polytope_sources VoroPP_2d.cc VoroPP_3d.cc)
endif()

# Check if we can incude the distributed tessellators:
if (POLYTOPE_ENABLE_MPI)
  list(APPEND polytope_headers
    DistributedTessellator.hh)
  list(APPEND polytope_sources
    DistributedTessellator.cc)
endif()

get_property(POLYTOPE_TPL_DEPENDS GLOBAL PROPERTY POLYTOPE_TPL_DEPENDS)
get_property(POLYTOPE_CXX_COMPILE_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_COMPILE_FLAGS)
get_property(POLYTOPE_CXX_LINK_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_LINK_FLAGS)
blt_add_library(NAME newpolytopeC
  HEADERS    ${polytope_headers}
  SOURCES    ${polytope_sources}
  DEPENDS_ON ${POLYTOPE_TPL_DEPENDS}
  SHARED     ${BUILD_SHARED_LIBS}
)

target_compile_options(newpolytopeC PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})
target_link_options(newpolytopeC PUBLIC ${POLYTOPE_CXX_LINK_FLAGS})

# We must set the polytope python target to be "polytope"
# so the module is named polytope. So our C++ target is called
# newpolytopeC. But we will export it as polytope
set_target_properties(newpolytopeC PROPERTIES OUTPUT_NAME "newpolytope")
set_target_properties(newpolytopeC PROPERTIES EXPORT_NAME "newpolytope")
# Library install targets
install(TARGETS newpolytopeC EXPORT polytope-targets DESTINATION lib)
install(FILES ${polytope_headers} DESTINATION include/polytope)

if(POLYTOPE_ENABLE_PYTHON)
  add_subdirectory(PYB11)
endif()
