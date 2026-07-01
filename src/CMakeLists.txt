
# Build a configuration header file from our options.
include_directories(.)
include_directories(${PROJECT_BINARY_DIR}/src)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/src/polytope.hh"
)

set(polytope_headers
  bisectSearch.hh
  checkDistributedTessellation.hh
  clipQuantizedTessellation.hh
  convexHull_2d.hh
  convexHull_3d.hh
  convexIntersect.hh
  convexWithin.hh
  deleteCells.hh
  DimensionTraits.hh
  DistributedTessellator.hh
  ErrorHandler.hh
  findBoundaryElements.hh
  Intersections.hh
  intersect.hh
  IntPointMap.hh
  KeyTraits.hh
  makeBoxPLC.hh
  MeshEditor.hh
  mortonOrderIndices.hh
  nearestPoint.hh
  OrphanageBase.hh
  PLC_CSG_2d.hh
  PLC_CSG_3d.hh
  PLC.hh
  Point.hh
  Polygon.hh
  polytope_geometric_utilities.hh
  polytope_internal.hh
  polytope_plc_canned_geometries.hh
  polytope_serialize.hh
  polytope_tessellator_utilities.hh
  polytope_write_OOGL.hh
  QuantizedCoordinates.hh
  QuantizedTessellation2d.hh
  QuantizedTessellation3d.hh
  QuantTessellation.hh
  ReducedPLC.hh
  removeElements.hh
  Segment.hh
  SerialDistributedTessellator.hh
  simplifyPLCfacets.hh
  snapToBoundary.hh
  Tessellation.hh
  Tessellator.hh
  TessellatorInline.hh
  timingUtilities.hh
  within.hh
)

set(polytope_sources
  KeyTraits.cc
  PLC_CSG.cc
  QuantizedTessellation2d.cc
  QuantizedTessellation3d.cc
  clipQuantizedTessellation2d.cc
  clipQuantizedTessellation3d.cc
  snapToBoundary.cc
  ErrorHandler.cc
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
    Clipper2d.hh
    PLC_Boost_2d.hh
    RegisterBoostPolygonTypes.hh
    BoostOrphanage.hh
    BoostTessellator.hh
    BoostTessellatorTraits.hh)
  list(APPEND polytope_sources
    BoostTessellator.cc
    MeshEditor.cc)
endif()

# Tetgen tessellator
if (POLYTOPE_ENABLE_TETGEN)
  list(APPEND polytope_headers TetgenTessellator.hh)
  list(APPEND polytope_sources TetgenTessellator.cc)
else()
  list(APPEND polytope_sources predicates.cc)
endif ()

# Triangle tessellator
if (POLYTOPE_ENABLE_TRIANGLE)
  list(APPEND polytope_headers TriangleTessellator.hh)
  list(APPEND polytope_sources TriangleTessellator.cc)
endif()

# VORO++ tessellator
if (POLYTOPE_ENABLE_VORO)
  # list(APPEND polytope_headers VoroPP_2d.hh VoroPP_3d.hh)
  # list(APPEND polytope_sources VoroPP_2d.cc VoroPP_3d.cc)
endif()

# Check if we can incude the distributed tessellators:
if (POLYTOPE_ENABLE_MPI)
  list(APPEND polytope_headers
    polytope_parallel_utilities.hh)
  list(APPEND polytope_sources
    DistributedTessellator.cc
    SerialDistributedTessellator.cc)
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

if(POLYTOPE_ENABLE_PYTHON)
  add_subdirectory(PYB11)
endif()
