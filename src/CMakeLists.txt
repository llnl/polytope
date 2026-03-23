
# Build a configuration header file from our options.
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/src/polytope.hh"
)
install(FILES ${PROJECT_BINARY_DIR}/src/polytope.hh
  DESTINATION include)

set(polytope_sources
  KeyTraits.cc
  predicates.cc
  PLC_CSG.cc
  QuantizedTessellation2d.cc
  QuantizedTessellation3d.cc
  clipQuantizedTessellation2d.cc
  clipQuantizedTessellation3d.cc
  snapToBoundary.cc
  ErrorHandler.cc
  polytope_internal_abort.cc)

file(GLOB polytope_headers CONFIGURE_DEPENDS "*.hh")
#list(APPEND polytope_headers "${PROJECT_BINARY_DIR}/src/polytope.hh")

# Include Silo stuff if we've found it.
if (POLYTOPE_ENABLE_SILO)
  list(APPEND polytope_sources
    SiloWriter_2d.cc
    SiloWriter_3d.cc
    SiloReader_2d.cc
    SiloReader_3d.cc
    SiloUtils.cc)
endif()

# Check if we can include each of the serial tessellators:
if (POLYTOPE_ENABLE_BOOST)
  list(APPEND polytope_sources
    BoostTessellator.cc
    MeshEditor.cc)
endif()

# TetgenTessellator
if (POLYTOPE_ENABLE_TETGEN)
  list(APPEND polytope_sources TetgenTessellator.cc)
endif ()

# Check if we can incude the distributed tessellators:
if (POLYTOPE_ENABLE_MPI)
  list(APPEND polytope_sources
    DistributedTessellator.cc
    SerialDistributedTessellator.cc)
endif()

get_property(POLYTOPE_TPL_DEPENDS GLOBAL PROPERTY POLYTOPE_TPL_DEPENDS)
blt_add_library(NAME polytopeC
  HEADERS    ${polytope_headers}
  SOURCES    ${polytope_sources}
  DEPENDS_ON ${POLYTOPE_TPL_DEPENDS}
  SHARED     ON
)

target_compile_options(polytopeC PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})

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
