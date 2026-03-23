# ----------------------------------------------------------------------------------
# CMake file for Polytope's Python bindings using PYB11Generator
#   -Builds and compiles the Python bindings according to UsePyBindGen cmake module
# ----------------------------------------------------------------------------------

# The serial tessellators to be wrapped
set(SERIAL_TESSELLATORS "")
if(POLYTOPE_ENABLE_BOOST)
  list(APPEND SERIAL_TESSELLATORS "BoostTessellator2d")
endif()
if(POLYTOPE_ENABLE_TETGEN)
  list(APPEND SERIAL_TESSELLATORS "TetgenTessellator3d")
endif()

# Generate the Tessellators.py script with the tessellators we're building
configure_file(Tessellators.py.in ${PROJECT_BINARY_DIR}/src/PYB11/Tessellators.py)

# Add our PYB11 module
PYB11Generator_add_module(polytope
  MODULE   polytope
  SOURCE   polytope_PYB11.py
  INCLUDES ${CMAKE_INSTALL_PREFIX}/include/polytope
  LINKS    polytopeC
  USE_BLT  ON
  INSTALL  ${POLYTOPE_SITE_PACKAGES_PATH}
  DEPENDS  polytopeC
)
