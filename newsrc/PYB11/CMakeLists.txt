# ----------------------------------------------------------------------------------
# Python bindings for the renovated Polytope C++ implementation using PYB11Generator.
# ----------------------------------------------------------------------------------

set(SERIAL_TESSELLATORS "")
if(POLYTOPE_ENABLE_BOOST)
  list(APPEND SERIAL_TESSELLATORS "BoostTessellator")
endif()
if(POLYTOPE_ENABLE_TRIANGLE)
  list(APPEND SERIAL_TESSELLATORS "TriangleTessellator")
endif()

configure_file(
  SerialTessellators.py.in
  ${PROJECT_BINARY_DIR}/newsrc/PYB11/SerialTessellators.py
)

PYB11Generator_add_module(polytope
  MODULE   polytope
  SOURCE   polytope_PYB11.py
  INCLUDES ${PROJECT_SOURCE_DIR}/newsrc
           ${PROJECT_BINARY_DIR}/newsrc
           ${CMAKE_CURRENT_SOURCE_DIR}
  LINKS    newpolytopeC
  USE_BLT  ON
  INSTALL  OFF
  DEPENDS  newpolytopeC
)

install(TARGETS polytope
  DESTINATION ${POLYTOPE_SITE_PACKAGES_PATH})

include(PolytopePythonEnvironment)

set(POLYTOPE_VIRT_DIR venv)

polytope_add_python_environment(
  build
  "${CMAKE_BINARY_DIR}/${POLYTOPE_VIRT_DIR}"
  "${LIBRARY_OUTPUT_PATH}"
  ALL
  TARGET_NAME polytope_build_virt_env
  DEPENDS Python3::Python newpolytopeC polytope
)

polytope_add_python_environment(
  install
  "${CMAKE_INSTALL_PREFIX}/${POLYTOPE_VIRT_DIR}"
  "${CMAKE_INSTALL_PREFIX}/${POLYTOPE_SITE_PACKAGES_PATH}"
  TARGET_NAME polytope_install_virt_env
  DEPENDS Python3::Python newpolytopeC polytope
)

# Add the virtual env target run to the install
install(CODE "
message(\"-- Setting up Polytope virtual environment.\")
execute_process(
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target polytope_install_virt_env
  RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR \"Failed to set up Polytope virtual environment\")
endif()
")

# Make a script that generates a wrapper file that allows us
# to run the virtual python environment directly

set(WRAPPER_TEMPLATE
  [=[#!/usr/bin/env bash
# Determine the path of this script file
POLYTOPE_BIN_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Adding -B to prevent byte-compiling
$POLYTOPE_BIN_DIR/../@POLYTOPE_VIRT_DIR@/bin/python -B "$@"
]=]
)

string(CONFIGURE
  "${WRAPPER_TEMPLATE}"
  WRAPPER_CONTENT
  @ONLY
)

set(POLYTOPE_PYTHON_WRAPPER_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/wrapper_script)

file(GENERATE
  OUTPUT "${POLYTOPE_PYTHON_WRAPPER_SCRIPT}"
  CONTENT "${WRAPPER_CONTENT}"
)

install(
  PROGRAMS ${POLYTOPE_PYTHON_WRAPPER_SCRIPT}
  DESTINATION bin
  RENAME python
)
