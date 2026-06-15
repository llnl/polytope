# CMake macros for adding tests into Polytope
# Consult test/CMakeLists.txt for usage
###############################################################

#--------------------------------------------------------------
# polytope_add_test
# Add a test to Polytope
#--------------------------------------------------------------
# TEST_WORK_DIR must be defined in calling CMake file
macro(polytope_add_test target)
  set(options )
  set(singleValueArgs NUMTASKS)
  set(multiValueArgs )

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED arg_NUMTASKS OR "${arg_NUMTASKS}" STREQUAL "")
    set(arg_NUMTASKS 1)
    message("--- Creating test ${target}")
  else()
    message("--- Creating test ${target}: N=${arg_NUMTASKS}")
  endif()

  get_property(POLYTOPE_TPL_DEPENDS GLOBAL PROPERTY POLYTOPE_TPL_DEPENDS)
  get_property(POLYTOPE_CXX_COMPILE_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_COMPILE_FLAGS)
  blt_add_executable(NAME ${target}
    SOURCES test_${target}.cc
    DEPENDS_ON polytopeC ${POLYTOPE_TPL_DEPENDS}
    OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin
  )

  target_compile_options(${target} PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})
  target_include_directories(${target} SYSTEM PRIVATE
    ${POLYTOPE_ROOT_DIR}/tests
    ${POLYTOPE_ROOT_DIR}/src
    ${PROJECT_BINARY_DIR}/src
  )

  blt_add_test(NAME ${target}_test
    COMMAND ${CMAKE_BINARY_DIR}/bin/${target}
    NUM_MPI_TASKS ${arg_NUMTASKS}
  )
  set_tests_properties(${target}_test PROPERTIES
    FIXTURES_REQUIRED "polytope_fixture"
    WORKING_DIRECTORY "${TEST_WORK_DIR}"
  )
endmacro()
macro(polytope_add_new_test newtarget)
  set(options )
  set(singleValueArgs NUMTASKS)
  set(multiValueArgs )
  set(target "new${newtarget}")

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED arg_NUMTASKS OR "${arg_NUMTASKS}" STREQUAL "")
    set(arg_NUMTASKS 1)
    message("--- Creating test ${target}")
  else()
    message("--- Creating test ${target}: N=${arg_NUMTASKS}")
  endif()

  get_property(POLYTOPE_TPL_DEPENDS GLOBAL PROPERTY POLYTOPE_TPL_DEPENDS)
  get_property(POLYTOPE_CXX_COMPILE_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_COMPILE_FLAGS)
  blt_add_executable(NAME ${target}
    SOURCES test_${newtarget}.cc
    DEPENDS_ON newpolytopeC ${POLYTOPE_TPL_DEPENDS}
    OUTPUT_DIR ${CMAKE_BINARY_DIR}/newbin
    OUTPUT_NAME ${newtarget}
  )

  target_compile_options(${target} PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})
  target_include_directories(${target} SYSTEM PRIVATE
    ${POLYTOPE_ROOT_DIR}/newtests
    ${POLYTOPE_ROOT_DIR}/newsrc
    ${PROJECT_BINARY_DIR}/newsrc
  )

  blt_add_test(NAME ${target}_test
    COMMAND ${CMAKE_BINARY_DIR}/newbin/${newtarget}
    NUM_MPI_TASKS ${arg_NUMTASKS}
  )
  set_tests_properties(${target}_test PROPERTIES
    FIXTURES_REQUIRED "newpolytope_fixture"
    WORKING_DIRECTORY "${NEW_TEST_WORK_DIR}"
  )
endmacro()
