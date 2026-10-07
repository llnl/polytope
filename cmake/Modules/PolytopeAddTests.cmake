# CMake macros for adding tests into Polytope
# Consult test/CMakeLists.txt for usage
###############################################################

#--------------------------------------------------------------
# polytope_add_test
# Add a test to Polytope
#--------------------------------------------------------------
macro(polytope_add_test target)
  set(options )
  set(singleValueArgs NUMTASKS)
  set(multiValueArgs )

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED arg_NUMTASKS OR "${arg_NUMTASKS}" STREQUAL "")
    set(test_name "${target}_test")
    set(arg_NUMTASKS 1)
    message("--- Creating test ${target}")
  else()
    set(test_name "${target}_${arg_NUMTASKS}_test")
    message("--- Creating test ${target}: N=${arg_NUMTASKS}")
  endif()

  if (NOT TARGET ${target})
    get_property(POLYTOPE_TPL_DEPENDS GLOBAL PROPERTY POLYTOPE_TPL_DEPENDS)
    get_property(POLYTOPE_CXX_COMPILE_FLAGS GLOBAL PROPERTY POLYTOPE_CXX_COMPILE_FLAGS)
    blt_add_executable(NAME ${target}
      SOURCES ${target}.cc
      DEPENDS_ON polytopeC ${POLYTOPE_TPL_DEPENDS}
      OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin
      OUTPUT_NAME ${target}
    )

    target_compile_options(${target} PRIVATE ${POLYTOPE_CXX_COMPILE_FLAGS})
    target_include_directories(${target} SYSTEM PRIVATE
      ${POLYTOPE_ROOT_DIR}/tests
      ${POLYTOPE_ROOT_DIR}/src
      ${POLYTOPE_ROOT_DIR}/src/Partitioners
      ${PROJECT_BINARY_DIR}/src
    )
  endif()
  blt_add_test(NAME ${test_name}
    COMMAND ${CMAKE_BINARY_DIR}/bin/${target}
    NUM_MPI_TASKS ${arg_NUMTASKS}
  )
  set_tests_properties(${test_name} PROPERTIES
    FIXTURES_REQUIRED "polytope_fixture"
    WORKING_DIRECTORY "${TEST_WORK_DIR}"
  )
endmacro()

macro(polytope_add_python_test target)
  set(options )
  set(singleValueArgs NUMTASKS)
  set(multiValueArgs )
  set(target "${target}_python")

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED arg_NUMTASKS OR "${arg_NUMTASKS}" STREQUAL "")
    set(test_name "${target}_test")
    set(arg_NUMTASKS 1)
    message("--- Creating test ${target}")
  else()
    set(test_name "${target}_{arg_NUMTASKS}_test")
    message("--- Creating test ${target}: N=${arg_NUMTASKS}")
  endif()

  blt_add_test(NAME ${test_name}
    COMMAND ${CMAKE_BINARY_DIR}/${POLYTOPE_VIRT_DIR}/bin/python ${CMAKE_CURRENT_SOURCE_DIR}/${target}.py
    NUM_MPI_TASKS ${arg_NUMTASKS}
  )

  set_tests_properties(${test_name} PROPERTIES
    FIXTURES_REQUIRED "polytope_fixture"
    WORKING_DIRECTORY "${TEST_WORK_DIR}"
  )
endmacro()
