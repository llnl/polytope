# Performance benchmarks are intentionally separate from CTest correctness
# tests and are not installed.

macro(cxx_benchmark name input_value mpi_tasks)
  blt_add_executable(NAME ${name}
    SOURCES ${name}.cc
    DEPENDS_ON polytopeC ${POLYTOPE_TPL_DEPENDS}
    OUTPUT_DIR ${CMAKE_BINARY_DIR}/benchmarks
  )

  target_compile_options(${name} PRIVATE
    ${POLYTOPE_CXX_COMPILE_FLAGS})
  target_link_options(${name} PUBLIC
    ${POLYTOPE_CXX_LINK_FLAGS})
  target_include_directories(${name} SYSTEM PRIVATE
    ${POLYTOPE_ROOT_DIR}/tests
    ${POLYTOPE_ROOT_DIR}/src
    ${POLYTOPE_ROOT_DIR}/src/Partitioners
    ${PROJECT_BINARY_DIR}/src
  )
  blt_add_benchmark(NAME ${name}CXX
    COMMAND ${CMAKE_BINARY_DIR}/benchmarks/${name} ${input_value}
    NUM_MPI_TASKS ${mpi_tasks}
  )
endmacro()

# Add the CXX benchmarks
cxx_benchmark(DistTessellator2D 1000000 36)
cxx_benchmark(SerialTessellator2D 10 1)

# Add the Python benchmarks
blt_add_benchmark(NAME DistTessellator2DPY
  COMMAND ${CMAKE_BINARY_DIR}/${POLYTOPE_VIRT_DIR}/bin/python -u ${CMAKE_CURRENT_SOURCE_DIR}/DistTessellator2D.py 1000000
  NUM_MPI_TASKS 36
)

blt_add_benchmark(NAME DistLloydTessellator2DPY
  COMMAND ${CMAKE_BINARY_DIR}/${POLYTOPE_VIRT_DIR}/bin/python -u ${CMAKE_CURRENT_SOURCE_DIR}/DistLloydTessellator2D.py 100000
  NUM_MPI_TASKS 36
)
