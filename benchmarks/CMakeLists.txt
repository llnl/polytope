# Performance benchmarks are intentionally separate from CTest correctness
# tests and are not installed.

blt_add_executable(NAME DistTessellator2D
  SOURCES DistTessellator2D.cc
  DEPENDS_ON polytopeC ${POLYTOPE_TPL_DEPENDS}
  OUTPUT_DIR ${CMAKE_BINARY_DIR}/benchmarks
)

target_compile_options(DistTessellator2D PRIVATE
  ${POLYTOPE_CXX_COMPILE_FLAGS})
target_link_options(DistTessellator2D PUBLIC
  ${POLYTOPE_CXX_LINK_FLAGS})
target_include_directories(DistTessellator2D SYSTEM PRIVATE
  ${POLYTOPE_ROOT_DIR}/tests
  ${POLYTOPE_ROOT_DIR}/src
  ${PROJECT_BINARY_DIR}/src
)

# Add the CXX benchmarks
blt_add_benchmark(NAME DistTessellator2DCXX
  COMMAND ${CMAKE_BINARY_DIR}/benchmarks/DistTessellator2D 1000000
  NUM_MPI_TASKS 36
)

# Add the Python benchmarks
blt_add_benchmark(NAME DistTessellator2DPY
  COMMAND ${CMAKE_BINARY_DIR}/${POLYTOPE_VIRT_DIR}/bin/python -u ${CMAKE_CURRENT_SOURCE_DIR}/DistTessellator2D.py 1000000
  NUM_MPI_TASKS 36
)
