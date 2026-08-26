# Performance benchmarks are intentionally separate from CTest correctness
# tests and are not installed.

blt_add_executable(NAME DistributedTessellationBenchmark
  SOURCES DistributedTessellationBenchmark.cc
  DEPENDS_ON polytopeC ${POLYTOPE_TPL_DEPENDS}
  OUTPUT_DIR ${CMAKE_BINARY_DIR}/benchmarks
)

target_compile_options(DistributedTessellationBenchmark PRIVATE
  ${POLYTOPE_CXX_COMPILE_FLAGS})
target_link_options(DistributedTessellationBenchmark PUBLIC
  ${POLYTOPE_CXX_LINK_FLAGS})
target_include_directories(DistributedTessellationBenchmark SYSTEM PRIVATE
  ${POLYTOPE_ROOT_DIR}/tests
  ${POLYTOPE_ROOT_DIR}/src
  ${PROJECT_BINARY_DIR}/src
)
