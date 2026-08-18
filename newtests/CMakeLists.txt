#-----------------------------------------------------------------------------
# Adding a C++ test to Polytope
#
# The polytope test directory assumes all tests follow the naming convention
#      "test_<name>.cc"
#
# To add a test, the user must register the test <name> and a list of the
# dependencies needed to build it:
#      polytope_add_test(name numtasks)

# where numtasks defaults to 1 if no argument is provided.
#
#-----------------------------------------------------------------------------

# The polytope_add_test macro
include(PolytopeAddTests)

#-----------------------------------------------------------------------------
# Polytope Setup/Cleanup
#-----------------------------------------------------------------------------

# Directory to run tests in to make cleanup easier
set(NEW_TEST_WORK_DIR ${CMAKE_BINARY_DIR}/newtests/test_silo_files)

# Setup tests
blt_add_test(NAME new_setup_tests
  COMMAND ${CMAKE_COMMAND} -E make_directory "${NEW_TEST_WORK_DIR}"
)

set_tests_properties(new_setup_tests PROPERTIES FIXTURES_SETUP "new_polytope_fixtures")

# Cleanup test files
blt_add_test(NAME new_cleanup_tests
 COMMAND ${CMAKE_COMMAND} -E rm -rf "${NEW_TEST_WORK_DIR}"
)

set_tests_properties(new_cleanup_tests PROPERTIES FIXTURES_CLEANUP "new_polytope_fixture")

#-----------------------------------------------------------------------------
# Polytope Utility Tests
#-----------------------------------------------------------------------------
# Cannot test any clipped geometries when hi bit 2D is enabled
if(NOT POLYTOPE_ENABLE_HIBIT2D)
  polytope_add_new_test(HashKey2D)
  polytope_add_new_test(SimpleTess2D)
  polytope_add_new_test(StarBoundary)
  polytope_add_new_test(SolidRotationAroundHoles)
  polytope_add_new_test(OrphanCases)
  polytope_add_new_test(BoostTessellator)
  polytope_add_new_test(RandomPoints)
  polytope_add_new_test(RotationTests)
  polytope_add_new_test(Circumcenter)
endif()
polytope_add_new_test(HashKey3D)
polytope_add_new_test(Quantizer)
polytope_add_new_test(QuantPLC2D)
#polytope_add_new_test(QuantPLC3D) # Currently broken
polytope_add_new_test(QuantConvexIntersect3D)
#polytope_add_new_test(clipPolygon2D) # Currently broken, only used by 3D
polytope_add_new_test(clipPolyhedron3D)
polytope_add_new_test(honeyComb)

#-----------------------------------------------------------------------------
# Tetgen Tessellator Tests
#-----------------------------------------------------------------------------

if(POLYTOPE_ENABLE_TETGEN)
  # 3D Tetgen Tessellator tests
  polytope_add_new_test(TetgenTessellator)
endif()

#-----------------------------------------------------------------------------
# Distributed Tessellator Tests
#-----------------------------------------------------------------------------
if(POLYTOPE_ENABLE_MPI AND NOT POLYTOPE_ENABLE_HIBIT2D)
  if(NOT POLYTOPE_ENABLE_HIBIT2D)
    polytope_add_new_test(DistributedRandomPoints NUMTASKS 4)
    polytope_add_new_test(DistributedOrderedPoints NUMTASKS 4)
  endif()
  polytope_add_new_test(DistributedIO2D NUMTASKS 4)
  polytope_add_new_test(DistributedVoronoi2D NUMTASKS 4)
  polytope_add_new_test(DistributedVisible2D NUMTASKS 6)
  polytope_add_new_test(DistributedRotationTests NUMTASKS 4)
endif()

#-----------------------------------------------------------------------------
# Python Binding Tests
#-----------------------------------------------------------------------------
if(POLYTOPE_ENABLE_PYTHON AND NOT POLYTOPE_ENABLE_LEGACY_PYTHON)
  polytope_add_python_test(Quantizer)
  polytope_add_python_test(HashKey)
  polytope_add_python_test(SerialTessellator2D)
  polytope_add_python_test(PLC)
  if(POLYTOPE_ENABLE_MPI)
    polytope_add_python_test(DistributedTessellator2D NUMTASKS 4)
  endif()
endif()
