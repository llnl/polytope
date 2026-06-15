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
#polytope_add_new_test( "PLC_CSG" )  #  <-- Suspending 'til we think it works
polytope_add_new_test(HashKey2D)
polytope_add_new_test(HashKey3D)
polytope_add_new_test(Quantizer)
polytope_add_new_test(QuantPLC2D)
polytope_add_new_test(QuantPLC3D)
polytope_add_new_test(QuantConvexIntersect3D)
polytope_add_new_test(clipPolygon2D)
polytope_add_new_test(clipPolyhedron3D)
polytope_add_new_test(StarBoundary)
polytope_add_new_test(SolidRotationAroundHoles)
polytope_add_new_test(walkBoxEdges)
#polytope_add_new_test(Cylinder)

#-----------------------------------------------------------------------------
# Serial Tessellator Tests
#-----------------------------------------------------------------------------

# Voro++ Tessellator tests
#
# NOTE: Temporarily suspending their build until more work can be
#       done making their output topologically consistent
#
#polytope_add_new_test( "VoroPP_2d"                                )
#polytope_add_new_test( "VoroPP_3d"                                )

# 2D Triangle and Boost Tessellator tests
# if(POLYTOPE_ENABLE_TETGEN)
#   # 3D Tetgen Tessellator tests
#   polytope_add_new_test(TetgenTessellator)
# endif()
if(POLYTOPE_ENABLE_TRIANGLE)
  # polytope_add_new_test(plot)
  # polytope_add_new_test(AspectRatio)
  if(ENABLE_MPI)
    polytope_add_new_test(DistributedCollinear NUMTASKS 4)
    polytope_add_new_test(SerialDistributedTriangle NUMTASKS 4)
    polytope_add_new_test(FailedCommunication NUMTASKS 4)
  endif()
endif()

# Distributed Tessellator Tests
#-----------------------------------------------------------------------------
if(POLYTOPE_ENABLE_MPI)
  polytope_add_new_test(DistributedUnitSquare NUMTASKS 4)
  polytope_add_new_test(DistributedRandomPoints NUMTASKS 4)
  polytope_add_new_test(DistributedRotationTests NUMTASKS 4)
  polytope_add_new_test(DistributedCentroidal NUMTASKS 4)
  #polytope_add_new_test( DistributedVoroPP_2d 4)
  #polytope_add_new_test( DistributedVoroPP_3d 4)
endif()
