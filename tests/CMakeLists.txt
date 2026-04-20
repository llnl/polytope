#-----------------------------------------------------------------------------
# Adding a C++ test to Polytope
#
# The polytope test directory assumes all tests follow the naming convention
#      "<name>.cc"
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
# Polytope Utility Tests
#-----------------------------------------------------------------------------
polytope_add_test(serialize)
polytope_add_test(convexHull_2d)
polytope_add_test(convexHull_3d)
polytope_add_test(nearestPoint_2d)
polytope_add_test(within_2d)
polytope_add_test(intersect_2d)
polytope_add_test(convexIntersect)
#polytope_add_test( "PLC_CSG" )  #  <-- Suspending 'til we think it works
polytope_add_test(hash)

#-----------------------------------------------------------------------------
# Serial Tessellator Tests
#-----------------------------------------------------------------------------

# Voro++ Tessellator tests
#
# NOTE: Temporarily suspending their build until more work can be
#       done making their output toplogically consistent
#
#polytope_add_test( "VoroPP_2d"                                )
#polytope_add_test( "VoroPP_3d"                                )

# 2D Triangle and Boost Tessellator tests
polytope_add_test(UnitSquare)
polytope_add_test(MeshEditor)
polytope_add_test(Degenerate)
polytope_add_test(BoundaryJitter)
polytope_add_test(RandomPoints)
polytope_add_test(Area)
polytope_add_test(OrphanedCell)
polytope_add_test(TwoGenerators)
polytope_add_test(Centroidal)
polytope_add_test(UnboundedToBounded)
polytope_add_test(Unbounded)
polytope_add_test(OrphanCases)
polytope_add_test(StarBoundary)
polytope_add_test(Circumcenters)
polytope_add_test(RotationTests)
polytope_add_test(SolidRotationAroundHoles)
polytope_add_test(TiltedLattice)
polytope_add_test(ProjectionIntersection)
if(POLYTOPE_ENABLE_BOOST)
  polytope_add_test(BoostTessellator)
  polytope_add_test(BoostGeometry)
endif()
if(POLYTOPE_ENABLE_TETGEN)
  # 3D Tetgen Tessellator tests
  polytope_add_test(testTetgenTessellator)
endif()
if(POLYTOPE_ENABLE_TRIANGLE)
  # polytope_add_test(plot)
  # polytope_add_test(AspectRatio)
  if(ENABLE_MPI)
    polytope_add_test(DistributedCollinear NUMTASKS 4)
    polytope_add_test(SerialDistributedTriangle NUMTASKS 4)
    polytope_add_test(FailedCommunication NUMTASKS 4)
  endif()
endif()

# Distributed Tessellator Tests
#-----------------------------------------------------------------------------
if(POLYTOPE_ENABLE_MPI)
  polytope_add_test(DistributedUnitSquare NUMTASKS 4)
  polytope_add_test(DistributedRandomPoints NUMTASKS 4)
  polytope_add_test(DistributedRotationTests NUMTASKS 4)
  polytope_add_test(DistributedCentroidal NUMTASKS 4)
  #polytope_add_test( DistributedVoroPP_2d 4)
  #polytope_add_test( DistributedVoroPP_3d 4)
endif()
