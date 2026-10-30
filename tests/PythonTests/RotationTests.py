"""
Exercise the four rotation fields on a Cartesian 2-D tessellation.
This test will also use Scipy if it is installed
"""

import sys

import ScipyTessellator2D
import polytope
from Boundary2d import Boundary2d
import polytope_test_utilities as ptu


def runTest(tessellator, flowType, nx):
    boundary = Boundary2d(0)
    boundary.mCenter = [0.5, 0.5]
    boundary.mDiff = 0.4
    boundary.setUnitSquare()
    blo, bhi = boundary.low, boundary.high
    dx = (bhi[0] - blo[0])/nx
    points = ptu.PointGenerator(2, boundary2d=boundary).generate_cartesian_points(nx)

    testName = f"PyRotationTests_{tessellator.name()}_{flowType}"
    description, dt, maxTime, scaleFactor = ptu.setConfig(flowType, bhi, blo, dx)
    ptu.rootprint(description)
    points = [0.5 + (coordinate - 0.5)*scaleFactor for coordinate in points]
    tessellator.setQBounds([0.0, 0.0, 1.0, 1.0])

    step, time = 0, 0.0

    while time < maxTime:
        if step % 5 == 0:
            ptu.rootprint(f"{round(100.0*time/maxTime)}%")
        mesh = tessellator.tessellate(points)
        ptu.outputMesh2d(mesh, testName, cycle=step, time=time)
        ptu.compareArea(mesh, tessellator.quantizer(), f"Test {flowType}")
        velocities = ptu.getVelocities(points, bhi, blo, flowType)
        halfTimePositions = [point + dt*velocity
                             for point, velocity in zip(points, velocities)]
        velocities = ptu.getVelocities(halfTimePositions, bhi, blo, flowType)
        points = [point + dt*velocity for point, velocity in zip(points, velocities)]
        time += dt
        step += 1

def test_rotation(nx=5):
    if nx <= 0:
        raise ValueError("Number of points specified must be > 0")
    ttypes = ("TriangleTessellator", "BoostTessellator")
    # These are already tested in the C++ test
    tessellators = [getattr(polytope, name)() for name in ttypes if hasattr(polytope, name)]
    # If Scipy is installed, use the wrappers provided by polytope
    try:
        import scipy
        tessellators.append(ScipyTessellator2D.ScipyVoronoi())
        tessellators.append(ScipyTessellator2D.ScipyDelaunay())
    except ModuleNotFoundError:
        print("No Scipy module found. Scipy tessellators not tested")
        pass
    for tessellator in tessellators:
        ptu.rootprint(f"\n{tessellator.name()}:\n")
        for flowType in range(1, 5):
            runTest(tessellator, flowType, nx)

if __name__ == "__main__":
    test_rotation(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
