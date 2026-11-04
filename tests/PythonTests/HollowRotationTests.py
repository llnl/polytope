"""
This creates a shrinking column of generators. Good test for unusual degeneracies.
"""

import sys
import polytope
from Boundary2D import Boundary2D
import polytope_test_utilities as ptu

def runTest(tessellator, ny, flowType, io=False):
    """
    This creates a top row of points across the domain.
    The next row has nx points going right to left.
    Each successive nrowc rows reduces ndx by dndx.
    """
    nx = ny - 2
    nrowc = 6
    dndx = 1
    ly = 0.9
    lx = ly
    dy = ly/ny
    dx = dy/2.
    y = ly
    points = []
    # Make the top row encompass almost the whole domain
    x0 = lx - dx
    x = x0
    while (x > dx):
        points.append([x, y])
        x -= dx
    y -= dy
    # Form column that shrinks every nrowc by dndx
    for j in range(ny):
        x = x0
        for i in range(nx):
            points.append([x, y])
            x -= dx
        y -= dy
        if (j%nrowc==0):
            nx -= dndx
    boundary = Boundary2D(0)
    boundary.mCenter = [0.5, 0.5]
    boundary.mDiff = ly
    boundary.setUnitSquare()
    blo, bhi = boundary.low, boundary.high

    testName = f"PyHollowRotationTests_{tessellator.name()}_{flowType}"
    description, dt, maxTime, scaleFactor = ptu.setConfig(flowType, bhi, blo, dx)
    ptu.rootprint(description)
    points = [0.5 + (x - 0.5)*scaleFactor for coordinate in points for x in coordinate]
    tessellator.setQBounds([0.0, 0.0, 1.0, 1.0])

    step, time = 0, 0.0
    while time < maxTime:
        if step % 5 == 0:
            ptu.rootprint(f"{round(100.0*time/maxTime)}%")
        mesh = tessellator.tessellate(points)
        if (step % 5 == 0 and io) or step == 0:
            ptu.outputMesh2d(mesh, testName, cycle=step, time=time)
        ptu.compareArea(mesh, tessellator.quantizer(), f"Test {tessellator.name()} step {step} {flowType}")
        velocities = ptu.getVelocities(points, bhi, blo, flowType)
        halfTimePositions = [point + dt*velocity
                             for point, velocity in zip(points, velocities)]
        velocities = ptu.getVelocities(halfTimePositions, bhi, blo, flowType)
        points = [point + dt*velocity for point, velocity in zip(points, velocities)]
        time += dt
        step += 1

def test_rotation(nx):
    if nx <= 0:
        raise ValueError("Number of points specified must be > 0")
    ttypes = ("TriangleTessellator", "BoostTessellator")
    # These are already tested in the C++ test
    tessellators = [getattr(polytope, name)() for name in ttypes if hasattr(polytope, name)]
    # Scipy tessellators only work with lower bit accuracy
    if (polytope.QuantizedKeyTraits2d().coordinateBits < 50):
        try:
            import ScipyTessellator2D
            tessellators.append(ScipyTessellator2D.ScipyVoronoi())
            tessellators.append(ScipyTessellator2D.ScipyDelaunay())
        except ModuleNotFoundError:
            print("No Scipy module found. Scipy tessellators not tested")
            pass
    for tessellator in tessellators:
        ptu.rootprint(f"\n{tessellator.name()}:\n")
        for flowType in range(1, 5):
            runTest(tessellator, nx, flowType)

if __name__ == "__main__":
    test_rotation(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
