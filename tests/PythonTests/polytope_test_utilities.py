import polytope
import random, builtins
import math
import sys
from time import perf_counter
from contextlib import contextmanager
import Boundary2d

# Print statement wrapper for printing only on the root rank
def rootprint(string):
    comm = polytope.Communicator.instance()
    if (comm.getRank() == comm.getRoot()):
        print(string)

def outputMesh2d(mesh,
                 filePrefix,
                 fields=None,
                 cycle=0,
                 time=0.,
                 numFiles=-1):
    writer = polytope.SiloWriter2d(mesh)
    writer.generateTestVars()
    if (fields):
        for key, vals in fields.item():
            writer.addField(polytope.FieldCentering.Cell,
                            key,
                            vals)
    writer.write(filePrefix, cycle, time, numFiles)

class PointGenerator:
    def __init__(self, dim = 2,
                 boundary2d = None,
                 lobounds = None,
                 hibounds = None):
        self.dim = dim
        self.boundary2d = boundary2d
        xmin = None
        xmax = None
        if boundary2d:
            lobounds = boundary2d.low
            hibounds = boundary2d.high
        if lobounds and hibounds:
            if dim == 3:
                xmin = polytope.Point3d(*lobounds)
                xmax = polytope.Point3d(*hibounds)
            else:
                xmin = polytope.Point2d(*lobounds)
                xmax = polytope.Point2d(*hibounds)
        if (not xmax or not xmin):
            raise ValueError("Must specify lobounds/hibounds, or a boundary2d")
        self.xmin = xmin
        self.xmax = xmax

    def generate_random_points(self, N, seed = -1):
        if (seed >= 0):
            random.seed(seed)
        L = (self.xmax - self.xmin)
        pout = []
        if (self.boundary2d and self.boundary2d.mType > 1):
            while(len(pout) < self.dim*N):
                point = [self.xmin[d] + random.random()*L[d] for d in range(self.dim)]
                if (self.boundary2d.testInside(point)):
                    pout.extend(point)
        else:
            for _ in range(N):
                pout.extend([self.xmin[d] + random.random()*L[d] for d in range(self.dim)])
        return pout

    def generate_normal_random_points(self, N, seed = -1,
                                      mu = 0.5, sigma = 1./6.17737):
        if (seed >= 0):
            random.seed(seed)
        L = (self.xmax - self.xmin)
        pout = []
        if (self.boundary2d and self.boundary2d.mType > 1):
            while(len(pout) < self.dim*N):
                point = [self.xmin[d] + random.normalvariate(mu, sigma)*L[d] for d in range(self.dim)]
                if (self.boundary2d.testInside(point)):
                    pout.extend(point)
        else:
            while(len(pout) < self.dim*N):
                rnums = [random.normalvariate(mu, sigma) for _ in range(self.dim)]
                if (min(rnums) >= 0. and max(rnums) <= 1.):
                    for d, rnum in enumerate(rnums):
                        pout.append(self.xmin[d] + rnum*L[d])
        return pout

    def generate_cartesian_points(self, Nx):
        L = (self.xmax - self.xmin)
        pout = []
        dx = [L[d]/Nx for d in range(self.dim)]
        for i in range(Nx):
            x = self.xmin[0] + (i + 0.5)*dx[0]
            for j in range(Nx):
                y = self.xmin[1] + (j + 0.5)*dx[1]
                if (self.dim == 3):
                    for k in range(Nx):
                        z = self.xmin[2] + (k + 0.5)*dx[2]
                        pout.extend([x, y, z])
                else:
                    pout.extend([x, y])
        return pout

@contextmanager
def timer(name):
    """
    Time a code section. Run as:

    with timer("my_functions"):
        function_one()
        function_two()
    """
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    comm.Barrier()
    tstart = perf_counter()

    try:
        yield
    finally:
        comm.Barrier()
        elapsed = perf_counter() - tstart
        rootprint(f"{name}: {elapsed:.6f} seconds")


# Helpers shared by the rotational tessellation tests. Coordinates are flat
# [x0, y0, x1, y1, ...] lists, as accepted by tessellate().
def getVelocities(points, bhi, blo, flowType):
    if len(points) % 2:
        raise ValueError("2-D point coordinates must come in pairs")
    if flowType not in (1, 2, 3, 4):
        raise ValueError("Unknown rotation flow type")

    length = (bhi[0] - blo[0], bhi[1] - blo[1])
    center = ((bhi[0] + blo[0])*0.5, (bhi[1] + blo[1])*0.5)
    maxRadius = math.hypot(*length)*0.5
    omega = 1.0/maxRadius
    velocities = []
    for x, y in zip(points[::2], points[1::2]):
        ux = (x - blo[0])/length[0]
        uy = (y - blo[1])/length[1]
        if flowType == 1:
            vx = -omega*(y - center[1])
            vy = omega*(x - center[0])
        elif flowType == 2:
            vx = length[0]*math.sin(math.pi*ux)*math.cos(math.pi*uy)
            vy = -length[1]*math.cos(math.pi*ux)*math.sin(math.pi*uy)
        elif flowType == 3:
            vx = 0.5*length[0]*math.cos(2.0*math.pi*(ux - 0.25))*math.sin(2.0*math.pi*(uy - 0.25))
            vy = -0.5*length[1]*math.sin(2.0*math.pi*(ux - 0.25))*math.cos(2.0*math.pi*(uy - 0.25))
        else:
            vx = -length[0]*math.sin(4.0*math.pi*(ux + 0.5))*math.sin(4.0*math.pi*(uy + 0.5))
            vy = -length[1]*math.cos(4.0*math.pi*(ux + 0.5))*math.cos(4.0*math.pi*(uy + 0.5))
        velocities.extend((vx, vy))
    return velocities


def setConfig(flowType, bhi, blo, dx):
    cfl = math.sqrt(2.0)
    length = (bhi[0] - blo[0], bhi[1] - blo[1])
    maxRadius = math.hypot(*length)*0.5
    if flowType == 1:
        return "Test 1: Solid Rotation Flow", cfl*dx, math.pi*maxRadius, 1.0/math.sqrt(2.0)
    if flowType == 2:
        return "Test 2: Single Vortex Flow", cfl*dx, 2.0, 1.0
    if flowType == 3:
        return "Test 3: Taylor-Green (4-Vortex) Flow", cfl*dx, 2.0, 1.0
    if flowType == 4:
        return "Test 4: Deformation (16-Vortex) Flow", 0.5*dx, 2.0, 1.0/math.sqrt(2.0)
    raise ValueError("Unknown rotation flow type")


def _signedCellArea(nodeIds, nodes):
    vertices = [nodes[i] for i in nodeIds]
    return 0.5*math.fsum(a[0]*b[1] - a[1]*b[0]
                         for a, b in zip(vertices, vertices[1:] + vertices[:1]))


def computeTessellationArea(mesh):
    nodes = mesh.nodesAsTuples
    return math.fsum(_signedCellArea(cell, nodes) for cell in mesh.zoneNodes)


def errorTolerance(quantizer, nCells):
    eps = sys.float_info.epsilon
    nAdds = float(nCells - 1)
    gamma = (nAdds*0.5*eps)/(1.0 - nAdds*0.5*eps)
    degeneracy = quantizer.degeneracy()
    return 2.0*(degeneracy.x + degeneracy.y) + 4.0*gamma + 32.0*eps


def compareArea(mesh, quantizer, testString=""):
    cells = mesh.zoneNodes
    nodes = mesh.nodesAsTuples
    areas = [_signedCellArea(cell, nodes) for cell in cells]
    area = math.fsum(areas)
    expectedArea = quantizer.area
    relErr = abs(expectedArea - area)/expectedArea
    tolerance = errorTolerance(quantizer, len(mesh.points))
    assert relErr < tolerance, (f"Error in area: error = {relErr} "
                                f"tolerance {tolerance} {testString}")
    for i, cellArea in enumerate(areas):
        assert cellArea > 0.0, f"Bad orientation at cell {i} {testString}"
