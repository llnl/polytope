"""
Runs a serial tessellation multiple times.
Usage:
PolytopePy/bin/python3 -u SerialTessellator2D.py Nr Np
where Nr is the number of times to run each tessellator
and Np is the number of generator points.
By default Nr = 5 and Np = 50000
"""

import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope
import sys


def _available_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name) for name in names if hasattr(polytope, name)]


def _assert_mesh_populated(mesh):
    assert not mesh.empty()
    assert len(mesh.points) > 0
    assert len(mesh.nodes) > 0
    assert len(mesh.faces) > 0
    assert len(mesh.cells) > 0
    assert len(mesh.faceCells) > 0
    assert len(mesh.nodesAsTuples) == len(mesh.nodes)
    assert len(mesh.zoneNodes) == len(mesh.cells)


def test_serial_2d_tessellators(Ngen, nruns):
    tessellator_types = _available_tessellators()
    assert tessellator_types
    seed = 1049600
    boundary = Boundary2d(10)
    generator = ptu.PointGenerator(2, boundary2d=boundary)
    points = generator.generate_normal_random_points(Ngen, seed)

    for tessellator_type in tessellator_types:
        tessellator = tessellator_type()
        tess_name = tessellator.name()
        assert tess_name

        print(f"Testing unbounded {tess_name}")
        with ptu.timer(f"unbounded_{tess_name}"):
            for i in range(nruns):
                mesh = tessellator.tessellate(points)
        _assert_mesh_populated(mesh)
        ptu.outputMesh2d(mesh=mesh,
                         filePrefix=f"PySerial{tess_name}",
                         cycle=0,
                         time=0.,
                         numFiles=1)

        if (tessellator.clippingEnabled()):
            print(f"Testing clipped {tess_name}")
            with ptu.timer(f"clipped_{tess_name}"):
                for i in range(nruns):
                    clipped_mesh = tessellator.tessellate(points, boundary.PLCpoints, boundary.PLC)
            _assert_mesh_populated(clipped_mesh)
            ptu.outputMesh2d(mesh=clipped_mesh,
                             filePrefix=f"PySerial{tess_name}",
                             cycle=1,
                             time=1.,
                             numFiles=1)

if __name__ == "__main__":
    N = int(50000)
    nruns = 5
    if (len(sys.argv) > 1):
        nruns = int(sys.argv[1])
        if (len(sys.argv) > 2):
            N = int(sys.argv[2])
    print(f"Tessellating {N} points {nruns} times")
    test_serial_2d_tessellators(N, nruns)
