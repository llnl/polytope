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


def test_serial_2d_tessellators(Ngen):
    tessellator_types = _available_tessellators()
    assert tessellator_types
    Q = polytope.Quantizer2d.instance()

    seed = 19001
    boundary = Boundary2d(10)
    points = ptu.generate_random_points(Ngen, seed, boundary)

    for tessellator_type in tessellator_types:
        tessellator = tessellator_type()
        tess_name = tessellator.name()
        assert tess_name

        print(f"Testing unbounded {tess_name}")
        mesh = polytope.Tessellation2d()
        tessellator.tessellate(points, mesh)
        _assert_mesh_populated(mesh)
        ptu.outputMesh2d(mesh=mesh,
                         filePrefix=f"PySerial{tess_name}",
                         cycle=0,
                         time=0.,
                         numFiles=1)

        print(f"Testing clipped {tess_name}")
        clipped_mesh = polytope.Tessellation2d()
        tessellator.tessellate(points, boundary.PLCpoints, boundary.PLC, clipped_mesh)
        _assert_mesh_populated(clipped_mesh)
        ptu.outputMesh2d(mesh=clipped_mesh,
                         filePrefix=f"PySerial{tess_name}",
                         cycle=1,
                         time=1.,
                         numFiles=1)

if __name__ == "__main__":
    N = int(50000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_serial_2d_tessellators(N)
