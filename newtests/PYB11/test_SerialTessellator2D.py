from polytope_test_utilities import *
import polytope


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


def test_serial_2d_tessellators():
    tessellator_types = _available_tessellators()
    assert tessellator_types
    Q = polytope.Quantizer2d.instance()

    points = [(0.25, 0.25), (0.75, 0.25), (0.5, 0.75)]
    plc_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]

    for tessellator_type in tessellator_types:
        tessellator = tessellator_type()
        assert tessellator.name()

        mesh = polytope.Tessellation2d()
        Q.initPoints(points)
        tessellator.tessellate(points, mesh)
        _assert_mesh_populated(mesh)

        plc_mesh = polytope.Tessellation2d()
        plc = polytope.PLC2d()
        plc.facets = make_square_facets()
        Q.initPoints(plc_points)
        tessellator.tessellate(points, plc_points, plc, plc_mesh)
        _assert_mesh_populated(plc_mesh)


if __name__ == "__main__":
    test_serial_2d_tessellators()
