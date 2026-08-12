import polytope


def _vector_unsigned(values):
    result = polytope.vector_of_unsigned()
    for value in values:
        result.append(value)
    return result


def _vector_vector_unsigned(rows):
    result = polytope.vector_of_vector_of_unsigned()
    for row in rows:
        result.append(_vector_unsigned(row))
    return result


def _available_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name) for name in names if hasattr(polytope, name)]


def _square_plc():
    plc = polytope.PLC2d()
    plc.facets = _vector_vector_unsigned(((0, 1), (1, 2), (2, 3), (3, 0)))
    assert plc.valid()
    return plc


def _assert_mesh_populated(mesh):
    assert not mesh.empty()
    assert len(mesh.points) > 0
    assert len(mesh.nodes) > 0
    assert len(mesh.faces) > 0
    assert len(mesh.cells) > 0
    assert len(mesh.faceCells) > 0
    assert len(mesh.boundaryNodes) > 0
    assert len(mesh.boundaryFaces) > 0
    assert len(mesh.nodesAsTuples) == len(mesh.nodes)
    assert len(mesh.zoneNodes) == len(mesh.cells)


def test_serial_2d_tessellators():
    tessellator_types = _available_tessellators()
    assert tessellator_types

    points = [(0.25, 0.25), (0.75, 0.25), (0.5, 0.75)]
    plc_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]

    for tessellator_type in tessellator_types:
        tessellator = tessellator_type()
        assert tessellator.name()
        assert tessellator.handlesPLCs()

        mesh = polytope.Tessellation2d()
        tessellator.tessellate(points, mesh)
        _assert_mesh_populated(mesh)

        plc_mesh = polytope.Tessellation2d()
        tessellator.tessellate(points, plc_points, _square_plc(), plc_mesh)
        _assert_mesh_populated(plc_mesh)


def test_reduced_plc_constructor():
    plc_points = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0]
    reduced = polytope.ReducedPLC2d(_square_plc(), plc_points)
    assert reduced.valid()
    assert len(reduced.points) == len(plc_points)


if __name__ == "__main__":
    test_serial_2d_tessellators()
    test_reduced_plc_constructor()
