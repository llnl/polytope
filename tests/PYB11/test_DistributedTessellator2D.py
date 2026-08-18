import polytope_test_utilities as ptu
import polytope


def _available_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name) for name in names if hasattr(polytope, name)]


def test_parallel_2d_tessellators():
    tessellator_types = _available_tessellators()
    assert tessellator_types
    Q = polytope.Quantizer2d.instance()
    # Number of points per rank
    N = 10
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    oseed = 1049600
    seed = oseed + rank

    plc_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]

    for tessellator_type in tessellator_types:
        tessellator = tessellator_type()
        tess_name = tessellator.name()
        assert tess_name

        mesh = polytope.Tessellation2d()
        Q.initPoints(plc_points)
        points = ptu.generate_random_points(N, seed)
        ptess = polytope.DistributedTessellator2d(tessellator)
        ptess.tessellate(points, mesh)
        locfields = ptu.make_test_fields(mesh)
        polytope.writeSilo(mesh=mesh,
                           filePrefix=f"PyParallel{tess_name}",
                           fields=locfields,
                           cycle=0,
                           time=0.)

        plc_mesh = polytope.Tessellation2d()
        plc = polytope.PLC2d()
        plc.facets = ptu.make_square_facets()
        ptess.tessellate(points, plc_points, plc, plc_mesh)
        locfields = ptu.make_test_fields(plc_mesh)
        polytope.writeSilo(mesh=plc_mesh,
                           filePrefix=f"PyParallel{tess_name}",
                           fields=locfields,
                           cycle=1,
                           time=1.)


if __name__ == "__main__":
    test_parallel_2d_tessellators()
