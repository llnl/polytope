import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope
import sys


def _available_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name) for name in names if hasattr(polytope, name)]


def test_distributed_2d_tessellators(Ngen):
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNRanks()
    tessellator_types = _available_tessellators()
    #assert tessellator_types
    Q = polytope.Quantizer2d.instance()

    # Seed for generating points randomly
    seed = 19001
    # Seed for partitioning points randomly
    partseed = seed + 10
    boundary = Boundary2d(10)
    # Generate Ngen points per rank
    Ntotal = Ngen*nranks
    all_points = ptu.generate_random_points(Ntotal, seed=seed, boundary2d=boundary)
    # Make a partitioner
    partitioner = polytope.QuasiVoronoiPartitioner2d(partseed)
    if (rank == root):
        print(f"Tessellating {Ntotal} generators")
    points = partitioner.computeLocalPartition(all_points)

    for tessellator_type in tessellator_types:
        serial_tessellator = tessellator_type()
        tess_name = serial_tessellator.name()
        assert tess_name

        if (rank == root):
            print(f"Testing unbounded {tess_name}")
        mesh = polytope.Tessellation2d()
        tessellator = polytope.DistributedTessellator2d(serial_tessellator)
        with ptu.timer("Unbounded tessellate"):
            tessellator.tessellate(points, mesh)
        locfields = ptu.make_test_fields(mesh)
        polytope.writeSilo(mesh=mesh,
                           filePrefix=f"PyDist{tess_name}",
                           fields=locfields,
                           cycle=0,
                           time=0.)

        if (rank == root):
            print(f"Testing clipped {tess_name}")
        clipped_mesh = polytope.Tessellation2d()
        with ptu.timer("Clipped tessellation"):
            tessellator.tessellate(points, boundary.PLCpoints, boundary.PLC, clipped_mesh)
        locfields = ptu.make_test_fields(clipped_mesh)
        polytope.writeSilo(mesh=clipped_mesh,
                           filePrefix=f"PyDist{tess_name}",
                           fields=locfields,
                           cycle=1,
                           time=1.)

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(10000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_distributed_2d_tessellators(N)
