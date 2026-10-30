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
    nranks = comm.getNRanks()
    tessellator_types = _available_tessellators()
    # Seed for generating points randomly
    seed = 19001
    # Seed for partitioning points randomly
    partseed = seed + 10
    boundary = Boundary2d(10)
    # Generate Ngen points per rank
    Ntotal = Ngen*nranks
    generator = ptu.PointGenerator(2, boundary2d=boundary)
    all_points = generator.generate_random_points(Ntotal, seed=seed)
    # Make a partitioner
    partitioner = polytope.QuasiVoronoiPartitioner2d(partseed)
    ptu.rootprint(f"Tessellating {Ntotal} generators")
    points = partitioner.computeLocalPartition(all_points)

    for tessellator_type in tessellator_types:
        serial_tessellator = tessellator_type()
        tess_name = serial_tessellator.name()
        assert tess_name

        ptu.rootprint(f"Testing unbounded {tess_name}")
        tessellator = polytope.DistributedTessellator2d(serial_tessellator)
        with ptu.timer("Unbounded tessellate"):
            mesh = tessellator.tessellate(points)
        ptu.outputMesh2d(mesh=mesh, 
                         filePrefix=f"PyDist{tess_name}",
                         cycle=0,
                         time=0.)

        if (tessellator.clippingEnabled()):
            ptu.rootprint(f"Testing clipped {tess_name}")
            with ptu.timer("Clipped tessellation"):
                clipped_mesh = tessellator.tessellate(points, boundary.PLCpoints, boundary.PLC)
            ptu.outputMesh2d(mesh=clipped_mesh,
                             filePrefix=f"PyDist{tess_name}",
                             cycle=1,
                             time=1.)

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(10000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_distributed_2d_tessellators(N)
