import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope
import sys, time

def test_partitioner(Ngen):
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    nranks = comm.getNRanks()
    # Seed for generating points randomly
    seed = 19001
    # Seed for partitioning points randomly
    partseed = 19
    boundary = Boundary2d(0)
    # Generate Ngen points per rank
    Ntotal = Ngen*nranks
    generator = ptu.PointGenerator(2, boundary2d=boundary)
    all_points = generator.generate_normal_random_points(Ntotal, seed=seed)
    serial_tessellator = polytope.BoostTessellator()
    tessellator = polytope.DistributedTessellator2d(serial_tessellator)
    # Make a partitioner
    partitioner = polytope.QuasiVoronoiPartitioner2d(partseed)
    ptu.rootprint(f"Tessellating {Ntotal} generators")
    for k, i in enumerate([0, 1, 10, 100]):
        partitioner.niter = i
        rank_points = partitioner.computeLocalPartition(all_points)
        with ptu.timer("Tessellate"):
            mesh = tessellator.tessellate(rank_points)
        ptu.outputMesh2d(mesh=mesh,
                         filePrefix="PartitionerTest",
                         cycle = k,
                         time = k)

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(10000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_partitioner(N)
