import sys
import polytope_test_utilities as ptu
from Boundary2D import Boundary2D
import polytope

# Number of iterations to test with the partitioner
part_iters = [10, 50, 100]

def get_tessellator():
    tessname = "TriangleTessellator"
    if (hasattr(polytope, tessname)):
        return getattr(polytope, tessname)()
    return getattr(polytope, "BoostTessellator")()

def test_distributed_2d_tessellators(Ngen):
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNRanks()
    Ntotal = Ngen*nranks

    # Seed for generating points in a normal distribution
    local_seed = 19001 + rank
    boundary = Boundary2D(0)
    generator = ptu.PointGenerator(2, boundary2d=boundary)
    local_points = generator.generate_normal_random_points(Ngen, seed=local_seed)
    # Make a partitioner
    partitioner = polytope.DistributedLloydPartitioner2d()
    ptu.rootprint(f"Tessellating {Ntotal} generators with {nranks} ranks")
    serial_tessellator = get_tessellator()
    tess_name = serial_tessellator.name()
    ptu.rootprint(f"Testing unbounded {tess_name}")
    tessellator = polytope.DistributedTessellator2d(serial_tessellator)

    for k, i in enumerate(part_iters):
        partitioner.niter = i
        with ptu.timer(f"partition_{i}_iter"):
            points = partitioner.partition(local_points)
        with ptu.timer(f"tessellate_{i}_iter"):
            mesh = tessellator.tessellate(points)
        ptu.outputMesh2d(mesh=mesh,
                         filePrefix="PyDistLloyd",
                         cycle=k,
                         time=float(k))

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(100000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_distributed_2d_tessellators(N)
