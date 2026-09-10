import sys, os
ptu_path = os.path.join(os.path.dirname(__file__), "../tests/PythonTests")
sys.path.append(ptu_path)
import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope

# Number of iterations to test with the partitioner
part_iters = [0, 10, 50, 100]

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
    boundary = Boundary2d(0)
    local_points = ptu.generate_normal_random_points(Ngen, seed=local_seed, boundary2d=boundary)
    # Make a partitioner
    partitioner = polytope.DistributedLloydPartitioner2d()
    if (rank == root):
        print(f"Tessellating {Ntotal} generators")
    serial_tessellator = get_tessellator()
    tess_name = serial_tessellator.name()
    if (rank == root):
        print(f"Testing unbounded {tess_name}")
    tessellator = polytope.DistributedTessellator2d(serial_tessellator)

    for k, i in enumerate([0, 10, 100]):
        partitioner.niter = i
        with ptu.timer(f"partition_{i}_iter"):
            points = partitioner.partition(local_points)
        mesh = polytope.Tessellation2d()
        with ptu.timer(f"tessellate_{i}_iter"):
            tessellator.tessellate(points, mesh)
        locfields = ptu.make_test_fields(mesh)
        polytope.writeSilo(mesh=mesh,
                           filePrefix="PyDistLloyd",
                           fields=locfields,
                           cycle=k,
                           time=float(k))

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(10000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_distributed_2d_tessellators(N)
