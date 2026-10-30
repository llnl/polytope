import sys, os, math, time
ptu_path = os.path.join(os.path.dirname(__file__), "../tests/PythonTests")
sys.path.append(ptu_path)
import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope

timer_log = "perf_benchmark_2d.log"

def get_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name)() for name in names if hasattr(polytope, name)]

def test_distributed_2d_tessellators(serial_tessellator, local_points):
    comm = polytope.Communicator.instance()
    tess_name = serial_tessellator.name()
    ptu.rootprint(f"Testing unbounded {tess_name}")
    tessellator = polytope.DistributedTessellator2d(serial_tessellator)
    comm.Barrier()
    startt = time.perf_counter()
    mesh = tessellator.tessellate(points)
    comm.Barrier()
    return time.perf_counter() - startt

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(100000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNRanks()
    coord_bits = round(math.log2(polytope.QuantizedKeyTraits2d().maxCoordinate()))
    if (rank == root):
        with open(timer_log, "w") as ff:
            ff.write(f"{nranks} Ranks, Coord bits {coord_bits}, {N} points per rank\n")
    boundary = Boundary2d(0)
    generator = ptu.PointGenerator(2, boundary2d=boundary)
    Ntotal = N*nranks
    local_seed = 19001 + rank
    # Test a uniform distribution of generators
    dist_name = "Uniform"
    local_points = generator.generate_random_points(N, seed=local_seed)
    partitioner = polytope.DistributedLloydPartitioner2d()
    points = partitioner.partition(local_points)
    if (rank == root):
        with open(timer_log, "a") as ff:
            ff.write(f"{dist_name} distribution\n")
    for tess in get_tessellators():
        tess_time = test_distributed_2d_tessellators(tess, points)
        if (rank == root):
            with open(timer_log, "a") as ff:
                ff.write(f"{tess.name()} time: {tess_time}\n")
    # Test a normal distribution of generators
    dist_name = "Normal"
    local_seed += 100
    local_points = generator.generate_normal_random_points(N, seed=local_seed)
    points = partitioner.partition(local_points)
    if (rank == root):
        with open(timer_log, "a") as ff:
            ff.write(f"{dist_name} distribution\n")
    for tess in get_tessellators():
        tess_time = test_distributed_2d_tessellators(tess, points)
        if (rank == root):
            with open(timer_log, "a") as ff:
                ff.write(f"{tess.name()} time: {tess_time}\n")
    ptu.rootprint(f"Times outputted to {timer_log}")
