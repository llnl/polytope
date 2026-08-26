import sys, time, os
ptu_path = os.path.join(os.path.dirname(__file__), "../tests/PythonTests")
sys.path.append(ptu_path)
import polytope_test_utilities as ptu
import polytope


def _available_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name) for name in names if hasattr(polytope, name)]


def test_parallel_2d_tessellators(N):
    tessellator_types = _available_tessellators()
    assert tessellator_types
    Q = polytope.Quantizer2d.instance()
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNProcs()
    seed = 1049600

    plc_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    Q.init(plc_points)
    if (rank == root):
            print(f"Degeneracy {Q.degeneracy()}")
            print(f"Generating {N} random points")
    gen_begin = time.perf_counter()
    points = ptu.generate_random_points(N, seed)
    gen_time = time.perf_counter() - gen_begin

    for tessellator_type in tessellator_types:
        tessellator = tessellator_type()
        tess_name = tessellator.name()
        ptess = polytope.DistributedTessellator2d(tessellator)
        # Determine how to distribute up the points
        r0 = int(nranks**0.5)
        r1 = int(nranks/r0)
        if (rank == root):
            print(f"Doing {tess_name} tessellation")
        parts = [polytope.LatticePartitioner2d([r0, r1])]
        for _ in range(2):
            for partitioner in parts:
                mesh = polytope.Tessellation2d()
                comm.Barrier()
                tess_begin = time.perf_counter()
                ptess.partitionAndTessellate(points, partitioner, mesh)
                comm.Barrier()
                tess_time = time.perf_counter() - tess_begin
                if (rank == root):
                    print(f"exchange points: {ptess.exchangePoints()}, {partitioner.name()}: {tess_time}")
            ptess.setExchangePoints(True)

if __name__ == "__main__":
    N = int(1E6)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_parallel_2d_tessellators(N)
