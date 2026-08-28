import sys, time, os
ptu_path = os.path.join(os.path.dirname(__file__), "../tests/PythonTests")
sys.path.append(ptu_path)
import polytope_test_utilities as ptu
import polytope

def _available_tessellators():
    names = ("BoostTessellator", "TriangleTessellator")
    return [getattr(polytope, name) for name in names if hasattr(polytope, name)]

def time_tessellation(allpoints, tessellator):
    Q = polytope.Quantizer2d.instance()
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNProcs()
    seed = 1049600
    partseed = 1042390
    ptess = polytope.DistributedTessellator2d(tessellator)
    # Determine how to load balance the generators for Lattice partitioner
    r0 = int(nranks**0.5)
    r1 = int(nranks/r0)
    if (rank == root):
        print(f"Using {tessellator.name()}")
    parts = [polytope.LatticePartitioner2d([r0, r1]),
             polytope.QuasiVoronoiPartitioner2d(partseed)]
    time_dicts = []
    for exchangetype in range(2):
        if exchangetype == 0:
            ptess.setExchangePoints(False)
        else:
            ptess.setExchangePoints(True)
        for encodingtype in range(2):
            if encodingtype == 0:
                Q.useMortonEncoding()
            else:
                Q.usePackedEncoding()
            for partitioner in parts:
                mesh = polytope.Tessellation2d()
                comm.Barrier()
                tess_begin = time.perf_counter()
                ptess.partitionAndTessellate(points, partitioner, mesh)
                comm.Barrier()
                tess_time = time.perf_counter() - tess_begin
                time_dict = {"MPI exchange": "points" if ptess.exchangePoints() else "hashes",
                             "partition": partitioner.name(),
                             "encoding": Q.keyName(),
                             "time": tess_time}
                time_dicts.append(time_dict)
    return time_dicts

if __name__ == "__main__":
    N = int(1E6)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
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

    for disttype in range(2):
        gen_begin = time.perf_counter()
        if (disttype == 0):
            points = ptu.generate_random_points(N, seed)
            distname = "uniform"
        else:
            points = ptu.generate_normal_random_points(N, seed=seed)
            distname = "normal"
        if (rank == root):
            print(f"Generators distributed in a {distname} distribution")
        gen_time = time.perf_counter() - gen_begin
        if (hasattr(polytope, "TriangleTessellator")):
            tessellator = polytope.TriangleTessellator()
        else:
            tessellator = polytope.BoostTessellator()
        time_dicts = time_tessellation(points, tessellator)
        if (rank == root):
            smallest = min(time_dicts, key=lambda time_dict: time_dict["time"])
            print(f"The quickest used the following configuration for a time of {smallest['time']}")
            print(smallest)
