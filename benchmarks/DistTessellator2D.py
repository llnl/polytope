import sys, os, math, time
ptu_path = os.path.join(os.path.dirname(__file__), "../tests/PythonTests")
sys.path.append(ptu_path)
import polytope_test_utilities as ptu
import polytope

# Name of timer file to write
timer_log = "timer.log"

def time_tessellation(allpoints, tessellator, partitioner):
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNRanks()
    seed = 1049600
    ptess = polytope.DistributedTessellator2d(tessellator)
    time_dicts = []
    for exchange_points in [False, True]:
        ptess.setExchangePoints(exchange_points)
        for encodingtype in range(2):
            if encodingtype == 0:
                ptess.quantizer().useMortonEncoding()
            else:
                ptess.quantizer().usePackedEncoding()
            comm.Barrier()
            tess_begin = time.perf_counter()
            mesh = ptess.partitionAndTessellate(points, partitioner)
            comm.Barrier()
            tess_time = time.perf_counter() - tess_begin
            time_dict = {"MPI exchange": "points" if ptess.exchangePoints() else "hashes",
                         "encoding": ptess.quantizer().keyName(),
                         "time": tess_time}
            time_dicts.append(time_dict)
    return time_dicts

if __name__ == "__main__":
    N = int(1E6)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    nranks = comm.getNRanks()
    seed = 1049600

    lov = 0.0
    hiv = 1.0
    lobounds = [lov, lov]
    hibounds = [hiv, hiv]
    bpoints = [(lov, lov), (hiv, lov), (hiv, hiv), (lov, hiv)]
    coord_bits = round(math.log2(polytope.QuantizedKeyTraits2d().maxCoordinate()))
    ptu.rootprint(f"Coordinate bits {coord_bits}")
    ptu.rootprint(f"Generating {N} random points with {nranks} ranks")
    if (rank == root):
        with open(timer_log, "w") as ff:
            ff.write(f"{nranks} Ranks, Coord bits {coord_bits}, {N} points\n")
    generator = ptu.PointGenerator(2, lobounds=lobounds, hibounds=hibounds)
    partseed = 1042390
    parts = [polytope.LatticePartitioner2d(),
             polytope.QuasiVoronoiPartitioner2d(partseed)]
    for disttype in range(2):
        gen_begin = time.perf_counter()
        if (disttype == 0):
            points = generator.generate_random_points(N, seed)
            distname = "uniform"
        else:
            points = generator.generate_normal_random_points(N, seed=seed)
            distname = "normal"
        ptu.rootprint(f"Generators distributed in a {distname} distribution")
        if (rank == root):
            with open(timer_log, "a") as ff:
                ff.write(f"{distname} distribution\n")
        gen_time = time.perf_counter() - gen_begin
        if (hasattr(polytope, "TriangleTessellator")):
            tessellator = polytope.TriangleTessellator()
        else:
            tessellator = polytope.BoostTessellator()
        ptu.rootprint(f"Using {tessellator.name()}")
        tessellator.setQBounds(bpoints)
        for partitioner in parts:
            ptu.rootprint(f"Using {partitioner.name()} partitioner")
            time_dicts = time_tessellation(points, tessellator, partitioner)
            if (rank == root):
                smallest = min(time_dicts, key=lambda time_dict: time_dict["time"])
                print(f"The quickest used the following configuration for a time of {smallest['time']}")
                print(smallest)
                with open(timer_log, "a") as ff:
                    ff.write(f"{partitioner.name()} partitioner\n")
                    for tt in time_dicts:
                        for key, value in tt.items():
                            ff.write(f" {key}: {value},")
                        ff.write("\n")
