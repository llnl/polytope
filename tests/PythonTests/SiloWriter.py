import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope
import sys, time, math

def test_siloIO(Ngen):
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    nranks = comm.getNRanks()
    # Seed for generating points randomly
    seed = 19001
    # Seed for partitioning points randomly
    partseed = 19
    boundary = Boundary2d(0)
    Q = polytope.Quantizer2d.instance()
    xlo = Q.domLo.x
    xlen = Q.domLength.x
    # Generate Ngen points per rank
    Ntotal = Ngen*nranks
    all_points = ptu.generate_normal_random_points(Ntotal, seed=seed, boundary2d=boundary)
    serial_tessellator = polytope.BoostTessellator()
    tessellator = polytope.DistributedTessellator2d(serial_tessellator)
    # Make a partitioner
    partitioner = polytope.QuasiVoronoiPartitioner2d(partseed)
    rank_points = partitioner.computeLocalPartition(all_points)
    mesh = polytope.Tessellation2d()
    tessellator.tessellate(rank_points, mesh)
    filePrefix = "SiloIOTest"
    writer = polytope.SiloWriter2d(mesh)
    # Normalized x coordinates
    normx = [(p[0] - xlo)/xlen for p in rank_points]
    field = [math.sin(2.*math.pi*x) for x in normx]
    writer.addField(polytope.FieldCentering.Cell, "xvelocity", field)
    writer.generateTestVars()
    matnames = ["O2", "H2O"]
    matvfs = []
    erx = [6.*x - 3. for x in normx]
    for x in erx:
        mat0 = math.erf(x)
        mat1 = 1. - mat0
        matvfs.append([mat0, mat1])
    writer.addMaterials(matnames, matvfs)
    writer.ovlType = True
    writer.write(filePrefix+"ovl", nranks)

if __name__ == "__main__":
    # Provide the number of generators per rank
    N = int(1000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_siloIO(N)
