import polytope
import random

def make_square_facets():
    return [(0, 1), (1, 2), (2, 3), (3, 0)]

def make_test_fields(tessellation):
    "Create a zone centered field dictionary of the generator locations"
    comm = polytope.Communicator.instance()
    nranks = comm.getNProcs()
    fieldnames = ["x", "y", "z"]
    centering = polytope.FieldCentering.Cell
    points = tessellation.pointsAsTuples
    result = dict()
    result[centering.name] = dict()
    for d in range(tessellation.numDims):
        result[centering.name][fieldnames[d]] = [x[d] for x in points]
    if (nranks > 1):
        rank = comm.getRank()
        result[centering.name]["ranks"] = [rank for _ in points]
    return result

def generate_random_points(N, seed = -1, dim = 2):
    if (dim == 2):
        Q = polytope.Quantizer2d.instance()
    else:
        Q = polytope.Quantizer3d.instance()
    xmin = Q.m_xlo
    xmax = Q.m_xhi
    if (seed >= 0):
        random.seed(seed)
    L = (xmax - xmin)
    pout = []
    for _ in range(N):
        for d in range(dim):
            pout.append(xmin[d] + random.random()*L[d])
    return pout
    
