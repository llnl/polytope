import polytope

def make_square_facets():
    return [(0, 1), (1, 2), (2, 3), (3, 0)]

def computeLocations(tessellation):
    "Create a zone centered field dictionary of the generator locations"
    fieldnames = ["x", "y", "z"]
    centering = polytope.FieldCentering.Cell
    points = tessellation.pointsAsTuples
    result = dict()
    result[centering.name] = dict()
    for d in range(tessellation.numDims):
        result[centering.name][fieldnames[d]] = [x[d] for x in points]
    return result
