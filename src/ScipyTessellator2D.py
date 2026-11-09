"""
This file provides python-based wrappers to use Scipy's Voronoi and Delaunay tessellators
with Polytope.
NOTE: Scipy tessellators do not work with the higher 2D bit count configuration.
"""
import polytope
import sys

class ScipyVoronoi(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)

    def name(self):
        return "ScipyVoronoi"

    def crossp(self, p0, p1):
        return p0[0]*p1[1] - p0[1]*p1[0]

    def finite_ridge_orientation(self, voronoi, gen0, gen1, v0, v1):
        "Use the cross product to determine which is the left generator"
        a, b = voronoi.vertices[v0], voronoi.vertices[v1]
        edge = b - a
        gen_vec = voronoi.points[gen1] - voronoi.points[gen0]
        cross = self.crossp(edge, gen_vec)
        if (cross == 0):
            raise RuntimeError("Invalid voronoi")
        if (cross > 0):
            return gen1, gen0
        return gen0, gen1

    def infinite_ridge_orientation(self, voronoi, points, gen0, gen1, vertex):
        "Find a third generator and set which generator is to the left"
        gens = (gen0, gen1)
        p0 = points[gen0]
        p1 = points[gen1]
        # Find a neighboring ridge sharing the finite vertex
        for (a, b), verts in voronoi.ridge_dict.items():
            # If this is the same edge, skip it
            if (a in gens and b in gens):
                continue
            # If this edge does not have the vertex, skip it
            if (vertex not in verts):
                continue
            if a in gens or b in gens:
                # Select the shared generator to be the left gen
                shared_point = gen0 if gen0 in (a, b) else gen1
                third_gen = b if a in gens else a
                d0 = [p1[d] - p0[d] for d in range(2)]
                d1 = [points[third_gen][d] - p0[d] for d in range(2)]
                cross = self.crossp(d0, d1)
                # Make sure the 3 generators are not collinear
                if (cross == 0):
                    continue
                return shared_point, third_gen
        raise ValueError("No suitable neighboring ridge found")

    def tessellateQuantizedImpl(self, quantMesh):
        from scipy.spatial import Voronoi
        points = quantMesh.getRealQPoints()
        # Recommend using Qbb and Qc for Qhull options. Do not use QJ
        voronoi = Voronoi(points, qhull_options="Qbb Qc")
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        assembler.useOrderedGens()
        assembler.fillTessNodes(voronoi.vertices)
        for (gen0, gen1), (v0, v1) in voronoi.ridge_dict.items():
            if (v0 >= 0 and v1 >= 0):
                if (v0 == v1):
                    continue
                ogens = self.finite_ridge_orientation(voronoi, gen0, gen1, v0, v1)
                assembler.addFiniteEdge(ogens[0], ogens[1], v0, v1)
            elif (v0 >= 0 or v1 >= 0):
                vertex = max([v0, v1])
                left_gen, third_gen = self.infinite_ridge_orientation(voronoi, points, gen0, gen1, vertex)
                other_gen = gen1 if gen0 == left_gen else gen0
                assembler.addRay(left_gen, other_gen, vertex, third_gen)
            else:
                assembler.addInfLines(gen0, gen1)
        return assembler

class ScipyDelaunay(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)

    def name(self):
        return "ScipyDelaunay"

    def tessellateQuantizedImpl(self, quantMesh):
        from scipy.spatial import Delaunay
        points = quantMesh.getRealQPoints()
        delaunay = Delaunay(points, qhull_options="Qbb Qc")
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        assembler.assembleDelaunay(delaunay.simplices, delaunay.neighbors)
        return assembler
