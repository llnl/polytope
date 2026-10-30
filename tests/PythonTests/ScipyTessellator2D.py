"""
This file provides python-based wrappers to use Scipy's Voronoi and Delaunay tessellators
with Polytope. NOTE: Scipy is not installed as part of the Polytope virtual python environment
by default. You must install it manually. In the install directory, run.

source PolytopePy/bin/activate
python3 -m pip install scipy
NOTE: More work is needed to make the assembler robust enough for these tessellators
"""
import polytope
import sys

class ScipyVoronoi(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)

    def name(self):
        return "ScipyVoronoi"

    def thirdGenerator(self, a, b, vertex, vertexGens, points):
        ax, ay = points[a]
        bx, by = points[b]
        third, side = None, None

        for c in vertexGens[vertex] - {a, b}:
            cx, cy = points[c]
            winding = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if winding == 0:
                continue
            if side is not None and side != (winding > 0):
                raise ValueError("Inconsistent generator sides")
            third, side = c, winding > 0

        if third is None:
            raise ValueError("No noncollinear third generator")
        return third

    def tessellateQuantizedImpl(self, quantMesh):
        from scipy.spatial import Voronoi
        points = quantMesh.getRealQPoints()
        # Highly recommend to use at least Qbb and Qc, possibly add Qz
        # QJ is also a requirement for the higher bit count
        voronoi = Voronoi(points, qhull_options="Qbb Qc")
        vertexGens = [set() for _ in voronoi.vertices]
        for gens, vertices in voronoi.ridge_dict.items():
            for v in vertices:
                if v >= 0:
                    vertexGens[v].update(gens)
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        assembler.fillTessNodes(voronoi.vertices)
        for (gen0, gen1), (v0, v1) in voronoi.ridge_dict.items():
            if (v0 >= 0 and v1 >= 0):
                #c = self.thirdGenerator(gen0, gen1, v0, vertexGens, points)
                assembler.addFiniteEdge(gen0, gen1, vertices[0], vertices[1])#, c)
            elif (v0 >= 0 or v1 >= 0):
                vertex = max([v0, v1])
                #c = self.thirdGenerator(gen0, gen1, vertex, vertexGens, points)
                g0 = polytope.Point2d(*points[gen0])
                g1 = polytope.Point2d(*points[gen1])
                diff = g0 - g1
                dr = polytope.Point2d(-diff.y, diff.x)
                assembler.addRay(gen0, gen1, vertex, dr)
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
