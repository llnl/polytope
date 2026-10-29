"""
This file provides python-based wrappers to use Scipy's Voronoi and Delaunay tessellators
with Polytope. NOTE: Scipy is not installed as part of the Polytope virtual python environment
by default. You must install it manually. In the install directory, run

source PolytopePy/bin/activate
python3 -m pip install scipy
"""
import polytope
import sys

class ScipyVoronoi(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)

    def name(self):
        return "ScipyVoronoi"

    def tessellateQuantizedImpl(self, quantMesh):
        from scipy.spatial import Voronoi
        points = quantMesh.getRealQPoints()
        # Highly recommend to use at least Qbb and Qc, possibly add Qz
        # QJ is also a requirement for the higher bit count
        voronoi = Voronoi(points, qhull_options="Qbb Qc QJ")
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        assembler.fillTessNodes(voronoi.vertices)
        for gens, vertices in voronoi.ridge_dict.items():
            gen0 = gens[0]
            gen1 = gens[1]
            if (max(vertices) < 0):
                assembler.addInfLines(gen0, gen1)
            else:
                if (min(vertices) >= 0):
                    assembler.addFiniteEdge(gen0, gen1, vertices[0], vertices[1])
                else:
                    vertex = max(vertices)
                    assembler.addRay(gen0, gen1, vertex)
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
