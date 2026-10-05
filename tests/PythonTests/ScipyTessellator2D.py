import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope
import sys

class ScipyVoronoi(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)

    def name(self):
        return "ScipyVoronoi"

    def third_generator(self, voronoi, gen0, gen1, vertex):
        """
        Given the vertex for an infinite ray, find a third generator to
        use to determine the direction.
        """
        third = next(
            point_index
            for point_index, region_index in enumerate(voronoi.point_region)
            if point_index not in (gen0, gen1)
            and vertex in voronoi.regions[region_index]
        )
        return third

    def tessellateQuantizedImpl(self, quantMesh):
        from scipy.spatial import Voronoi
        Q = polytope.Quantizer2d.instance()
        points = quantMesh.getRealQPoints()
        # Highly recommend to use at least Qbb and Qc, possibly add Qz
        voronoi = Voronoi(points, qhull_options="Qbb Qc")
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        assembler.fillTessNodes(voronoi.vertices)
        for gens, vertices in voronoi.ridge_dict.items():
            gen0 = gens[0]
            gen1 = gens[1]
            if (min(vertices) >= 0):
                assembler.addFiniteEdge(gen0, gen1, vertices[0], vertices[1])
            elif (max(vertices) < 0):
                assembler.addInfLines(gen0, gen1)
            else:
                vertex = max(vertices)
                gen2 = self.third_generator(voronoi, gen0, gen1, vertex)
                assembler.addRay(gen0, gen1, vertex, gen2)
        return assembler

class ScipyDelaunay(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)

    def name(self):
        return "ScipyDelaunay"

    def tessellateQuantizedImpl(self, quantMesh):
        from scipy.spatial import Delaunay
        Q = polytope.Quantizer2d.instance()
        points = quantMesh.getRealQPoints()
        delaunay = Delaunay(points, qhull_options="Qbb Qc")
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        assembler.assembleDelaunay(delaunay.simplices, delaunay.neighbors)
        return assembler

def test_serial_2d_tessellators(Ngen):
    Q = polytope.Quantizer2d.instance()

    seed = 19001
    # Generate points around a star in the middle
    boundary = Boundary2d(10, 0.1)
    # Generate random points
    points = ptu.generate_random_points(Ngen, seed, boundary)
    # Generate a Voronoi using Boost
    tessellator = polytope.BoostTessellator()
    tess_name = tessellator.name()
    print(f"Testing unbounded {tess_name}")
    boost_mesh = polytope.Tessellation2d()
    tessellator.tessellate(points, boost_mesh)
    ptu.outputMesh2d(mesh=boost_mesh,
                     filePrefix=f"boost",
                     cycle=0,
                     time=0.,
                     numFiles=1)

    # If Scipy is available, use it to generate a Voronoi and a Delaunay
    try:
        import scipy
        tessellator = ScipyVoronoi()
        tess_name = tessellator.name()
        print(f"Testing unbounded {tess_name}")
        scipy_mesh = polytope.Tessellation2d()
        tessellator.tessellate(points, scipy_mesh)
        ptu.outputMesh2d(mesh=scipy_mesh,
                         filePrefix=f"scipy",
                         cycle=0,
                         time=0.,
                         numFiles=1)
        assert polytope.isApproxTess2d(scipy_mesh, boost_mesh, 1.E-12)
        tessellator = ScipyDelaunay()
        tess_name = tessellator.name()
        print(f"Testing unbounded {tess_name}")
        scipyd_mesh = polytope.Tessellation2d()
        tessellator.tessellate(points, scipyd_mesh)
        ptu.outputMesh2d(mesh=scipyd_mesh,
                         filePrefix=f"scipyd",
                         cycle=0,
                         time=0.,
                         numFiles=1)
        assert polytope.isApproxTess2d(scipyd_mesh, boost_mesh, 1.E-12)
    except ModuleNotFoundError:
        print("No Scipy module found. Skipping test.")
        pass

if __name__ == "__main__":
    N = int(1000)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
    test_serial_2d_tessellators(N)
