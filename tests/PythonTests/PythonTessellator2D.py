# This test ensures python wrapped tessellators work without having to
# install a separate tessellator
import polytope_test_utilities as ptu
from Boundary2d import Boundary2d
import polytope
import sys

class FakeTessellator(polytope.SerialTessellator2d):
    def __init__(self):
        polytope.SerialTessellator2d.__init__(self)
        # Expected nodes 
        self.orig_nodes = [[2., 1.5], [6., 1.5], [4., 2.5]]
        self.nodes = []

    def name(self):
        return "FakeTessellator"

    def setNodes(self, boost_nodes):
        """
        Get nodes from the Boost mesh that are closest to the original nodes.
        Must retain the order of the original.
        """
        for target in self.orig_nodes:
            closest = min(
                boost_nodes,
                key=lambda point: sum(
                    (x - y) ** 2
                    for x, y in zip(target, point)
                )
            )
            self.nodes.append(closest)

    def tessellateQuantizedImpl(self, quantMesh):
        Q = polytope.Quantizer2d.instance()
        assembler = polytope.VoronoiAssembler2d(quantMesh)
        qnodes = [Q.quantizeReal(x) for x in self.nodes]
        assembler.fillTessNodes(qnodes)

        # The input generators form three Voronoi vertices:
        #   0: generators 0, 1, 3
        #   1: generators 1, 2, 4
        #   2: generators 1, 3, 4
        # Connect the finite edges first, then add the rays on the convex hull
        # of the generator set. The third generator establishes each ray's
        # outward direction.
        assembler.addFiniteEdge(1, 3, 0, 2)
        assembler.addFiniteEdge(1, 4, 2, 1)
        assembler.addRay(0, 1, 0, 3)
        assembler.addRay(0, 3, 0, 1)
        assembler.addRay(1, 2, 1, 4)
        assembler.addRay(2, 4, 1, 1)
        assembler.addRay(3, 4, 2, 1)
        return assembler

def test_python_2d_tessellators():
    Q = polytope.Quantizer2d.instance()

    boundary = Boundary2d(0, 0.01)
    boundary.mDiff = 10.
    boundary.mCenter = [0., 0.]
    boundary.setDefaultBoundary(0)
    # Use a specific set of generators that have an easy Voronoi
    points = [
        (0, 0),  # generator 0
        (4, 0),  # generator 1
        (8, 0),  # generator 2
        (2, 4),  # generator 3
        (6, 4),  # generator 4
    ]
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

    # Use our FakeTessellator
    tessellator = FakeTessellator()
    tessellator.setNodes(boost_mesh.nodes)
    tess_name = tessellator.name()
    print(f"Testing unbounded {tess_name}")
    mesh = polytope.Tessellation2d()
    tessellator.tessellate(points, mesh)
    ptu.outputMesh2d(mesh=mesh,
                     filePrefix=f"pyfake",
                     cycle=0,
                     time=0.,
                     numFiles=1)
    assert mesh == boost_mesh

if __name__ == "__main__":
    test_python_2d_tessellators()
