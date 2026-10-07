from PYB11Generator import *

@PYB11template("int Dimension")
class VoronoiAssembler:
    """Collect raw Voronoi edges and assemble a quantized tessellation.

    The assembler retains a reference to the supplied ``QuantTessellation``;
    keep that object alive for at least as long as the assembler.
    """

    @PYB11keepalive(1, 2)
    def pyinit(self,
               input="QuantTessellation<%(Dimension)s>&"):
        "Construct an assembler for a quantized tessellation."

    @PYB11implementation("""[](VoronoiAssembler<%(Dimension)s>& self,
                               const py::object& triangles,
                               const py::object& neighbors) {
                                 auto trivec = pybind11_helpers::copyPyToTriList(triangles, "assembleDelaunay");
                                 auto nvec = pybind11_helpers::copyPyToTriList(neighbors, "assembleDelaunay");
                                 self.assembleDelaunay(trivec, nvec);
                                }""")
    def assembleDelaunay(self,
                         triangles="const py::object&",
                         neighbors="const py::object&"):
        return "void"

    @PYB11implementation("""[](VoronoiAssembler<%(Dimension)s>& self,
                               const py::object& nodes) {
                                 auto coords = pybind11_helpers::copyCoords<%(Dimension)s, double>(nodes);
                                 self.fillTessNodes(coords);
                               }""")
    def fillTessNodes(self,
                      nodes="const py::object&"):
        "Set tessellator nodes from flattened coordinates or coordinate tuples."
        return "void"

    def addFiniteEdge(self,
                      gen0="const int",
                      gen1="const int",
                      vertex0="const int",
                      vertex1="const int"):
        "Add a finite edge given the generator points and start and stop vertices"
        return "void"

    def addInfLines(self,
                    gen0="const int",
                    gen1="const int"):
        return "void"

    def addRay(self,
               gen0="const int",
               gen1="const int",
               vertex0="const int",
               dir="const Point<%(Dimension)s, double>&"):
        "Add a ray with a given direction."
        return "void"

    @PYB11pycppname("addRay")
    def addThirdPointRay(self,
                         gen0="const int",
                         gen1="const int",
                         vertex0="const int",
                         gen2="const int"):
        "Add a ray with third generator to determine direction."
        return "void"

    def addStartRay(self,
                    gen0="const int",
                    gen1="const int",
                    vertex0="const int"):
        "Add a ray that goes to infinity and starts at a point."
        return "void"

    def addEndRay(self,
                  gen0="const int",
                  gen1="const int",
                  vertex1="const int"):
        "Add a ray that starts at infinity and ends at a point."
        return "void"

    def assemble(self):
        return "void"

VoronoiAssembler2d = PYB11TemplateClass(VoronoiAssembler, template_parameters="2")
