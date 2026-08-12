from PYB11Generator import *

@PYB11template("int Dimension", "RealType")
class Tessellator:
    "Abstract base class for Voronoi and Voronoi-like tessellators."

    @PYB11implementation("""[](Tessellator<%(Dimension)s, %(RealType)s>& self,
                               const py::object& points,
                               Tessellation<%(Dimension)s, %(RealType)s>& mesh) {
                                 const auto coords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(points);
                                 self.tessellate(coords, mesh);
                               }""")
    def tessellate(self,
                   points="const py::object&",
                   mesh="Tessellation<%(Dimension)s, %(RealType)s>&"):
        "Generate a tessellation from flattened coordinates or coordinate tuples."
        return "void"

    @PYB11implementation("""[](Tessellator<%(Dimension)s, %(RealType)s>& self,
                               const py::object& points,
                               const py::object& PLCpoints,
                               const PLC<%(Dimension)s>& geometry,
                               Tessellation<%(Dimension)s, %(RealType)s>& mesh) {
                                 const auto pointCoords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(points);
                                 const auto plcCoords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(PLCpoints);
                                 self.tessellate(pointCoords, plcCoords, geometry, mesh);
                               }""")
    @PYB11pycppname("tessellate")
    def tessellatePLC(self,
                      points="const py::object&",
                      PLCpoints="const py::object&",
                      geometry="const PLC<%(Dimension)s>&",
                      mesh="Tessellation<%(Dimension)s, %(RealType)s>&"):
        "Generate a tessellation inside a PLC boundary."
        return "void"

    @PYB11const
    def handlesPLCs(self):
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def name(self):
        return "std::string"

    @PYB11const
    def degeneracy(self):
        return "Point<%(Dimension)s, double>"

Tessellator2d = PYB11TemplateClass(Tessellator, template_parameters=("2", "double"))
Tessellator3d = PYB11TemplateClass(Tessellator, template_parameters=("3", "double"))
