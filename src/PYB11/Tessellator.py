from PYB11Generator import *

@PYB11template("int Dimension", "RealType")
class Tessellator:
    "Abstract base class for Voronoi and Voronoi-like tessellators."

    @PYB11implementation("""[](Tessellator<%(Dimension)s, %(RealType)s>& self,
                               const py::object& points) {
                                 const auto coords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(points);
                                 return self.tessellate(coords);
                               }""")
    def tessellate(self,
                   points="const py::object&"):
        "Generate a tessellation from flattened coordinates or coordinate tuples."
        return "Tessellation<%(Dimension)s, %(RealType)s>"

    @PYB11implementation("""[](Tessellator<%(Dimension)s, %(RealType)s>& self,
                               const py::object& points,
                               const py::object& PLCpoints,
                               const PLC<%(Dimension)s>& geometry) {
                                 const auto pointCoords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(points);
                                 const auto plcCoords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(PLCpoints);
                                 return self.tessellate(pointCoords, plcCoords, geometry);
                               }""")
    @PYB11pycppname("tessellate")
    def tessellatePLC(self,
                      points="const py::object&",
                      PLCpoints="const py::object&",
                      geometry="const PLC<%(Dimension)s>&"):
        "Generate a tessellation inside a PLC boundary."
        return "Tessellation<%(Dimension)s, %(RealType)s>"

    @PYB11returnpolicy("reference_internal")
    def quantizer(self):
        return "Quantizer<%(Dimension)s>&"

    @PYB11implementation("""[](Tessellator<%(Dimension)s, %(RealType)s>& self,
                               const py::object& points) {
                                 const auto coords = pybind11_helpers::copyCoords<%(Dimension)s, %(RealType)s>(points);
                                 return self.setQBounds(coords);
                               }""")
    def setQBounds(self,
                   points="const py::object&"):
        "Include these points when determining quantizer bounds, without using them as generators"
        return "void"

    def setPadding(self,
                   padding="const double&"):
        "Set the percent to pad the quantized space relative to the points."
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def name(self):
        return "std::string"

    @PYB11const
    def degeneracy(self):
        return "Point<%(Dimension)s, double>"

    @PYB11const
    def clippingEnabled(self):
        return "bool"

Tessellator2d = PYB11TemplateClass(Tessellator, template_parameters=("2", "double"))
Tessellator3d = PYB11TemplateClass(Tessellator, template_parameters=("3", "double"))
