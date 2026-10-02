from PYB11Generator import *

@PYB11template("int Dimension")
class PLC:
    "A piecewise linear complex in 3D, or planar straight line graph in 2D."

    def pyinit(self):
        "Default constructor"

    def clear(self):
        "Clear facets and holes."
        return "void"

    @PYB11const
    def empty(self):
        return "bool"

    @PYB11const
    def valid(self):
        return "bool"

    @PYB11implementation("[](const PLC<%(Dimension)s>& self) { std::stringstream ss; ss << self; return ss.str(); }")
    def __str__(self):
        return "std::string"

    facets = PYB11property(getterraw="[](PLC<%(Dimension)s>& self) -> std::vector<std::vector<unsigned>>& { return self.facets; }",
                           setterraw="[](PLC<%(Dimension)s>& self, const py::object& value) { self.facets = pybind11_helpers::copyFacetList(value, \"facets\"); }",
                           returnpolicy="reference_internal")
    holes = PYB11property(getterraw="[](PLC<%(Dimension)s>& self) -> std::vector<std::vector<std::vector<unsigned>>>& { return self.holes; }",
                          setterraw="[](PLC<%(Dimension)s>& self, const py::object& value) { self.holes = pybind11_helpers::copyHoleList(value, \"holes\"); }",
                          returnpolicy="reference_internal")

PLC2d = PYB11TemplateClass(PLC, template_parameters="2")
PLC3d = PYB11TemplateClass(PLC, template_parameters="3")
