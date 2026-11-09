from PYB11Generator import *
from Tessellator import Tessellator
from Partitioner import Partitioner

@PYB11template("int Dimension")
@PYB11template_dict({"RealType": "double"})
class DistributedTessellator(Tessellator):
    "Distributed tessellator using a supplied serial tessellator."

    @PYB11keepalive(1, 2)
    def pyinit(self,
               serialTessellator="Tessellator<%(Dimension)s, double>&"):
        "Construct with a serial tessellator."

    @PYB11virtual
    @PYB11const
    def name(self):
        return "std::string"

    def setExchangePoints(self, value="const bool"):
        "Select direct quantized-point exchange."
        return "void"

    @PYB11const
    def exchangePoints(self):
        "Return the currently synchronized generator exchange representation."
        return "bool"

    @PYB11implementation("""[](DistributedTessellator<%(Dimension)s>& self,
                               const py::object& points,
                               const Partitioner<%(Dimension)s>& partitioner) {
                                 const auto coords = pybind11_helpers::copyCoords<%(Dimension)s, double>(points);
                                 return self.partitionAndTessellate(coords, partitioner);
                               }""")
    def partitionAndTessellate(self,
                               points="const py::object&",
                               partitioner="const Partitioner<%(Dimension)s>&"):
        "Partition replicated generators, then generate this rank's tessellation."
        return "Tessellation<%(Dimension)s, double>"

    @PYB11implementation("""[](DistributedTessellator<%(Dimension)s>& self,
                               const py::object& points,
                               const py::object& PLCpoints,
                               const PLC<%(Dimension)s>& geometry,
                               const Partitioner<%(Dimension)s>& partitioner) {
                                 const auto pointCoords = pybind11_helpers::copyCoords<%(Dimension)s, double>(points);
                                 const auto plcCoords = pybind11_helpers::copyCoords<%(Dimension)s, double>(PLCpoints);
                                 return self.partitionAndTessellate(pointCoords, plcCoords, geometry, partitioner);
                               }""")
    @PYB11pycppname("partitionAndTessellate")
    def partitionAndTessellatePLC(self,
                                  points="const py::object&",
                                  PLCpoints="const py::object&",
                                  geometry="const PLC<%(Dimension)s>&",
                                  partitioner="const Partitioner<%(Dimension)s>&"):
        "Partition replicated generators inside a PLC, then tessellate."
        return "Tessellation<%(Dimension)s, double>"

DistributedTessellator2d = PYB11TemplateClass(DistributedTessellator, template_parameters="2")
DistributedTessellator3d = PYB11TemplateClass(DistributedTessellator, template_parameters="3")
