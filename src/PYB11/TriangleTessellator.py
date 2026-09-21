from PYB11Generator import *
from SerialTessellator import SerialTessellator

@PYB11template()
@PYB11template_dict({"Dimension": "2", "RealType": "double"})
class TriangleTessellator(SerialTessellator):
    "2D Voronoi tessellator backed by Triangle."

    def pyinit(self):
        "Default constructor"

    @PYB11virtual
    @PYB11const
    def name(self):
        return "std::string"

    @PYB11virtual
    @PYB11const
    def tessellateQuantizedImpl(self,
                                input="const QuantTessellation<%(Dimension)s>&"):
        return "std::vector<std::vector<VoronoiPrimitive<%(Dimension)s>>>"
