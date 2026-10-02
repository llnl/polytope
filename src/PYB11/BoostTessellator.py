from PYB11Generator import *
from SerialTessellator import SerialTessellator

@PYB11template()
@PYB11template_dict({"Dimension": "2", "RealType": "double"})
class BoostTessellator(SerialTessellator):
    "2D Voronoi tessellator backed by Boost.Polygon."

    def pyinit(self):
        "Default constructor"

    @PYB11virtual
    @PYB11const
    def name(self):
        return "std::string"

    @PYB11virtual
    @PYB11const
    def tessellateQuantizedImpl(self,
                                input="QuantTessellation<%(Dimension)s>&"):
        return "VoronoiAssembler<%(Dimension)s>"
