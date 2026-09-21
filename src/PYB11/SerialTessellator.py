from PYB11Generator import *
from Tessellator import Tessellator

@PYB11template("int Dimension")
@PYB11template_dict({"RealType": "double"})
class SerialTessellator(Tessellator):
    """Abstract serial tessellator base.
    Concrete implementations generate Voronoi primitives and use the shared
    Voronoi assembler to produce quantized tessellations.
    """

    def pyinit(self):
        "Default constructor"

    @PYB11pure_virtual
    @PYB11const
    def tessellateQuantizedImpl(self,
                                input="const QuantTessellation<%(Dimension)s>&"):
        return "std::vector<std::vector<VoronoiPrimitive<%(Dimension)s>>>"

SerialTessellator2d = PYB11TemplateClass(SerialTessellator, template_parameters="2")
