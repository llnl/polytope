from PYB11Generator import *
from Tessellator import Tessellator

@PYB11template("int Dimension")
@PYB11template_dict({"RealType": "double"})
class SerialTessellator(Tessellator):
    """Abstract serial tessellator base.
    Concrete implementations populate and return a Voronoi assembler, which
    produces the quantized tessellation.
    """

    def pyinit(self):
        "Default constructor"

    @PYB11pure_virtual
    @PYB11const
    def tessellateQuantizedImpl(self,
                                input="QuantTessellation<%(Dimension)s>&"):
        return "VoronoiAssembler<%(Dimension)s>"

SerialTessellator2d = PYB11TemplateClass(SerialTessellator, template_parameters="2")
