from PYB11Generator import *
from Tessellator import Tessellator

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

DistributedTessellator2d = PYB11TemplateClass(DistributedTessellator, template_parameters="2")
DistributedTessellator3d = PYB11TemplateClass(DistributedTessellator, template_parameters="3")
