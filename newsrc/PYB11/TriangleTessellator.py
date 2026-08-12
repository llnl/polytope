from PYB11Generator import *
from Tessellator import Tessellator

@PYB11template()
@PYB11template_dict({"Dimension": "2", "RealType": "double"})
class TriangleTessellator(Tessellator):
    "2D Voronoi tessellator backed by Triangle."

    def pyinit(self):
        "Default constructor"

    @PYB11const
    def name(self):
        return "std::string"
