from PYB11Generator import *
from Tessellator import Tessellator

@PYB11template()
@PYB11template_dict({"Dimension": "2", "RealType": "double"})
class BoostTessellator(Tessellator):
    "2D Voronoi tessellator backed by Boost.Polygon."

    def pyinit(self):
        "Default constructor"

    @PYB11virtual
    @PYB11const
    def name(self):
        return "std::string"
