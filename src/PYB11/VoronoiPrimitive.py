from PYB11Generator import *

@PYB11template("int Dimension")
class VoronoiPrimitive:
    """Raw Voronoi primitive emitted by a serial tessellation backend.

    The currently supported specialization is :class:`VoronoiPrimitive2d`.
    Its generator pair is exposed as ``gp``; finite endpoints are ``rp0`` and
    ``rp1`` and their corresponding infinity flags are ``inf0`` and ``inf1``.
    """

    def pyinit(self,
               gp0="const int",
               gp1="const int"):
        "Construct a primitive between two generator indices."

    rp0 = PYB11readwrite()
    rp1 = PYB11readwrite()
    gp3 = PYB11readwrite()
    inf0 = PYB11readwrite()
    inf1 = PYB11readwrite()

VoronoiPrimitive2d = PYB11TemplateClass(VoronoiPrimitive, template_parameters="2")
