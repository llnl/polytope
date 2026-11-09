from PYB11Generator import *
from PLC import PLC


@PYB11template("int Dimension")
class QuantPLC(PLC):
    """A PLC with quantized coordinates and spatial-query operations.

    The caller owns the quantized vertex vector and passes it to queries.
    The supplied Quantizer must outlive this object.
    """

    PYB11typedefs = """
  using QuantPLCType = QuantPLC<%(Dimension)s>;
  using RealPoint = typename QuantPLCType::RealPoint;
  using QuantizedPoint = polytope::QuantizedPoint<%(Dimension)s>;
"""

    @PYB11keepalive(1, 2)
    def pyinit(self, quantizer="const Quantizer<%(Dimension)s>&"):
        "Construct an empty QuantPLC referencing quantizer."

    @PYB11keepalive(1, 2)
    def pyinitFromPLC(self,
                      quantizer="const Quantizer<%(Dimension)s>&",
                      points="const std::vector<QuantizedPoint>&",
                      plc="const PLC<%(Dimension)s>&"):
        "Construct from a PLC and caller-owned quantized vertices."

    @PYB11keepalive(1, 2)
    def pyinitFromPoints(self,
                         quantizer="const Quantizer<%(Dimension)s>&",
                         points="const std::vector<QuantizedPoint>&"):
        "Construct from caller-owned quantized vertices."

    @PYB11pycppname("init")
    def initFromPLC(self,
                    points="const std::vector<QuantizedPoint>&",
                    plc="const PLC<%(Dimension)s>&"):
        "Initialize from a PLC and quantized vertices."
        return "void"

    @PYB11pycppname("init")
    def initFromPoints(self,
                       points="const std::vector<QuantizedPoint>&"):
        "Initialize from quantized vertices."
        return "void"

    @PYB11static
    @PYB11implementation("""[](const Quantizer<%(Dimension)s>& quantizer, const py::object& points) {
                                 return QuantPLC<%(Dimension)s>::quantizePoints(
                                   quantizer, pybind11_helpers::copyCoords<%(Dimension)s, double>(points));
                               }""")
    def quantizePoints(quantizer="const Quantizer<%(Dimension)s>&", points="const py::object&"):
        "Quantize real coordinates into a separate vertex vector."
        return "std::vector<QuantizedPoint>"

    def reduce(self, points="std::vector<QuantizedPoint>&"):
        "Discard points not referenced by facets or holes."
        return "void"

    def makeConvex(self, points="const std::vector<QuantizedPoint>&"):
        "Replace this PLC with the convex hull of the supplied vertices."
        return "void"

    def orderFacets(self):
        "Order facets into connected boundary loops."
        return "void"

    @PYB11const
    def isValid(self, points="const std::vector<QuantizedPoint>&"):
        return "bool"

    @PYB11const
    @PYB11implementation("""[](const QuantPLC<%(Dimension)s>& self,
                               const py::object& point,
                               const std::vector<QuantizedPoint>& points) {
                                 return self.within(pybind11_helpers::pyToPoint<%(Dimension)s, double>(point), points);
                               }""")
    def within(self,
               point="const py::object&",
               points="const std::vector<QuantizedPoint>&"):
        "Return whether a real coordinate tuple lies within the PLC."
        return "bool"

    @PYB11const
    @PYB11pycppname("within")
    def withinQuantized(self,
                        point="const QuantizedPoint&",
                        points="const std::vector<QuantizedPoint>&"):
        "Return whether a quantized point lies within the PLC."
        return "bool"

    @PYB11static
    def convexPLCIntersection(a="const QuantPLC<%(Dimension)s>&",
                              aPoints="const std::vector<QuantizedPoint>&",
                              b="const QuantPLC<%(Dimension)s>&",
                              bPoints="const std::vector<QuantizedPoint>&"):
        "Return whether two convex QuantPLCs intersect."
        return "bool"

    m_reduced = PYB11readwrite()
    m_convex = PYB11readwrite()
    m_loBounds = PYB11readwrite(returnpolicy="reference_internal")
    m_hiBounds = PYB11readwrite(returnpolicy="reference_internal")


QuantPLC2d = PYB11TemplateClass(QuantPLC, template_parameters="2")
QuantPLC3d = PYB11TemplateClass(QuantPLC, template_parameters="3")
