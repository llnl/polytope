from PYB11Generator import *
from PLC import PLC


@PYB11template("int Dimension")
class QuantPLC(PLC):
    """A PLC with quantized coordinates and spatial-query operations.

    The :class:`Quantizer2d` or :class:`Quantizer3d` singleton must be
    initialized before constructing or initializing a QuantPLC from real
    coordinates.
    """

    PYB11typedefs = """
  using QuantPLCType = QuantPLC<%(Dimension)s>;
  using RealPoint = typename QuantPLCType::RealPoint;
  using QuantizedPoint = polytope::QuantizedPoint<%(Dimension)s>;
"""

    def pyinit(self):
        "Default constructor."

    @PYB11implementation("""[](const py::object& points,
                               const PLC<%(Dimension)s>& plc) {
                                 return QuantPLC<%(Dimension)s>(
                                   pybind11_helpers::copyCoords<%(Dimension)s, double>(points),
                                   plc);
                               }""")
    def pyinitFromPLC(self,
                      points="const py::object&",
                      plc="const PLC<%(Dimension)s>&"):
        "Construct from a PLC and flattened coordinates or coordinate tuples."

    @PYB11implementation("""[](const py::object& points) {
                                 return QuantPLC<%(Dimension)s>(
                                   pybind11_helpers::copyCoords<%(Dimension)s, double>(points));
                               }""")
    def pyinitFromPoints(self,
                         points="const py::object&"):
        "Construct a convex PLC from flattened coordinates or coordinate tuples."

    @PYB11implementation("""[](QuantPLC<%(Dimension)s>& self,
                               const py::object& points,
                               const PLC<%(Dimension)s>& plc) {
                                 self.init(
                                   pybind11_helpers::copyCoords<%(Dimension)s, double>(points),
                                   plc);
                               }""")
    @PYB11pycppname("init")
    def initFromPLC(self,
                    points="const py::object&",
                    plc="const PLC<%(Dimension)s>&"):
        "Initialize from a PLC and flattened coordinates or coordinate tuples."
        return "void"

    @PYB11implementation("""[](QuantPLC<%(Dimension)s>& self,
                               const py::object& points) {
                                 self.init(
                                   pybind11_helpers::copyCoords<%(Dimension)s, double>(points));
                               }""")
    @PYB11pycppname("init")
    def initFromPoints(self,
                       points="const py::object&"):
        "Initialize a convex PLC from flattened coordinates or coordinate tuples."
        return "void"

    def reduce(self):
        "Discard points not referenced by facets or holes."
        return "void"

    def makeConvex(self):
        "Replace this PLC with its convex hull."
        return "void"

    def orderFacets(self):
        "Order facets into connected boundary loops."
        return "void"

    @PYB11const
    def isValid(self):
        return "bool"

    @PYB11const
    @PYB11implementation("""[](const QuantPLC<%(Dimension)s>& self,
                               const py::object& point) {
                                 return self.within(pybind11_helpers::pyToPoint<%(Dimension)s, double>(point));
                               }""")
    def within(self,
               point="const py::object&"):
        "Return whether a real coordinate tuple lies within the PLC."
        return "bool"

    @PYB11const
    @PYB11pycppname("within")
    def withinQuantized(self,
                        point="const QuantizedPoint&"):
        "Return whether a quantized point lies within the PLC."
        return "bool"

    @PYB11const
    def getRealQPoints(self):
        "Return quantized vertices as real-valued point objects."
        return "std::vector<RealPoint>"

    @PYB11const
    def getRealPoints(self):
        "Return the dequantized vertices."
        return "std::vector<RealPoint>"

    @PYB11static
    def convexPLCIntersection(a="const QuantPLC<%(Dimension)s>&",
                              b="const QuantPLC<%(Dimension)s>&"):
        "Return whether two convex QuantPLCs intersect."
        return "bool"

    points = PYB11readwrite(returnpolicy="reference_internal")
    m_reduced = PYB11readwrite()
    m_convex = PYB11readwrite()
    m_loBounds = PYB11readwrite(returnpolicy="reference_internal")
    m_hiBounds = PYB11readwrite(returnpolicy="reference_internal")


QuantPLC2d = PYB11TemplateClass(QuantPLC, template_parameters="2")
QuantPLC3d = PYB11TemplateClass(QuantPLC, template_parameters="3")
