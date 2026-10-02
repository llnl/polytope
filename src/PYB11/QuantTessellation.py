from PYB11Generator import *

@PYB11template("int Dimension")
class QuantTessellation:
    """Generator data supplied to a serial quantized tessellation backend.

    This binding intentionally exposes the generator and quantization queries
    useful to ``SerialTessellator::tessellateQuantizedImpl``.  Mesh assembly,
    clipping, and distributed-mesh operations remain C++ implementation
    details.
    """

    PYB11typedefs = """
  using QuantTessellationType = QuantTessellation<%(Dimension)s>;
  using RealPoint = typename QuantTessellationType::RealPoint;
  using QuantizedPoint = polytope::QuantizedPoint<%(Dimension)s>;
  using QuantCoord = polytope::QuantizedCoordinate<%(Dimension)s>;
"""

    def pyinit(self):
        "Construct an empty quantized tessellation."

    @PYB11implementation("""[](const py::object& points) {
                                 return QuantTessellation<%(Dimension)s>(
                                   pybind11_helpers::copyCoords<%(Dimension)s, double>(points));
                               }""")
    def pyinitFromPoints(self,
                         points="const py::object&"):
        "Construct from flattened coordinates or coordinate tuples."

    @PYB11implementation("""[](QuantTessellation<%(Dimension)s>& self,
                               const py::object& points) {
                                 self.init(
                                   pybind11_helpers::copyCoords<%(Dimension)s, double>(points));
                               }""")
    @PYB11pycppname("init")
    def initFromPoints(self,
                       points="const py::object&"):
        "Replace generators from flattened coordinates or coordinate tuples."
        return "void"

    def clear(self):
        return "void"

    def disableSort(self):
        "Preserve caller-provided generator order."
        return "void"

    def enableSort(self):
        "Enable deterministic generator sorting by quantized key."
        return "void"

    @PYB11const
    def sortEnabled(self):
        return "bool"

    def sortByHash(self):
        "Sort generators by their encoded quantized keys."
        return "void"

    @PYB11const
    @PYB11implementation("""[](const QuantTessellation<%(Dimension)s>& self) {
                                 return pybind11_helpers::pointsAsTuples<%(Dimension)s, QuantCoord>(self.getQuantizedPoints());
                               }""")
    def getQuantizedPoints(self):
        "Return generator coordinates in quantized space."
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const QuantTessellation<%(Dimension)s>& self) {
                                 return pybind11_helpers::pointsAsTuples<%(Dimension)s, double>(self.getRealQPoints());
                               }""")
    def getRealQPoints(self):
        "Return nested list of quantized generator coordinates cast to doubles, not dequantized."
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const QuantTessellation<%(Dimension)s>& self) {
                                 return pybind11_helpers::pointsAsTuples<%(Dimension)s, double>(self.getRealPoints());
                               }""")
    def getRealPoints(self):
        "Return generator coordinates dequantized to physical space."
        return "py::list"

    @PYB11const
    def keyEncoding(self):
        "Return the key encoding captured when this instance was created."
        return "KeyEncoding"

    numGenerators = PYB11property(
        getterraw="[](const QuantTessellation<%(Dimension)s>& self) { return self.points.size(); }",
        doc="Number of quantized generator points.")

    points = PYB11readwrite(returnpolicy="reference_internal")
    nodes = PYB11readwrite(returnpolicy="reference_internal")
    cells = PYB11readwrite(returnpolicy="reference_internal")
    faces = PYB11readwrite(returnpolicy="reference_internal")

QuantTessellation2d = PYB11TemplateClass(QuantTessellation, template_parameters="2")
