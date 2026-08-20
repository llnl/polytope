from PYB11Generator import *

@PYB11template("int Dimension")
class HashKey:
    "Morton hash traits for quantized coordinates."

    PYB11typedefs = """
  using HashKeyType = HashKey<%(Dimension)s>;
  using IntPoint = typename HashKeyType::IntPoint;
"""
    numDims = PYB11property(constexpr=True, static=True, doc="Number of dimensions")
    bits1D = PYB11property(constexpr=True, static=True, doc="Number of bits in a single direction")

    def pyinit(self):
        "Default constructor"
        return

    @PYB11static
    @PYB11implementation("[]() { return py::int_(HashKeyType::coordMax()); }")
    def coordMax(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[]() { return pybind11_helpers::coordHashToPy<%(Dimension)s>(HashKeyType::hashMax()); }")
    def hashMax(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[](const IntPoint& point) { return pybind11_helpers::coordHashToPy<%(Dimension)s>(HashKeyType::hash(point)); }")
    def hash(self, point="const IntPoint&"):
        return "py::object"

    @PYB11static
    @PYB11implementation("[](const py::object& key) { return HashKeyType::unhash(pybind11_helpers::pyToCoordHash<%(Dimension)s>(key)); }")
    def unhash(self, key="const py::object&"):
        return "IntPoint"

HashKey2d = PYB11TemplateClass(HashKey, template_parameters="2")
HashKey3d = PYB11TemplateClass(HashKey, template_parameters="3")
