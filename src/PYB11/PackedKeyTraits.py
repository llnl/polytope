from PYB11Generator import *

@PYB11template("int Dimension")
class PackedKeyTraits:
    "Consecutive-field packed key traits for quantized coordinates."

    PYB11typedefs = """
  using PackedKeyTraitsType = PackedKeyTraits<%(Dimension)s>;
  using PointType = QuantizedPoint<%(Dimension)s>;
"""
    numDim = PYB11property(constexpr=True, static=True, doc="Number of dimensions")
    coordinateBits = PYB11property(constexpr=True, static=True, doc="Number of bits in a single direction")

    def pyinit(self):
        "Default constructor"
        return

    @PYB11static
    @PYB11implementation("[]() { return py::int_(PackedKeyTraitsType::maxCoordinate()); }")
    def maxCoordinate(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[]() { return pybind11_helpers::keyToPy<%(Dimension)s>(PackedKeyTraitsType::maxKey()); }")
    def maxKey(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[](const PointType& point) { return pybind11_helpers::keyToPy<%(Dimension)s>(PackedKeyTraitsType::encode(point)); }")
    def encode(self, point="const PointType&"):
        return "py::object"

    @PYB11static
    @PYB11implementation("[](const py::object& key) { return PackedKeyTraitsType::decode(pybind11_helpers::pyToKey<%(Dimension)s>(key)); }")
    def decode(self, key="const py::object&"):
        return "PointType"

PackedKeyTraits2d = PYB11TemplateClass(PackedKeyTraits, template_parameters="2")
PackedKeyTraits3d = PYB11TemplateClass(PackedKeyTraits, template_parameters="3")
