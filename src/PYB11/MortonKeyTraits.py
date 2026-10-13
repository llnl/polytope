from PYB11Generator import *

@PYB11template("int Dimension")
class MortonKeyTraits:
    "Morton key traits for quantized coordinates."

    PYB11typedefs = """
  using MortonKeyTraitsType = MortonKeyTraits<%(Dimension)s>;
  using PointType = QuantizedPoint<%(Dimension)s>;
"""
    numDim = PYB11property(constexpr=True, static=True, doc="Number of dimensions")
    coordinateBits = PYB11property(constexpr=True, static=True, doc="Number of bits in a single direction")

    def pyinit(self):
        "Default constructor"
        return

    @PYB11static
    @PYB11implementation("[]() { return py::int_(MortonKeyTraitsType::maxCoordinate()); }")
    def maxCoordinate(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[]() { return pybind11_helpers::keyToPy<%(Dimension)s>(MortonKeyTraitsType::maxKey()); }")
    def maxKey(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[](const PointType& point) { return pybind11_helpers::keyToPy<%(Dimension)s>(MortonKeyTraitsType::encode(point)); }")
    def encode(self, point="const PointType&"):
        return "py::object"

    @PYB11static
    @PYB11implementation("[](const py::object& key) { return MortonKeyTraitsType::decode(pybind11_helpers::pyToKey<%(Dimension)s>(key)); }")
    def decode(self, key="const py::object&"):
        return "PointType"

MortonKeyTraits2d = PYB11TemplateClass(MortonKeyTraits, template_parameters="2")
MortonKeyTraits3d = PYB11TemplateClass(MortonKeyTraits, template_parameters="3")
