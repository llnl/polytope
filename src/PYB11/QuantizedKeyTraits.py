from PYB11Generator import *

@PYB11template("int Dimension")
class QuantizedKeyTraits:
    "Codec-neutral key storage traits for quantized coordinates."

    PYB11typedefs = """
  using QuantizedKeyTraitsType = QuantizedKeyTraits<%(Dimension)s>;
"""
    numDim = PYB11property(constexpr=True, static=True, doc="Number of dimensions")
    coordinateBits = PYB11property(constexpr=True, static=True, doc="Number of bits in a single direction")

    def pyinit(self):
        "Default constructor"
        return

    @PYB11static
    @PYB11implementation("[]() { return py::int_(QuantizedKeyTraitsType::maxCoordinate()); }")
    def maxCoordinate(self):
        return "py::object"

    @PYB11static
    @PYB11implementation("[]() { return pybind11_helpers::keyToPy<%(Dimension)s>(QuantizedKeyTraitsType::maxKey()); }")
    def maxKey(self):
        return "py::object"

QuantizedKeyTraits2d = PYB11TemplateClass(QuantizedKeyTraits, template_parameters="2")
QuantizedKeyTraits3d = PYB11TemplateClass(QuantizedKeyTraits, template_parameters="3")
