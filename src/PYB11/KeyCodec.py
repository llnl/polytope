from PYB11Generator import *

@PYB11template("int Dimension")
class KeyCodec:
    "Runtime-selectable encoding for quantized-coordinate keys."

    PYB11typedefs = """
  using KeyCodecType = KeyCodec<%(Dimension)s>;
  using PointType = QuantizedPoint<%(Dimension)s>;
"""

    def pyinit(self,
               encoding=("KeyEncoding", "KeyEncoding::Morton")):
        return

    def encoding(self):
        return "KeyEncoding"

    def setEncoding(self,
                    encoding="KeyEncoding"):
        return "void"

    @PYB11const
    def keyName(self):
        return "const std::string"

    @PYB11const
    @PYB11implementation("[](const KeyCodecType& self, const PointType& point) { return pybind11_helpers::keyToPy<%(Dimension)s>(self.encode(point)); }")
    def encode(self,
               point="const PointType&"):
        return "py::object"

    @PYB11const
    @PYB11implementation("[](const KeyCodecType& self, const py::object& key) { return self.decode(pybind11_helpers::pyToKey<%(Dimension)s>(key)); }")
    def decode(self,
               key="const py::object&"):
        return "PointType"

KeyCodec2d = PYB11TemplateClass(KeyCodec, template_parameters="2")
KeyCodec3d = PYB11TemplateClass(KeyCodec, template_parameters="3")
