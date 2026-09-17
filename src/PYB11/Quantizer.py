from PYB11Generator import *

@PYB11singleton
@PYB11template("int Dimension")
class Quantizer:
    "Singleton for quantizing and dequantizing points."

    PYB11typedefs = """
  using QuantizerType = Quantizer<%(Dimension)s>;
  using RealType = typename QuantizerType::RealType;
  using PointType = QuantizedPoint<%(Dimension)s>;
  using RealPoint = typename QuantizerType::RealPoint;
"""

    @PYB11static
    @PYB11returnpolicy("reference")
    def instance(self):
        return "QuantizerType&"

    @PYB11pycppname("init")
    def initBounds(self,
                   xlo="const RealPoint&",
                   xhi="const RealPoint&",
                   pad=("const RealType&", "-1.0")):
        return "void"

    @PYB11pycppname("init")
    @PYB11implementation("""[](QuantizerType& self,
                               const py::object& points,
                               const RealType& pad) {
                                 const auto coords = pybind11_helpers::copyCoords<%(Dimension)s, RealType>(points);
                                 self.init(coords, pad);
                               }""")
    def initPoints(self,
                   points="const py::object&",
                   pad=("const RealType&", "-1.0")):
        return "void"

    @PYB11const
    def quantize(self,
                 x="const RealPoint&"):
        return "PointType"

    @PYB11const
    def dequantize(self,
                   x="const PointType&"):
        return "RealPoint"

    @PYB11const
    @PYB11implementation("[](const QuantizerType& self, const PointType& x) { return pybind11_helpers::keyToPy<%(Dimension)s>(self.encode(x)); }")
    def encode(self,
               x="const PointType&"):
        return "py::object"

    @PYB11const
    @PYB11implementation("[](const QuantizerType& self, const RealPoint& x) { return pybind11_helpers::keyToPy<%(Dimension)s>(self.quantizeAndEncode(x)); }")
    def quantizeAndEncode(self,
                          x="const RealPoint&"):
        return "py::object"

    @PYB11const
    @PYB11implementation("[](const QuantizerType& self, const py::object& h) { return self.decode(pybind11_helpers::pyToKey<%(Dimension)s>(h)); }")
    def decode(self,
               h="const py::object&"):
        return "PointType"

    @PYB11const
    @PYB11implementation("[](const QuantizerType& self, const py::object& h) { return self.decodeAndDequantize(pybind11_helpers::pyToKey<%(Dimension)s>(h)); }")
    def decodeAndDequantize(self,
                            h="const py::object&"):
        return "RealPoint"

    @PYB11const
    def keyEncoding(self):
        return "KeyEncoding"

    def useMortonEncoding(self):
        return "void"

    def usePackedEncoding(self):
        return "void"

    @PYB11const
    def keyName(self):
        return "const std::string"

    @PYB11const
    def degeneracy(self):
        return "RealPoint"

    @PYB11const
    def inBounds(self,
                 point="const RealPoint&"):
        return "bool"

    @PYB11const
    def inQBounds(self,
                  point="const PointType&"):
        return "bool"

    area = PYB11property(getter="area", doc="Padded physical area (2D) of quantized space")
    volume = PYB11property(getter="volume", doc="Padded physical volume (3D) of quantized space")
    padding = PYB11property(getter="getPadding", setter="setPadding",
                            doc="Padding for quantized space")

    domLength = PYB11property(getter="domLength", doc="Padded domain length")
    domLo = PYB11property(getter="domLo", doc="Padded lower bound")
    domHi = PYB11property(getter="domHi", doc="Padded upper bound")
    dx = PYB11property(getter="dx", doc="Physical grid spacing of padded quantized")
    domLoOrig = PYB11property(getter="domLoOrig", doc="Original lower bound without padding")
    domHiOrig = PYB11property(getter="domHiOrig", doc="Original upper bound without padding")
    maxBound = PYB11readonly()
    minBound = PYB11readonly()
    m_init = PYB11readonly(doc="Whether quantizer has been initialized")

Quantizer2d = PYB11TemplateClass(Quantizer, template_parameters="2")
Quantizer3d = PYB11TemplateClass(Quantizer, template_parameters="3")
