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

    def extend(self,
               extendPad="const RealType"):
        "Modify the padding of the Quantizer by a certain percent of the domain length"
        return "void"

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

    m_lx_o = PYB11readwrite(returnpolicy="reference_internal")
    m_xlo_o = PYB11readwrite(returnpolicy="reference_internal")
    m_dx_o = PYB11readwrite(returnpolicy="reference_internal")
    m_xlo = PYB11readwrite(returnpolicy="reference_internal")
    m_xhi = PYB11readwrite(returnpolicy="reference_internal")
    m_pad = PYB11readwrite()
    # maxCoord = PYB11readwrite(returnpolicy="reference_internal")
    # minCoord = PYB11readwrite(returnpolicy="reference_internal")
    maxBound = PYB11readwrite(returnpolicy="reference_internal")
    minBound = PYB11readwrite(returnpolicy="reference_internal")
    rmaxBound = PYB11readwrite(returnpolicy="reference_internal")
    rminBound = PYB11readwrite(returnpolicy="reference_internal")
    m_init = PYB11readwrite()

Quantizer2d = PYB11TemplateClass(Quantizer, template_parameters="2")
Quantizer3d = PYB11TemplateClass(Quantizer, template_parameters="3")
