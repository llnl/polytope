from PYB11Generator import *

@PYB11template("int Dimension")
class Partitioner:
    "Abstract partitioner for replicated quantized generator points."

    PYB11typedefs = """
  using PointType = QuantizedPoint<%(Dimension)s>;
"""

    @PYB11virtual
    @PYB11const
    def name(self):
        return "std::string"

    @PYB11const
    @PYB11implementation("""[](const Partitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyPyToVector<PointType>(points, "points");
                                 return self.computePartition(generators);
                               }""")
    def computePartition(self,
                         points="const py::object&"):
        "Return this rank's subset of identically ordered quantized generators."
        return "std::vector<PointType>"


@PYB11template("int Dimension")
class RandomPartitioner(Partitioner):
    "Deterministically assign quantized generators to ranks using a seed."

    def pyinit(self,
               seed="const std::uint64_t"):
        "Construct with a deterministic ownership seed."

@PYB11template("int Dimension")
class QuasiVoronoiPartitioner(Partitioner):
    "Assigns a random section of generators to each rank."

    def pyinit(self,
               seed="const unsigned",
               maxNRank=("const unsigned", "Communicator::getNProcs()")):
        "Construct with a seed and maximum number of ranks to use."

@PYB11template("int Dimension")
class LatticePartitioner(Partitioner):
    "Partition quantized generators into a Cartesian MPI-rank lattice."

    PYB11typedefs = "using RanksPerAxis = std::array<unsigned, %(Dimension)s>;"

    def pyinit(self,
               ranksPerAxis="const RanksPerAxis&"):
        "Construct from ranks per axis; bounds come from the initialized Quantizer."


Partitioner2d = PYB11TemplateClass(Partitioner, template_parameters="2")
Partitioner3d = PYB11TemplateClass(Partitioner, template_parameters="3")
RandomPartitioner2d = PYB11TemplateClass(RandomPartitioner, template_parameters="2")
RandomPartitioner3d = PYB11TemplateClass(RandomPartitioner, template_parameters="3")
QuasiVoronoiPartitioner2d = PYB11TemplateClass(QuasiVoronoiPartitioner, template_parameters="2")
QuasiVoronoiPartitioner3d = PYB11TemplateClass(QuasiVoronoiPartitioner, template_parameters="3")
LatticePartitioner2d = PYB11TemplateClass(LatticePartitioner, template_parameters="2")
LatticePartitioner3d = PYB11TemplateClass(LatticePartitioner, template_parameters="3")
