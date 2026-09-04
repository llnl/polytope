from PYB11Generator import *

@PYB11template("int Dimension")
class Partitioner:
    "Abstract partitioner for replicated generator points."

    PYB11typedefs = """
  using PointType = Point<%(Dimension)s, double>;
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
        "Return generators grouped by logical partition."
        return "std::vector<std::vector<PointType>>"

    @PYB11const
    @PYB11implementation("""[](const Partitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyPyToVector<PointType>(points, "points");
                                 return self.computeOwners(generators);
                               }""")
    def computeOwners(self,
                      points="const py::object&"):
        "Return the logical partition owner for every input generator."
        return "std::vector<unsigned>"

    @PYB11const
    def numPartitions(self):
        "Return the number of logical partitions."
        return "unsigned"

    def setNumPartitions(self,
                         numPartitions="const unsigned"):
        "Set the number of logical partitions."

    @PYB11const
    @PYB11implementation("""[](const Partitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyPyToVector<PointType>(points, "points");
                                 return self.computeLocalPartition(generators);
                               }""")
    def computeLocalPartition(self,
                              points="const py::object&"):
        "Return this rank's subset of identically ordered generators."
        return "std::vector<PointType>"

    @PYB11pycppname("computeLocalPartition")
    @PYB11const
    @PYB11implementation("""[](const Partitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyPyToVector<double>(points, "points");
                                 return self.computeLocalPartition(generators);
                               }""")
    def computeLocalPartition2(self,
                              points="const py::object&"):
        "Return this rank's subset of identically ordered generators."
        return "std::vector<PointType>"

@PYB11template("int Dimension")
class RandomPartitioner(Partitioner):
    "Deterministically assign generators to logical partitions using a seed."

    def pyinit(self,
               seed="const std::uint64_t",
               numPartitions=("const unsigned", "Communicator::getNRanks()")):
        "Construct with a deterministic ownership seed."

@PYB11template("int Dimension")
class QuasiVoronoiPartitioner(Partitioner):
    "Assigns a random section of generators to each rank."

    def pyinit(self,
               seed="const unsigned",
               numPartitions=("const unsigned", "Communicator::getNRanks()")):
        "Construct with a seed and number of logical partitions."

@PYB11template("int Dimension")
class LatticePartitioner(Partitioner):
    "Partition generators into a Cartesian lattice using the Quantizer bounds."

    PYB11typedefs = "using RanksPerAxis = std::array<unsigned, %(Dimension)s>;"

    def pyinit(self,
               ranksPerAxis="const RanksPerAxis&",
               numPartitions=("const unsigned", "Communicator::getNRanks()")):
        "Construct from ranks per axis; bounds come from the initialized Quantizer."

    def pyinit2(self,
                numPartitions=("const unsigned", "Communicator::getNRanks()")):
        "Compute an optimal number of ranks per axis."


Partitioner2d = PYB11TemplateClass(Partitioner, template_parameters="2")
Partitioner3d = PYB11TemplateClass(Partitioner, template_parameters="3")
RandomPartitioner2d = PYB11TemplateClass(RandomPartitioner, template_parameters="2")
RandomPartitioner3d = PYB11TemplateClass(RandomPartitioner, template_parameters="3")
QuasiVoronoiPartitioner2d = PYB11TemplateClass(QuasiVoronoiPartitioner, template_parameters="2")
QuasiVoronoiPartitioner3d = PYB11TemplateClass(QuasiVoronoiPartitioner, template_parameters="3")
LatticePartitioner2d = PYB11TemplateClass(LatticePartitioner, template_parameters="2")
LatticePartitioner3d = PYB11TemplateClass(LatticePartitioner, template_parameters="3")
