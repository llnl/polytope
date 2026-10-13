from PYB11Generator import *

@PYB11template("int Dimension")
class Partitioner:
    "Abstract collective partitioner."

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
                                 const auto generators = pybind11_helpers::copyCoords<%(Dimension)s, double>(points);
                                 return pybind11_helpers::pointsAsTuples<%(Dimension)s, double>(self.partition(generators));
                               }""")
    def partition(self,
                  points="const py::object&"):
        "Collectively return this rank's local generators."
        return "py::list"

    nparts = PYB11property(getter="getNumPartitions",
                           doc="Number of logical output partitions")

@PYB11template("int Dimension")
class ReplicatedPartitioner(Partitioner):
    "Base for partitioners that require identical generator lists on all ranks."

    PYB11typedefs = """
  using PointType = Point<%(Dimension)s, double>;
"""

    @PYB11const
    @PYB11implementation("""[](const ReplicatedPartitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyCoords<%(Dimension)s, double>(points);
                                 return pybind11_helpers::nestedPointsAsTuples<%(Dimension)s, double>(self.computePartition(generators));
                               }""")
    def computePartition(self,
                         points="const py::object&"):
        "Return generators grouped by logical partition."
        return "py::list"

    @PYB11const
    @PYB11implementation("""[](const ReplicatedPartitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyCoords<%(Dimension)s, double>(points);
                                 return self.computeOwners(generators);
                               }""")
    def computeOwners(self,
                      points="const py::object&"):
        "Return the logical partition owner for every input generator."
        return "std::vector<unsigned>"

    @PYB11const
    @PYB11implementation("""[](const ReplicatedPartitioner<%(Dimension)s>& self,
                               const py::object& points) {
                                 const auto generators = pybind11_helpers::copyCoords<%(Dimension)s, double>(points);
                                 return pybind11_helpers::pointsAsTuples<%(Dimension)s, double>(self.computeLocalPartition(generators));
                               }""")
    def computeLocalPartition(self,
                              points="const py::object&"):
        "Return this rank's subset of identically ordered generators."
        return "py::list"

    nparts = PYB11property(getter="getNumPartitions", setter="setNumPartitions",
                           doc="Number of partitions to distribute points over")

@PYB11template("int Dimension")
class QuasiVoronoiPartitioner(ReplicatedPartitioner):
    "Assigns a random section of generators to each rank."

    def pyinit(self,
               seed="const unsigned",
               numPartitions=("const unsigned", "Communicator::getNRanks()"),
               Niter=("const unsigned", "100")):
        "Construct with a seed and number of logical partitions."

    niter = PYB11property(getter="getNumIter", setter="setNumIter",
                          doc="Maximum number of iterations to run Lloyd's algorithm")

@PYB11template("int Dimension")
class LatticePartitioner(ReplicatedPartitioner):
    "Partition generators into a Cartesian lattice using the Quantizer bounds."

    PYB11typedefs = "using RanksPerAxis = std::array<unsigned, %(Dimension)s>;"

    def pyinit(self,
               ranksPerAxis="const RanksPerAxis&",
               numPartitions=("const unsigned", "Communicator::getNRanks()")):
        "Construct from ranks per axis; bounds come from the initialized Quantizer."

    def pyinit2(self,
                numPartitions=("const unsigned", "Communicator::getNRanks()")):
        "Compute an optimal number of ranks per axis."

@PYB11template("int Dimension")
class DistributedLloydPartitioner(Partitioner):
    "Collectively redistribute rank-local generators using lattice-seeded Lloyd iterations."

    PYB11typedefs = "using PointType = Point<%(Dimension)s, double>;"

    def pyinit(self,
               Niter=("const unsigned", "100")):
        "Construct with a number of Lloyd iterations."

    @PYB11const
    def name(self):
        return "std::string"

    niter = PYB11property(getter="getNumIter", setter="setNumIter",
                          doc="Number of lattice-seeded Lloyd iterations")


Partitioner2d = PYB11TemplateClass(Partitioner, template_parameters="2")
Partitioner3d = PYB11TemplateClass(Partitioner, template_parameters="3")
ReplicatedPartitioner2d = PYB11TemplateClass(ReplicatedPartitioner, template_parameters="2")
ReplicatedPartitioner3d = PYB11TemplateClass(ReplicatedPartitioner, template_parameters="3")
QuasiVoronoiPartitioner2d = PYB11TemplateClass(QuasiVoronoiPartitioner, template_parameters="2")
QuasiVoronoiPartitioner3d = PYB11TemplateClass(QuasiVoronoiPartitioner, template_parameters="3")
LatticePartitioner2d = PYB11TemplateClass(LatticePartitioner, template_parameters="2")
LatticePartitioner3d = PYB11TemplateClass(LatticePartitioner, template_parameters="3")
DistributedLloydPartitioner2d = PYB11TemplateClass(DistributedLloydPartitioner, template_parameters="2")
DistributedLloydPartitioner3d = PYB11TemplateClass(DistributedLloydPartitioner, template_parameters="3")
