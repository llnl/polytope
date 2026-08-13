from PYB11Generator import *

@PYB11template("int Dimension", "RealType")
class Tessellation:
    "A topologically-consistent arbitrary polygonal/polyhedral mesh."

    PYB11typedefs = """
  typedef typename Cell<%(Dimension)s, %(RealType)s>::CellType CellType;
"""

    numDims = PYB11property(constexpr=True, static=True, doc="Number of dimensions")

    def pyinit(self):
        "Default constructor"

    def clear(self):
        return "void"

    @PYB11const
    def empty(self):
        return "bool"

    def computeNodeCells(self):
        return "std::vector<std::set<unsigned>>"

    def computeCellToNodes(self):
        return "std::vector<std::set<unsigned>>"

    def computeFaceCells(self):
        return "void"

    @PYB11const
    def getCell(self, cellIndex="const unsigned"):
        return "CellType"

    @PYB11implementation("[](const Tessellation<%(Dimension)s, %(RealType)s>& self) { std::stringstream ss; ss << self; return ss.str(); }")
    def __str__(self):
        return "std::string"

    points = PYB11readwrite(returnpolicy="reference_internal")
    nodes = PYB11readwrite(returnpolicy="reference_internal")
    cells = PYB11readwrite(returnpolicy="reference_internal")
    faces = PYB11readwrite(returnpolicy="reference_internal")
    boundaryNodes = PYB11readwrite(returnpolicy="reference_internal")
    boundaryFaces = PYB11readwrite(returnpolicy="reference_internal")
    faceCells = PYB11readwrite(returnpolicy="reference_internal")
    convexHull = PYB11readwrite(returnpolicy="reference_internal")
    cellRank = PYB11readwrite(returnpolicy="reference_internal")
    neighborDomains = PYB11readwrite(returnpolicy="reference_internal")
    sharedNodes = PYB11readwrite(returnpolicy="reference_internal")
    sharedFaces = PYB11readwrite(returnpolicy="reference_internal")

    pointsAsTuples = PYB11property(getterraw="[](const Tessellation<%(Dimension)s, %(RealType)s>& self) { return pybind11_helpers::pointsAsTuples<%(Dimension)s, %(RealType)s>(self.points); }")
    nodesAsTuples = PYB11property(getterraw="[](const Tessellation<%(Dimension)s, %(RealType)s>& self) { return pybind11_helpers::pointsAsTuples<%(Dimension)s, %(RealType)s>(self.nodes); }")

    zoneNodes = PYB11property(getterraw="""[](const Tessellation<%(Dimension)s, %(RealType)s>& self) -> std::vector<std::vector<int>> {
                                             const auto nzones = self.cells.size();
                                             std::vector<std::vector<int>> result(nzones);
                                             if (%(Dimension)s == 2) {
                                               for (auto izone = 0u; izone < nzones; ++izone) {
                                                 std::transform(self.cells[izone].begin(), self.cells[izone].end(), std::back_inserter(result[izone]),
                                                                [&](const int iface) { return iface < 0 ? self.faces[~iface][1] : self.faces[iface][0]; });
                                               }
                                             } else {
                                               for (auto izone = 0u; izone < nzones; ++izone) {
                                                 for (auto iface: self.cells[izone]) {
                                                   iface = iface < 0 ? ~iface : iface;
                                                   std::copy(self.faces[iface].begin(), self.faces[iface].end(), std::back_inserter(result[izone]));
                                                 }
                                                 std::sort(result[izone].begin(), result[izone].end());
                                                 result[izone].erase(std::unique(result[izone].begin(), result[izone].end()), result[izone].end());
                                               }
                                             }
                                             return result;
                                           }""")

Tessellation2d = PYB11TemplateClass(Tessellation, template_parameters=("2", "double"))
Tessellation3d = PYB11TemplateClass(Tessellation, template_parameters=("3", "double"))
