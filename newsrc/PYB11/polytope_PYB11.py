"Python bindings for the renovated Polytope C++ interface."

from PYB11Generator import *

PYB11includes = ['"polytope.hh"',
                 '"polytope_pybind11_helpers.hh"',
                 '"HashKey.hh"',
                 '"Quantizer.hh"',
                 '"Point.hh"',
                 '"Cell.hh"',
                 '"PLC.hh"',
                 '"Tessellation.hh"',
                 '"Tessellator.hh"']

PYB11namespaces = ["polytope"]

FieldCentering = PYB11enum(("Node", "Edge", "Face", "Cell"),
                           namespace="polytope",
                           doc="Centering locations for mesh fields.")

vector_of_unsigned = PYB11_bind_vector("unsigned", opaque=True, local=True)
vector_of_int = PYB11_bind_vector("int", opaque=True, local=True)
vector_of_int64 = PYB11_bind_vector("int64_t", opaque=True, local=True)
vector_of_double = PYB11_bind_vector("double", opaque=True, local=True)
vector_of_string = PYB11_bind_vector("std::string", opaque=True, local=True)
vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<unsigned>", opaque=True, local=True)
vector_of_vector_of_int = PYB11_bind_vector("std::vector<int>", opaque=True, local=True)
vector_of_vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<std::vector<unsigned>>", opaque=True, local=True)
vector_of_vector_of_vector_of_int = PYB11_bind_vector("std::vector<std::vector<int>>", opaque=True, local=True)
vector_of_set_of_unsigned = PYB11_bind_vector("std::set<unsigned>", opaque=True, local=True)
vector_of_vector_of_set_of_unsigned = PYB11_bind_vector("std::vector<std::set<unsigned>>", opaque=True, local=True)

from Point import *
from HashKey import *
from Quantizer import *

vector_of_Point2d = PYB11_bind_vector("Point<2, double>", opaque=True, local=True)
vector_of_Point3d = PYB11_bind_vector("Point<3, double>", opaque=True, local=True)
vector_of_CoordPoint2d = PYB11_bind_vector("Point<2, polytope::HashKey<2>::IntType>", opaque=True, local=True)
vector_of_CoordPoint3d = PYB11_bind_vector("Point<3, polytope::HashKey<3>::IntType>", opaque=True, local=True)
vector_of_HashPoint2d = PYB11_bind_vector("Point<2, polytope::HashKey<2>::CoordHash>", opaque=True, local=True)
vector_of_HashPoint3d = PYB11_bind_vector("Point<3, polytope::HashKey<3>::CoordHash>", opaque=True, local=True)
vector_of_vector_of_Point2d = PYB11_bind_vector("std::vector<Point<2, double>>", opaque=True, local=True)
vector_of_vector_of_Point3d = PYB11_bind_vector("std::vector<Point<3, double>>", opaque=True, local=True)
vector_of_vector_of_vector_of_Point3d = PYB11_bind_vector("std::vector<std::vector<Point<3, double>>>", opaque=True, local=True)

from PLC import *
from Tessellation import *
from Tessellator import *
from SerialTessellators import *

@PYB11template("int Dimension")
@PYB11implementation("""[](const Tessellation<%(Dimension)s, double>& mesh,
                           const std::string& filePrefix,
                           py::dict fieldsDict,
                           const std::string& directory,
                           int cycle,
                           double time,
                           int numFiles) {
#ifdef POLYTOPE_ENABLE_SILO
                             using Writer = SiloWriter<%(Dimension)s, Tessellation<%(Dimension)s, double>>;
                             using FieldMap = typename Writer::FieldMap;
                             typename Writer::FieldTypeMap fields;

                             auto copyFieldMap = [](py::dict fieldDict) {
                               FieldMap fieldMap;
                               for (const auto fieldItem: fieldDict) {
                                 const auto fieldName = fieldItem.first.cast<std::string>();
                                 fieldMap[fieldName] = fieldItem.second.cast<std::vector<double>>();
                               }
                               return fieldMap;
                             };

                             for (const auto fieldTypeItem: fieldsDict) {
                               const auto fieldType = fieldTypeItem.first.cast<FieldCentering>();
                               const auto fieldDict = fieldTypeItem.second.cast<py::dict>();
                               fields[fieldType] = copyFieldMap(fieldDict);
                             }

                             if (fields.empty() and directory.empty() and cycle == 0 and time == 0.0) {
                               Writer::write(mesh, filePrefix, numFiles);
                             } else {
                               Writer::write(mesh, fields, filePrefix, directory, cycle, time, numFiles);
                             }
#else
                             throw std::runtime_error("Polytope built without SILO support");
#endif
                           }""")
def writeSilo(mesh="const Tessellation<%(Dimension)s, double>&",
              filePrefix="const std::string&",
              fields=("py::dict", "py::dict()"),
              directory=("const std::string&", "\"\""),
              cycle=("int", "0"),
              time=("double", "0.0"),
              numFiles=("int", "-1")):
    "Write a tessellation and optional centered fields to a Silo file. The fields dict is keyed by FieldCentering."
    return "void"

writeSilo2d = PYB11TemplateFunction(writeSilo, template_parameters="2", pyname="writeSilo")
writeSilo3d = PYB11TemplateFunction(writeSilo, template_parameters="3", pyname="writeSilo")
