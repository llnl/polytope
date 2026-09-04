from PYB11Generator import *

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
                            fieldMap[fieldName] = pybind11_helpers::copyPyToVector<double>(fieldItem.second, fieldName);
                          }
                          return fieldMap;
                        };

                        for (const auto fieldTypeItem: fieldsDict) {
                          const auto fieldType = FieldCenteringMap[fieldTypeItem.first.cast<std::string>()];
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
