from PYB11Generator import *

FieldCentering = PYB11enum(("Node", "Edge", "Face", "Cell"),
                           export_values=True)

@PYB11template("int Dimension", "typename TessType")
class SiloWriter:
    "Write a tessellation to a Silo file."

    @PYB11keepalive(1, 2)
    def pyinit(self,
               mesh="const %(TessType)s&"):
        "Construct a writer for mesh."

    @PYB11template("FieldType")
    @PYB11implementation("""[](SiloWriter<%(Dimension)s>& self,
                               const FieldCentering& centering,
                               const std::string& name,
                               const py::object& values) {
                                 const auto vec_vals = pybind11_helpers::copyPyToVector<%(FieldType)s>(values, "vec_vals");
                                 self.addField<%(FieldType)s>(name, vec_vals);
                               };""")
    def addField(self,
                 centering="const FieldCentering&",
                 name="const std::string&",
                 values="const py::object&"):
        return "void"
        
    def write(self,
              filePrefix="const std::string&",
              directory="const std::string&",
              cycle="int",
              time="double",
              numFiles=("int", "-1")):
        "Write with explicit directory, cycle, and time metadata."
        return "void"

    @PYB11pycppname("write")
    def writeWithoutCycle(self,
                          filePrefix="const std::string&",
                          directory="const std::string&",
                          numFiles=("int", "-1")):
        "Write without cycle or time metadata."
        return "void"

    @PYB11pycppname("write")
    def writeWithoutDir(self,
                        filePrefix="const std::string&",
                        cycle="int",
                        time="double",
                        numFiles=("int", "-1")):
        "Write without directory"
        return "void"

    @PYB11pycppname("write")
    def writeSimple(self,
                    filePrefix="const std::string&",
                    numFiles=("int", "-1")):
        "Write without directory, cycle, or time metadata."
        return "void"

    def generateTestVars(self):
        return "void"

    ovlType = PYB11property(getter="getOvlType", setter="setOvlType",
                            doc="Whether to write an Overlink file type")

SiloWriter2d = PYB11TemplateClass(
    SiloWriter,
    template_parameters=("2", "Tessellation<2, double>"))
SiloWriter3d = PYB11TemplateClass(
    SiloWriter,
    template_parameters=("3", "Tessellation<3, double>"))
