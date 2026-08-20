"Python bindings for the renovated Polytope C++ interface."

from PYB11Generator import *

PYB11includes = ['"polytope.hh"',
                 '"polytope_pybind11_helpers.hh"',
                 '"MortonKeyTraits.hh"',
                 '"Quantizer.hh"',
                 '"Point.hh"',
                 '"Cell.hh"',
                 '"PLC.hh"',
                 '"Tessellation.hh"',
                 '"Tessellator.hh"']

PYB11modulepreamble = """
// Initialize MPI
Communicator::init();

// Call these routines when module is exited
auto atexit = py::module_::import("atexit");
atexit.attr("register")(py::cpp_function([]() {
   Communicator::finalize();
}));
"""

PYB11namespaces = ["polytope"]

FieldCentering = PYB11enum(("Node", "Edge", "Face", "Cell"),
                           namespace="polytope",
                           export_values=True,
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
from MortonKeyTraits import *
from Quantizer import *
from Communicator import *

vector_of_Point2d = PYB11_bind_vector("Point<2, double>", opaque=True, local=True)
vector_of_Point3d = PYB11_bind_vector("Point<3, double>", opaque=True, local=True)
vector_of_CoordinatePoint2d = PYB11_bind_vector("polytope::QuantizedPoint<2>", opaque=True, local=True)
vector_of_CoordinatePoint3d = PYB11_bind_vector("polytope::QuantizedPoint<3>", opaque=True, local=True)
vector_of_KeyPoint2d = PYB11_bind_vector("Point<2, polytope::MortonKey<2>>", opaque=True, local=True)
vector_of_KeyPoint3d = PYB11_bind_vector("Point<3, polytope::MortonKey<3>>", opaque=True, local=True)
vector_of_vector_of_Point2d = PYB11_bind_vector("std::vector<Point<2, double>>", opaque=True, local=True)
vector_of_vector_of_Point3d = PYB11_bind_vector("std::vector<Point<3, double>>", opaque=True, local=True)
vector_of_vector_of_vector_of_Point3d = PYB11_bind_vector("std::vector<std::vector<Point<3, double>>>", opaque=True, local=True)

from PLC import *
from Tessellation import *
from Tessellator import *
from SerialTessellators import *
from DistributedTessellator import *
from SiloUtils import *
