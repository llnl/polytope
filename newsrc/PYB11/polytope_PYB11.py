"Python bindings for the renovated Polytope C++ interface."

from PYB11Generator import *

PYB11includes = ['"polytope.hh"',
                 '"polytope_pybind11_helpers.hh"',
                 '"Point.hh"',
                 '"Cell.hh"',
                 '"PLC.hh"',
                 '"Tessellation.hh"',
                 '"Tessellator.hh"']

PYB11namespaces = ["polytope"]

vector_of_unsigned = PYB11_bind_vector("unsigned", opaque=True, local=True)
vector_of_int = PYB11_bind_vector("int", opaque=True, local=True)
vector_of_double = PYB11_bind_vector("double", opaque=True, local=True)
vector_of_string = PYB11_bind_vector("std::string", opaque=True, local=True)
vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<unsigned>", opaque=True, local=True)
vector_of_vector_of_int = PYB11_bind_vector("std::vector<int>", opaque=True, local=True)
vector_of_vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<std::vector<unsigned>>", opaque=True, local=True)
vector_of_vector_of_vector_of_int = PYB11_bind_vector("std::vector<std::vector<int>>", opaque=True, local=True)
vector_of_set_of_unsigned = PYB11_bind_vector("std::set<unsigned>", opaque=True, local=True)
vector_of_vector_of_set_of_unsigned = PYB11_bind_vector("std::vector<std::set<unsigned>>", opaque=True, local=True)

from Point import *

vector_of_Point2D = PYB11_bind_vector("Point<2, double>", opaque=True, local=True)
vector_of_Point3D = PYB11_bind_vector("Point<3, double>", opaque=True, local=True)
vector_of_vector_of_Point2D = PYB11_bind_vector("std::vector<Point<2, double>>", opaque=True, local=True)
vector_of_vector_of_Point3D = PYB11_bind_vector("std::vector<Point<3, double>>", opaque=True, local=True)
vector_of_vector_of_vector_of_Point3D = PYB11_bind_vector("std::vector<std::vector<Point<3, double>>>", opaque=True, local=True)

from PLC import *
from Tessellation import *
from Tessellator import *
from SerialTessellators import *
