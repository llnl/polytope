from PYB11Generator import *

@PYB11template("CoordType")
@PYB11cppname("Point2")
class Point2:
    "A 2D coordinate point with an index."

    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                value="const %(CoordType)s"):
        "Construct with x and y set to one value."

    def pyinit2(self,
                x="const %(CoordType)s",
                y="const %(CoordType)s",
                index=("const unsigned", "0")):
        "Construct from x and y."

    @PYB11const
    def iszero(self):
        return "bool"

    def zero(self):
        return "void"

    def one(self):
        return "void"

    @PYB11const
    def maxAxis(self):
        return "int"

    @PYB11const
    def allLess(self, rhs="const Point2<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def allLessEqual(self, rhs="const Point2<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def allGreater(self, rhs="const Point2<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def allGreaterEqual(self, rhs="const Point2<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def minElements(self, rhs="const Point2<%(CoordType)s>&"):
        return "Point2<%(CoordType)s>"

    @PYB11const
    def maxElements(self, rhs="const Point2<%(CoordType)s>&"):
        return "Point2<%(CoordType)s>"

    @PYB11implementation("[](const Point2<%(CoordType)s>& self) { std::stringstream ss; ss << self; return ss.str(); }")
    def __str__(self):
        return "std::string"

    @PYB11implementation("[](const Point2<%(CoordType)s>& self, size_t i) { if (i >= 2) throw py::index_error(); return self[i]; }")
    def __getitem__(self, i="size_t"):
        return "%(CoordType)s"

    x = PYB11readwrite()
    y = PYB11readwrite()
    index = PYB11readwrite()

@PYB11template("CoordType")
@PYB11cppname("Point3")
class Point3:
    "A 3D coordinate point with an index."

    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                value="const %(CoordType)s"):
        "Construct with x, y, and z set to one value."

    def pyinit2(self,
                x="const %(CoordType)s",
                y="const %(CoordType)s",
                z="const %(CoordType)s",
                index=("const unsigned", "0")):
        "Construct from x, y, and z."

    @PYB11const
    def iszero(self):
        return "bool"

    def zero(self):
        return "void"

    def one(self):
        return "void"

    @PYB11const
    def maxAxis(self):
        return "int"

    @PYB11const
    def allLess(self, rhs="const Point3<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def allLessEqual(self, rhs="const Point3<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def allGreater(self, rhs="const Point3<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def allGreaterEqual(self, rhs="const Point3<%(CoordType)s>&"):
        return "bool"

    @PYB11const
    def minElements(self, rhs="const Point3<%(CoordType)s>&"):
        return "Point3<%(CoordType)s>"

    @PYB11const
    def maxElements(self, rhs="const Point3<%(CoordType)s>&"):
        return "Point3<%(CoordType)s>"

    @PYB11implementation("[](const Point3<%(CoordType)s>& self) { std::stringstream ss; ss << self; return ss.str(); }")
    def __str__(self):
        return "std::string"

    @PYB11implementation("[](const Point3<%(CoordType)s>& self, size_t i) { if (i >= 3) throw py::index_error(); return self[i]; }")
    def __getitem__(self, i="size_t"):
        return "%(CoordType)s"

    x = PYB11readwrite()
    y = PYB11readwrite()
    z = PYB11readwrite()
    index = PYB11readwrite()

for ndim in [2, 3]:
    for (name, ctype) in [("", "double"),
                          ("Coord", f"polytope::HashKey<{ndim}>::IntType"),
                          ("Hash", f"polytope::HashKey<{ndim}>::CoordHash")]:
        exec(f'''
{name}Point{ndim}d = PYB11TemplateClass(Point{ndim}, template_parameters="{ctype}")
''')
