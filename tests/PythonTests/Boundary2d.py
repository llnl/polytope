"""Python counterpart of the 2-D fixtures in ``Boundary2D.hh``."""

from enum import IntEnum
from math import cos, pi, sin, sqrt
from random import uniform

import polytope


class BoundaryType(IntEnum):
    square, circle, donut, mwithholes, funkystar, circlewithstarhole, cardioid, trogdor, starwithhole, trogdor2, squarewithstarhole, squarewithtrihole = range(12)


class Boundary2d:
    """Construct the standard PLC fixtures and query them through QuantPLC2d."""

    def __init__(self, input_type):
        self.mDiff = 0.5
        self.mCenter = [0.0, 0.0]
        self.m_pad = 0.1
        self.mType = input_type
        self.clear()
        self.setDefaultBoundary(self.mType)

    def clear(self):
        self.PLC = polytope.PLC2d()
        self.PLCpoints, self._facets, self._holes = [], [], []
        self.QPLC = polytope.QuantPLC2d()
        self.low = self.high = None
        self.outerSet = False

    def finalize(self):
        if not self.PLCpoints:
            raise ValueError("Cannot finalize an empty boundary")
        polytope.Quantizer2d.instance().init(self.PLCpoints, self.m_pad)
        self.QPLC.init(self.PLC, self.PLCpoints)
        self.low = [min(self.PLCpoints[0::2]), min(self.PLCpoints[1::2])]
        self.high = [max(self.PLCpoints[0::2]), max(self.PLCpoints[1::2])]

    @staticmethod
    def signedDoubleArea(vertices):
        if len(vertices) < 6 or len(vertices) % 2:
            raise ValueError("A polygon requires at least three (x, y) vertices")
        n = len(vertices)//2
        return sum(vertices[2*i]*vertices[2*((i + 1) % n) + 1] -
                   vertices[2*i + 1]*vertices[2*((i + 1) % n)] for i in range(n))

    area = signedDoubleArea       # Compatibility: this is twice the signed area.

    def setOuterBoundary(self, vertices):
        vertices = list(vertices)
        if self.signedDoubleArea(vertices) <= 0.0:
            raise ValueError("Outer-boundary vertices must be counter-clockwise")
        self.clear()
        self.PLCpoints = vertices
        n = len(vertices)//2
        self._facets = [[i, (i + 1) % n] for i in range(n)]
        self.PLC.facets = self._facets
        self.outerSet = True

    def addHole(self, vertices):
        vertices = list(vertices)
        if not self.outerSet:
            raise ValueError("Set the outer boundary before adding a hole")
        if self.signedDoubleArea(vertices) >= 0.0:
            raise ValueError("Hole vertices must be clockwise")
        start, n = len(self.PLCpoints)//2, len(vertices)//2
        self.PLCpoints.extend(vertices)
        self._holes.append([[start + i, start + (i + 1) % n] for i in range(n)])
        self.PLC.holes = self._holes

    def addCircularHole(self, radius, center, nfacets):
        if radius <= 0.0 or nfacets < 3:
            raise ValueError("A circular hole requires a positive radius and at least three facets")
        vertices = []
        for i in range(nfacets):
            theta = 2.0*pi*(1.0 - i/(nfacets + 1.0))
            vertices += [center[0] + radius*cos(theta), center[1] + radius*sin(theta)]
        self.addHole(vertices)

    def addBoxHole(self, low, high):
        if low[0] >= high[0] or low[1] >= high[1]:
            raise ValueError("Box-hole lower bounds must be less than upper bounds")
        self.addHole([low[0], low[1], low[0], high[1], high[0], high[1], high[0], low[1]])

    def setCustomBoundary(self, numVertices, vertices):
        if numVertices != len(vertices)//2:
            raise ValueError("numVertices does not match the supplied coordinates")
        self.setOuterBoundary(vertices)
        self.finalize()

    def setDefaultBoundary(self, boundaryType):
        factories = (self.setUnitSquare, self.setUnitCircle, self.setDonut, self.setMWithHoles,
                     self.setFunkyStar, self.setCircleWithStarHole, self.setCardioid,
                     self.setTrogdor, self.setStarWithHole, self.setTrogdor2,
                     self.setSquareWithStarHole, self.setSquareWithTriHole)
        factories[BoundaryType(boundaryType)]()

    def _circleVertices(self, radius, nfacets=90):
        result = []
        for i in range(nfacets):
            theta = 2.0*pi*i/(nfacets + 1.0)
            result += [self.mCenter[0] + radius*cos(theta), self.mCenter[1] + radius*sin(theta)]
        return result

    def _starVertices(self, nPoints, outerRadius, clockwise):
        if nPoints < 3:
            raise ValueError("A star requires at least three points")
        theta0 = 2.0*pi/nPoints
        innerRadius = outerRadius*sin(theta0/4.0)/sin(3.0*theta0/4.0)
        result, direction = [], -1.0 if clockwise else 1.0
        for i in range(nPoints):
            theta = pi/2.0 + direction*i*theta0
            result += [self.mCenter[0] + outerRadius*cos(theta), self.mCenter[1] + outerRadius*sin(theta)]
            theta += direction*theta0/2.0
            result += [self.mCenter[0] + innerRadius*cos(theta), self.mCenter[1] + innerRadius*sin(theta)]
        return result

    def _setOuter(self, vertices, boundaryType):
        self.setOuterBoundary(vertices)
        self.mType = boundaryType
        self.finalize()

    def setUnitSquare(self):
        x1, y1 = self.mCenter[0] - self.mDiff, self.mCenter[1] - self.mDiff
        x2, y2 = self.mCenter[0] + self.mDiff, self.mCenter[1] + self.mDiff
        self._setOuter([x1, y1, x2, y1, x2, y2, x1, y2], BoundaryType.square)

    def setUnitCircle(self):
        self._setOuter(self._circleVertices(2.0*self.mDiff), BoundaryType.circle)

    def setDonut(self, innerRadius=0.25):
        if not 0.0 < innerRadius < 1.0:
            raise ValueError("innerRadius must lie strictly between zero and one")
        self.setOuterBoundary(self._circleVertices(2.0*self.mDiff))
        self.addCircularHole(innerRadius, self.mCenter, 90)
        self.mType = BoundaryType.donut
        self.finalize()

    def setMWithHoles(self):
        self.setOuterBoundary([0., 0., 2., 0., 2., 2., 1., 1., 0., 2.])
        self.addBoxHole([.25, .25], [.75, .75])
        self.addBoxHole([1.25, .25], [1.75, .75])
        self.mType = BoundaryType.mwithholes
        self.finalize()

    def setFunkyStar(self):
        vertices = []
        for i in range(9, -1, -1):
            radius = 2.0 + .5*sin(12.0*pi*i/9.0)
            vertices += [radius*sin(pi*i/9.0), radius*cos(pi*i/9.0)]
        self._setOuter(vertices, BoundaryType.funkystar)

    def setCircleWithStarHole(self, nPoints=5):
        self.setOuterBoundary(self._circleVertices(2.0*self.mDiff))
        self.addHole(self._starVertices(nPoints, .75, clockwise=True))
        self.mType = BoundaryType.circlewithstarhole
        self.finalize()

    def setCardioid(self, z=2.0):
        if z <= 0.0:
            raise ValueError("Cardioid coefficient must be positive")
        vertices = []
        for i in range(90):
            theta = 2.0*pi*i/91.0
            vertices += [z*cos(theta) - cos(2.0*theta), z*sin(theta) - sin(2.0*theta)]
        self._setOuter(vertices, BoundaryType.cardioid)

    @staticmethod
    def _reverseVertices(vertices):
        return [coordinate for point in reversed(list(zip(vertices[::2], vertices[1::2])))
                for coordinate in point]

    def setTrogdor(self):
        points = [2., 9., 4., 8.9, 5., 9.2, 6.5, 8.8, 7., 8., 6.5, 7., 5., 6.3,
                  4., 5.5, 3.7, 4.6, 4., 3.5, 5.5, 2.8, 6.8, 3.4, 5.5, 2.5,
                  4., 2.6, 3., 3., 2.5, 4., 2.4, 4.5, 2.5, 5.3, 3., 5.9, 4.5,
                  7., 4.9, 7.3, 5.1, 7.7, 5., 8., 4.5, 8., 3.5, 7.5, 2., 7.,
                  2., 7.4, 3.3, 8., 1.75, 8., 1.75, 9.2]
        self._setOuter(self._reverseVertices(points), BoundaryType.trogdor)

    def setStarWithHole(self):
        self.setOuterBoundary(self._starVertices(5, 1.0, clockwise=False))
        self.addHole([.05, -.05, .10, .10, .20, -.30, -.25, -.15])
        self.mType = BoundaryType.starwithhole
        self.finalize()

    def setTrogdor2(self):
        points = [5.2, 2., 7., 1.5, 7.2, 2.7, 8., 1.2, 8.5, 1.2, 8.5, 3., 9.5, 1.6,
                  10.4, 2., 9.8, 3.8, 11.4, 3.3, 11.9, 5.6, 11., 7., 9.8, 7.8, 7., 9.,
                  5.6, 10.5, 4.9, 9.6, 3.4, 8.8, 4., 8.2, 4.1, 6.5, 5., 6.6, 5.2, 6.1,
                  5., 6., 5.1, 5.5, 5.4, 5.5, 5.3, 5., 4., 5., 3., 6., 1.5, 7., 1.1, 8.,
                  .7, 8.7, 1.3, 9.6, .8, 10.2, .9, 11.7, 2., 12.5, 3.2, 12.3, 3.8, 13.,
                  4.6, 13.1, 5.2, 14., 4., 14.3, 3., 14., 2.8, 15., 1.9, 15.3, 1.3, 16.2,
                  .5, 16.6, 2.2, 17.5, 3.6, 15.9, 5.5, 14.3, 6.5, 14.8, 7.5, 17.2, 7.9,
                  19., 9.8, 18.5, 9.7, 17.8, 9.8, 17., 9.2, 16.8, 8.8, 16., 8., 15.9,
                  7.2, 15.1, 9.6, 15., 12.2, 14., 14.5, 14., 14.9, 14.8, 15.5, 14.8,
                  16., 14.2, 15.8, 12.7, 12.5, 12.7, 13.1, 11.8, 15.5, 11., 15.4, 10.5,
                  13.5, 10.4, 12., 11.3, 9.9, 12.1, 8.8, 12.3, 8.7, 12., 9., 11., 11., 9.2,
                  14., 7., 14.8, 5., 14., 3., 12., 1.1, 10.2, .5, 8.7, .4, 7.2, .8]
        self._setOuter(self._reverseVertices(points), BoundaryType.trogdor2)

    def setSquareWithStarHole(self, nPoints=5):
        self.mDiff = 1.0
        self.setOuterBoundary([self.mCenter[0] - 1., self.mCenter[1] - 1., self.mCenter[0] + 1., self.mCenter[1] - 1., self.mCenter[0] + 1., self.mCenter[1] + 1., self.mCenter[0] - 1., self.mCenter[1] + 1.])
        self.addHole(self._starVertices(nPoints, .75, clockwise=True))
        self.mType = BoundaryType.squarewithstarhole
        self.finalize()

    def setSquareWithTriHole(self):
        self.mDiff = 1.0
        self.setOuterBoundary([self.mCenter[0] - 1., self.mCenter[1] - 1., self.mCenter[0] + 1., self.mCenter[1] - 1., self.mCenter[0] + 1., self.mCenter[1] + 1., self.mCenter[0] - 1., self.mCenter[1] + 1.])
        self.addHole([.6, -.8, .4, -.8, .4, .8])
        self.mType = BoundaryType.squarewithtrihole
        self.finalize()

    def testInside(self, x, y):
        return self.QPLC.within((x, y))

    def getBoundingRadius(self):
        return max(sqrt((x - self.mCenter[0])**2 + (y - self.mCenter[1])**2)
                   for x, y in zip(self.PLCpoints[::2], self.PLCpoints[1::2]))

    def getPointInside(self):
        if self.low is None:
            raise ValueError("Finalize the boundary before sampling it")
        while True:
            point = uniform(self.low[0], self.high[0]), uniform(self.low[1], self.high[1])
            if self.testInside(*point):
                return point

    # Legacy construction spellings retained for existing Python callers.
    def initBox(self, low=(0., 0.), high=(1., 1.)):
        self.mCenter = [.5*(low[0] + high[0]), .5*(low[1] + high[1])]
        self._setOuter([low[0], low[1], high[0], low[1], high[0], high[1], low[0], high[1]], BoundaryType.square)

    def initCircle(self, radius=1., center=(0., 0.), nfacets=90):
        self.mCenter = list(center)
        self._setOuter(self._circleVertices(radius, nfacets), BoundaryType.circle)

    def initTorus(self, outerRadius=1., innerRadius=.25, center=(0., 0.), nfacets=90):
        self.mCenter = list(center)
        self.setOuterBoundary(self._circleVertices(outerRadius, nfacets))
        self.addCircularHole(innerRadius, center, nfacets)
        self.mType = BoundaryType.donut
        self.finalize()

    initMWithHoles = setMWithHoles
    initCardioid = setCardioid
    initTrogdor = setTrogdor2


Boundary2D = Boundary2d
