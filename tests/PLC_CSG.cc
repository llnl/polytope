// Unit tests for Computational Solid Geometry (CSG) operations on PLC's.

#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>

#include "polytope.hh"
#include "PLC_CSG_2d.hh"
#include "PLC_CSG_3d.hh"
#include "simplifyPLCfacets.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "polytope_test_utilities.hh"
#include "polytope_write_OOGL.hh"
#include "polytope_plc_canned_geometries.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

namespace {

//------------------------------------------------------------------------------
// Rotate a set of points (3D).
//------------------------------------------------------------------------------
template<typename RealType>
void
rotatePoints(std::vector<RealType>& points,
             const RealType theta,
             const RealType phi,
             const RealType psi) {
  POLY_ASSERT(points.size() % 3 == 0);
  const RealType R[3][3] = {{ cos(theta)*cos(psi), cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), sin(phi)*sin(psi) - cos(phi)*sin(theta)*cos(psi)},
                            {-cos(theta)*sin(psi), cos(phi)*cos(psi) - sin(phi)*sin(theta)*sin(psi), sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)},
                            { sin(theta),         -sin(phi)*cos(theta),                              cos(phi)*cos(theta)}};
  const unsigned n = points.size()/3;
  for (unsigned i = 0; i != n; ++i) {
    const RealType x = points[3*i], y = points[3*i+1], z = points[3*i+2];
    points[3*i  ] = R[0][0]*x + R[0][1]*y + R[0][2]*z;
    points[3*i+1] = R[1][0]*x + R[1][1]*y + R[1][2]*z;
    points[3*i+2] = R[2][0]*x + R[2][1]*y + R[2][2]*z;
  }
}

//------------------------------------------------------------------------------
// Emergency dump.
//------------------------------------------------------------------------------
template<typename PLCType>
std::string
escapePod(const std::string nameEnd,
          const PLCType& plc) {
    std::stringstream os;
    os << "test_PLC_CSG_" << nameEnd;
    writePLCtoOFF(plc, plc.points, os.str());
    return " : attempted to write to file " + os.str();
}

} // anonymous namespace

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
  typedef Point3<double> PointType;
  typedef Point3<int64_t> IntPointType;

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  //----------------------------------------------------------------------
  // Convert from PLC->CSG_internal_2d::Segments->PLC
  //----------------------------------------------------------------------
  {
    cout << "2D PLC <-> Segments test." << endl;
    double low[2] = {0.0, 0.0}, high[2] = {1.0, 1.0};
    const ReducedPLC<2, double> box1 = plc_box<2, double>(low, high);
    const std::vector<CSG::CSG_internal_2d::Segment<double> > segments = CSG::CSG_internal_2d::ReducedPLCtoSegments(box1);
    POLY_CHECK2(segments.size() == 4, "Num segments = " << segments.size() << ", expected 4.");
    const ReducedPLC<2, double> box2 = CSG::CSG_internal_2d::ReducedPLCfromSegments(segments);
    POLY_CHECK(box2.facets.size() == 4);
    POLY_CHECK(box2.points.size() == 2*8);
    const ReducedPLC<2, double> box3 = polytope::simplifyPLCfacets(box2, box2.points, low, high, 1.0e-10);
    POLY_CHECK(box3.facets.size() == 4);
    POLY_CHECK(box3.points.size() == 2*4);
  }

  //----------------------------------------------------------------------
  // Same as above, but force the creation of collinear edges and make
  // sure they're removed in the end.
  // Convert from PLC->CSG_internal_2d::Segments->PLC
  //----------------------------------------------------------------------
  {
    cout << "2D PLC <-> Segments test with collinear points." << endl;
    double low[2] = {0.0, 0.0}, high[2] = {1.0, 1.0};
    ReducedPLC<2, double> box1;
    box1.points.push_back(low[0]);                 box1.points.push_back(low[1]);
    box1.points.push_back(0.5*(low[0] + high[0])); box1.points.push_back(low[1]);
    box1.points.push_back(high[0]);                box1.points.push_back(low[1]);
    box1.points.push_back(high[0]);                box1.points.push_back(0.5*(low[1] + high[1]));
    box1.points.push_back(high[0]);                box1.points.push_back(high[1]);
    box1.points.push_back(0.5*(low[0] + high[0])); box1.points.push_back(high[1]);
    box1.points.push_back(low[0]);                 box1.points.push_back(high[1]);
    box1.points.push_back(low[0]);                 box1.points.push_back(0.5*(low[1] + high[1]));
    box1.facets.resize(8);
    box1.facets[0].push_back(0); box1.facets[0].push_back(1);
    box1.facets[1].push_back(1); box1.facets[1].push_back(2);
    box1.facets[2].push_back(2); box1.facets[2].push_back(3);
    box1.facets[3].push_back(3); box1.facets[3].push_back(4);
    box1.facets[4].push_back(4); box1.facets[4].push_back(5);
    box1.facets[5].push_back(5); box1.facets[5].push_back(6);
    box1.facets[6].push_back(6); box1.facets[6].push_back(7);
    box1.facets[7].push_back(7); box1.facets[7].push_back(0);
    const std::vector<CSG::CSG_internal_2d::Segment<double> > segments = CSG::CSG_internal_2d::ReducedPLCtoSegments(box1);
    POLY_CHECK2(segments.size() == 8, "Num segments = " << segments.size() << ", expected 8.");
    const ReducedPLC<2, double> box2 = CSG::CSG_internal_2d::ReducedPLCfromSegments(segments);
    POLY_CHECK(box2.facets.size() == 8);
    POLY_CHECK(box2.points.size() == 2*16);
    const ReducedPLC<2, double> box3 = polytope::simplifyPLCfacets(box2, box2.points, low, high, 1.0e-10);
    POLY_CHECK(box3.facets.size() == 4);
    POLY_CHECK(box3.points.size() == 2*4);
  }

  //----------------------------------------------------------------------
  // Operate on two boxes offset in x.
  //----------------------------------------------------------------------
  {
    cout << "2D box interactions test." << endl;
    double low1[2] = {0.0, 0.0}, high1[2] = {1.0, 1.0},
           low2[2] = {0.5, 0.0}, high2[2] = {1.5, 1.0};
    const ReducedPLC<2, double> box1 = plc_box<2, double>(low1, high1),
                                box2 = plc_box<2, double>(low2, high2);
    const ReducedPLC<2, double> box_union = CSG::csg_union(box1, box2);
    POLY_CHECK(box_union.points.size() % 2 == 0);
    const ReducedPLC<2, double> box_union_simplify = polytope::simplifyPLCfacets(box_union, box_union.points, low1, high2, 1.0e-10);
    const ReducedPLC<2, double> box_intersect = CSG::csg_intersect(box1, box2);
    const ReducedPLC<2, double> box_intersect_simplify = polytope::simplifyPLCfacets(box_intersect, box_intersect.points, low1, high2, 1.0e-10);
    POLY_CHECK(box_intersect_simplify.facets.size() == 4);
    POLY_CHECK(box_intersect_simplify.points.size() == 2*4);
  }

  //----------------------------------------------------------------------
  // Various combinations of a box and circle.
  //----------------------------------------------------------------------
  {
    cout << "2D box-circle interactions test." << endl;
    const Point2<double> low(0.0, 0.0), high(1.0, 1.0), high2(2.0, 2.0);
    const ReducedPLC<2, double> square = plc_box<2, double>(&low.x, &high.x),
                                circle = plc_circle<double>(high, 0.5, 20);
    const ReducedPLC<2, double> sc_union = CSG::csg_union(square, circle),
                       sc_union_simplify = polytope::simplifyPLCfacets(sc_union, sc_union.points, &low.x, &high2.x, 1.0e-10);
    // escapePod("square_circle_union", sc_union);
    // escapePod("square_circle_union_simplify", sc_union_simplify);
    const ReducedPLC<2, double> sc_intersect = CSG::csg_intersect(square, circle),
                       sc_intersect_simplify = polytope::simplifyPLCfacets(sc_intersect, sc_intersect.points, &low.x, &high2.x, 1.0e-10);
    // escapePod("square_circle_intersect", sc_intersect);
    // escapePod("square_circle_intersect_simplify", sc_intersect_simplify);
    const ReducedPLC<2, double> sc_subtract = CSG::csg_subtract(square, circle),
                       sc_subtract_simplify = polytope::simplifyPLCfacets(sc_subtract, sc_subtract.points, &low.x, &high2.x, 1.0e-10);
    // escapePod("square_circle_subtract", sc_subtract);
    // escapePod("square_circle_subtract_simplify", sc_subtract_simplify);
  }

  //----------------------------------------------------------------------
  // Various combinations of a box and circle (int64 versions).
  //----------------------------------------------------------------------
  {
    cout << "2D box-circle interactions test (int64_t based)." << endl;
    const Point2<int64_t> low(0, 0), high(1 << 10, 1 << 10), high2(1 << 11, 1 << 11);
    const ReducedPLC<2, int64_t> square = plc_box<2, int64_t>(&low.x, &high.x),
                                 circle = plc_circle<int64_t>(high, (1 << 9), 20);
    const ReducedPLC<2, int64_t> sc_union = CSG::csg_union(square, circle),
                        sc_union_simplify = polytope::simplifyPLCfacets(sc_union, sc_union.points, &low.x, &high2.x, int64_t(1));
    // escapePod("square_circle_union", sc_union);
    // escapePod("square_circle_union_simplify", sc_union_simplify);
    const ReducedPLC<2, int64_t> sc_intersect = CSG::csg_intersect(square, circle),
                       sc_intersect_simplify = polytope::simplifyPLCfacets(sc_intersect, sc_intersect.points, &low.x, &high2.x, int64_t(1));
    // escapePod("square_circle_intersect", sc_intersect);
    // escapePod("square_circle_intersect_simplify", sc_intersect_simplify);
    const ReducedPLC<2, int64_t> sc_subtract = CSG::csg_subtract(square, circle),
                       sc_subtract_simplify = polytope::simplifyPLCfacets(sc_subtract, sc_subtract.points, &low.x, &high2.x, int64_t(1));
    // escapePod("square_circle_subtract", sc_subtract);
    // escapePod("square_circle_subtract_simplify", sc_subtract_simplify);
  }

  //----------------------------------------------------------------------
  // Convert from PLC->CSG_internal_3d::Polygons->PLC
  //----------------------------------------------------------------------
  {
    cout << "3D PLC <-> Polygons test." << endl;
    typedef Point3<double> RealPoint;
    double low[3] = {0.0, 0.0, 0.0}, high[3] = {1.0, 1.0, 1.0};
    const ReducedPLC<3, double> box1 = plc_box<3, double>(low, high);
    std::vector<CSG::CSG_internal_3d::Polygon<double> > polys = CSG::CSG_internal_3d::ReducedPLCtoPolygons(box1);
    POLY_CHECK2(polys.size() == 12, "Num polys = " << polys.size() << ", expected 12.");

    // Deliberately split one triangle facet into three, creating some inconsistent hanging nodes.
    std::vector<CSG::CSG_internal_3d::Vertex<double> > triangle(3);
    polys.pop_back(); polys.pop_back();  // Pop off the last two triangles.
    POLY_ASSERT(box1.facets.back().size() == 4);
    const unsigned k1 = box1.facets.back()[0],
                   k2 = box1.facets.back()[1],
                   k3 = box1.facets.back()[2],
                   k4 = box1.facets.back()[3];
    triangle[0] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k1], box1.points[3*k1+1], box1.points[3*k1+2]));
    triangle[1] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k2], box1.points[3*k2+1], box1.points[3*k2+2]));
    triangle[2] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k3], box1.points[3*k3+1], box1.points[3*k3+2]));
    polys.push_back(triangle);
    RealPoint diag1((box1.points[3*k1  ] + 2.0*box1.points[3*k3  ])/3.0,
                    (box1.points[3*k1+1] + 2.0*box1.points[3*k3+1])/3.0,
                    (box1.points[3*k1+2] + 2.0*box1.points[3*k3+2])/3.0),
              diag2((2.0*box1.points[3*k1  ] + box1.points[3*k3  ])/3.0,
                    (2.0*box1.points[3*k1+1] + box1.points[3*k3+1])/3.0,
                    (2.0*box1.points[3*k1+2] + box1.points[3*k3+2])/3.0);
    triangle[0] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k3], box1.points[3*k3+1], box1.points[3*k3+2]));
    triangle[1] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k4], box1.points[3*k4+1], box1.points[3*k4+2]));
    triangle[2] = CSG::CSG_internal_3d::Vertex<double>(diag1);
    polys.push_back(triangle);
    triangle[0] = CSG::CSG_internal_3d::Vertex<double>(diag1);
    triangle[1] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k4], box1.points[3*k4+1], box1.points[3*k4+2]));
    triangle[2] = CSG::CSG_internal_3d::Vertex<double>(diag2);
    polys.push_back(triangle);
    triangle[0] = CSG::CSG_internal_3d::Vertex<double>(diag2);
    triangle[1] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k4], box1.points[3*k4+1], box1.points[3*k4+2]));
    triangle[2] = CSG::CSG_internal_3d::Vertex<double>(RealPoint(box1.points[3*k1], box1.points[3*k1+1], box1.points[3*k1+2]));
    polys.push_back(triangle);

    const ReducedPLC<3, double> box2 = CSG::CSG_internal_3d::ReducedPLCfromPolygons(polys);
    std::cerr << "box2 : " << box2 << std::endl;
    POLY_CHECK(box2.facets.size() == 14);
    const ReducedPLC<3, double> box3 = polytope::simplifyPLCfacets(box2, box2.points, low, high, 1.0e-10);
    POLY_CHECK(box3.facets.size() == 6);
    POLY_CHECK(box3.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Operate on two boxes offset in x.
  //----------------------------------------------------------------------
  {
    cout << "3D box interactions test." << endl;
    double low1[3] = {0.0, 0.0, 0.0}, high1[3] = {1.0, 1.0, 1.0},
           low2[3] = {0.5, 0.0, 0.0}, high2[3] = {1.5, 1.0, 1.0};
    const ReducedPLC<3, double> box1 = plc_box<3, double>(low1, high1),
                                box2 = plc_box<3, double>(low2, high2);
    const ReducedPLC<3, double> box_union = CSG::csg_union(box1, box2);
    POLY_CHECK(box_union.points.size() % 3 == 0);
    escapePod("box1", box1);
    escapePod("box2", box2);
    cerr << escapePod("box_union_test", box_union) << endl;
    const ReducedPLC<3, double> box_union_simplify = polytope::simplifyPLCfacets(box_union, box_union.points, low1, high2, 1.0e-10);
    cerr << escapePod("box_union_test_simplify", box_union_simplify) << endl;
    const ReducedPLC<3, double> box_intersect = CSG::csg_intersect(box1, box2);
    // cerr << escapePod("box_intersect_test", box_intersect) << endl;
    const ReducedPLC<3, double> box_intersect_simplify = polytope::simplifyPLCfacets(box_intersect, box_intersect.points, low1, high2, 1.0e-10);
    // cerr << escapePod("box_intersect_test_simplify", box_intersect_simplify) << endl;
    POLY_CHECK(box_intersect_simplify.facets.size() == 6);
    POLY_CHECK(box_intersect_simplify.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Generate a complicated sphere with holes cut out of it.
  //----------------------------------------------------------------------
  {
    cout << "3D sphere-cylinder interactions test." << endl;
    const PointType origin(0.0, 0.0, 0.0);
    const unsigned nphi = 20;
    clock_t t0 = clock();
    double low[3] = {-1.0, -1.0, -1.0}, high[3] = {1.0, 1.0, 1.0};
    const ReducedPLC<3, double> a = plc_box<3, double>(low, high),
                                b = plc_sphere<double>(origin, 1.35, nphi),
                                c = plc_cylinder(origin, 0.7, 2.0, nphi);
    ReducedPLC<3, double> d(c), e(c);
    rotatePoints(d.points, 0.0, M_PI/2.0, 0.0);
    rotatePoints(e.points, 0.0, M_PI/2.0, M_PI/2.0);
    clock_t t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds to generate input PLCs." << endl;
    t0 = clock();
    const ReducedPLC<3, double> holey_sphere = CSG::csg_subtract(CSG::csg_intersect(a, b), 
                                                                 CSG::csg_union(CSG::csg_union(c, d), e));
    t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds to perform CSG operations." << endl;
    escapePod("holey_sphere", holey_sphere);
  }

  // ********************************************************************************
  // Repeat above tests for int types.
  //----------------------------------------------------------------------
  // Convert from PLC->CSG_internal_3d::Polygons->PLC
  //----------------------------------------------------------------------
  {
    cout << "3D PLC <-> Polygons test (in64_t)." << endl;
    int64_t low[3] = {0, 0, 0}, high[3] = {1 << 10, 1 << 10, 1 << 10};
    const ReducedPLC<3, int64_t> box1 = plc_box<3, int64_t>(low, high);
    const std::vector<CSG::CSG_internal_3d::Polygon<int64_t> > polys = CSG::CSG_internal_3d::ReducedPLCtoPolygons(box1);
    POLY_CHECK2(polys.size() == 12, "Num polys = " << polys.size() << ", expected 12.");
    const ReducedPLC<3, int64_t> box2 = CSG::CSG_internal_3d::ReducedPLCfromPolygons(polys);
    const ReducedPLC<3, int64_t> box3 = polytope::simplifyPLCfacets(box2, box2.points, low, high, int64_t(1));
    POLY_CHECK(box3.facets.size() == 6);
    POLY_CHECK(box3.points.size() == 3*8);
  }

  //----------------------------------------------------------------------
  // Operate on two boxes offset in x.
  //----------------------------------------------------------------------
  {
    cout << "3D box interactions test (int64_t)." << endl;
    int64_t low1[3] = {0, 0, 0},          high1[3] = {1 << 10, 1 << 10, 1 << 10},
            low2[3] = {high1[0]/2, 0, 0}, high2[3] = {low2[0] + high1[0], high1[1], high1[2]};
    const ReducedPLC<3, int64_t> box1 = plc_box<3, int64_t>(low1, high1),
                                 box2 = plc_box<3, int64_t>(low2, high2);
    const ReducedPLC<3, int64_t> box_union = CSG::csg_union(box1, box2);
    POLY_CHECK(box_union.points.size() % 3 == 0);
    escapePod("box1_int", box1);
    escapePod("box2_int", box2);
    cerr << escapePod("box_union_test_int", box_union) << endl;
    const ReducedPLC<3, int64_t> box_union_simplify = polytope::simplifyPLCfacets(box_union, box_union.points, low1, high2, int64_t(1));
    cerr << escapePod("box_union_test_simplify_int", box_union_simplify) << endl;
    const ReducedPLC<3, int64_t> box_intersect = CSG::csg_intersect(box1, box2);
    cerr << escapePod("box_intersect_test_int", box_intersect) << endl;
    const ReducedPLC<3, int64_t> box_intersect_simplify = polytope::simplifyPLCfacets(box_intersect, box_intersect.points, low1, high2, int64_t(1));
    cerr << escapePod("box_intersect_test_simplify_int", box_intersect_simplify) << endl;
    POLY_CHECK(box_intersect_simplify.facets.size() == 6);
    POLY_CHECK(box_intersect_simplify.points.size() == 3*8);
  }

  cout << "PASS" << endl;

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif
  return 0;
}
