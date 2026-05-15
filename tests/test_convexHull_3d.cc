// Try tessellating a simple lattice of generators in a box in parallel.
// We use randomly chosen seed locations to divide up the generators
// between processors.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "convexHull_3d.hh"
#include "polytope_test_utilities.hh"
#include "simplifyPLCfacets.hh"
#include "ReducedPLC.hh"
#include "polytope_write_OOGL.hh"
#include "Generators.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;

namespace {

//------------------------------------------------------------------------------
// Check if a point (p) is contained in a convex PLC surface.
// Returns the facet id for the first facet the point is outside, or the number
// of facets if it's contained.
//------------------------------------------------------------------------------
template<typename RealType>
unsigned
convexContains(const polytope::PLC<3>& surface,
               const RealType* points,
               const RealType* p,
               const RealType tolerance) {
  bool containmentTest = true;
  unsigned ifacet = 0, n, iv, i, j, k;
  RealType xn, yn, zn, dx1, dy1, dz1, dx2, dy2, dz2, norm;
  while (containmentTest and ifacet < surface.facets.size()) {
    n = surface.facets[ifacet].size();
    POLY_CHECK(n >= 3);
    xn = 0.0; yn = 0.0; zn = 0.0;

    // Compute the facet normal.
    i = surface.facets[ifacet][0];
    for (iv = 1; iv != n - 1; ++iv) {
      j = surface.facets[ifacet][iv];
      k = surface.facets[ifacet][iv + 1];
      dx1 = points[3*j]     - points[3*i];
      dy1 = points[3*j + 1] - points[3*i + 1];
      dz1 = points[3*j + 2] - points[3*i + 2];
      dx2 = points[3*k]     - points[3*i];
      dy2 = points[3*k + 1] - points[3*i + 1];
      dz2 = points[3*k + 2] - points[3*i + 2];
      xn += dy1*dz2 - dz1*dy2;
      yn += dz1*dx2 - dx1*dz2;
      zn += dx1*dy2 - dy1*dx2;
    }
    norm = sqrt(xn*xn + yn*yn + zn*zn);
    xn /= norm;
    yn /= norm;
    zn /= norm;

    // Now test the containment.
    dx1 = p[0] - points[3*i];
    dy1 = p[1] - points[3*i + 1];
    dz1 = p[2] - points[3*i + 2];
    containmentTest = ((dx1*xn + dy1*yn + dz1*zn) < tolerance);
    if (containmentTest) {
      ++ifacet;
    } else {
      cerr << "Facet distance : " << dx1*xn + dy1*yn + dz1*zn << endl
           << "Facet point    : " << points[3*i] << " " << points[3*i + 1] << " " << points[3*i + 2] << endl
           << "Facet normal   : " << xn << " " << yn << " " << zn << endl
           << "Test point     : " << p[0] << " " << p[1] << " " << p[2] << endl;
      return ifacet;
    }
  }

  return ifacet;
}

//------------------------------------------------------------------------------
// Emergency dump.
//------------------------------------------------------------------------------
std::string
escapePod(const std::string nameEnd,
          const polytope::PLC<3>& plc,
          const std::vector<double>& points) {
    std::stringstream os;
    os << "test_PLC_convexHull_" << nameEnd;
    writePLCtoOFF(plc, points, os.str());
    return " : attempted to write to file " + os.str();
}

}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Tessellate a cube.
  {
    cout << "Cube test." << endl;
    double low[3] = {0.0, 0.0, 0.0};
    double high[3] = {1.0, 1.0, 1.0};
    vector<double> points(3*8);
    points[3*0+0] = low[0];  points[3*0+1] = low[1];  points[3*0+2] = high[2];
    points[3*1+0] = high[0]; points[3*1+1] = low[1];  points[3*1+2] = high[2];
    points[3*2+0] = low[0];  points[3*2+1] = high[1]; points[3*2+2] = high[2];
    points[3*3+0] = high[0]; points[3*3+1] = high[1]; points[3*3+2] = high[2];
    points[3*4+0] = low[0];  points[3*4+1] = low[1];  points[3*4+2] = low[2];
    points[3*5+0] = high[0]; points[3*5+1] = low[1];  points[3*5+2] = low[2];
    points[3*6+0] = low[0];  points[3*6+1] = high[1]; points[3*6+2] = low[2];
    points[3*7+0] = high[0]; points[3*7+1] = high[1]; points[3*7+2] = low[2];

    // Get the hull.
    cout << "Generating convex hull... ";
    clock_t t0 = clock();
    const double tolerance = 1.0e-5;
    polytope::PLC<3> box_hull = polytope::convexHull_3d(points, low, tolerance);
    clock_t t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

    // escapePod("box", box_hull, points);
    // polytope::ReducedPLC<3, double> simple_box_hull = simplifyPLCfacets(box_hull,
    //                                                                     points,
    //                                                                     low,
    //                                                                     high,
    //                                                                     1.0e-10);
    // escapePod("simple_box", simple_box_hull, points);
  }

  // Tessellate random points.
  {
    cout << "Generating random points....";
    clock_t t0 = clock();
    const unsigned n = 100000;
    double low[3] = {0.0, 0.0, 0.0};
    double high[3] = {1.0, 1.0, 1.0};
    vector<double> points;
    for (unsigned i = 0; i != n; ++i) {
      points.push_back(low[0] + (high[0] - low[0])*random01());
      points.push_back(low[1] + (high[1] - low[1])*random01());
      points.push_back(low[2] + (high[2] - low[2])*random01());
    }
    POLY_CHECK(points.size() == 3*n);
    clock_t t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

    // Get the hull.
    cout << "Generating convex hull... ";
    t0 = clock();
    const double tolerance = 1.0e-5;
    polytope::PLC<3> hull = polytope::convexHull_3d(points, low, tolerance);
    t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

    // Check that the hull contains all the input vertices.
    cout << "Checking hull for containment... ";
    t0 = clock();
    unsigned badFacet = 0, nFacets = hull.facets.size();
    for (unsigned i = 0; i != n; ++i) {
      badFacet = convexContains(hull, &points.front(), &points[3*i], tolerance);
      if (badFacet < nFacets) cerr << "Failed containment for facet : " << badFacet << endl;
      POLY_CHECK(badFacet == nFacets);
    }
    t1 = clock();
    cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

    cerr << "Initial random hull has (nverts, nfacets) = (" << points.size()/3 <<  ", " << hull.facets.size() << ")" << endl;
    // escapePod("random", hull, points);
    // polytope::ReducedPLC<3, double> simple_random_hull = simplifyPLCfacets(hull,
    //                                                                        points,
    //                                                                        low,
    //                                                                        high,
    //                                                                        1.0e-10);
    // cerr << "Simplified random hull has (nverts, nfacets) = (" << simple_random_hull.points.size()/3 <<  ", " << simple_random_hull.facets.size() << ")" << endl;
    // escapePod("simple_random", simple_random_hull, simple_random_hull.points);
  }

  

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif

  return 0;
}
