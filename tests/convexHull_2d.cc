// 2D convex hull unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "convexHull_2d.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;

//------------------------------------------------------------------------------
// Check if a point (p) is contained in a convex PLC surface.
// Returns the facet id for the first facet the point is outside, or the number
// of facets if it's contained.
//------------------------------------------------------------------------------
template<typename RealType>
unsigned
convexContains(const polytope::PLC<2, RealType>& surface,
               const RealType* points,
               const RealType* p,
               const RealType tolerance) {
  bool containmentTest = true;
  unsigned ifacet = 0, n, i, j;
  RealType xn, yn, dx1, dy1, norm;
  while (containmentTest and ifacet < surface.facets.size()) {
    n = surface.facets[ifacet].size();
    POLY_CHECK(n == 2);
    xn = 0.0; yn = 0.0;

    // Compute the facet normal.
    i = surface.facets[ifacet][0];
    j = surface.facets[ifacet][1];
    dx1 = points[2*j    ] - points[2*i];
    dy1 = points[2*j + 1] - points[2*i + 1];
    xn = dy1;
    yn = -dx1;
    norm = sqrt(xn*xn + yn*yn);
    xn /= norm;
    yn /= norm;

    // Now test the containment.
    dx1 = p[0] - points[2*i];
    dy1 = p[1] - points[2*i + 1];
    containmentTest = ((dx1*xn + dy1*yn) < tolerance);
    if (containmentTest) {
      ++ifacet;
    } else {
      cerr << "Facet distance : " << dx1*xn + dy1*yn << endl
           << "Facet point    : " << points[2*i] << " " << points[2*i + 1] << endl
           << "Facet normal   : " << xn << " " << yn << endl
           << "Test point     : " << p[0] << " " << p[1] << endl;
    }
  }

  return ifacet;
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Generate some random seed points.
  cout << "Generating random points....";
  clock_t t0 = clock();
  const unsigned n = 100000;
  double low[2] = {0.0, 0.0};
  double high[2] = {1.0, 1.0};
  vector<double> points;
  for (unsigned i = 0; i != n; ++i) {
    points.push_back(low[0] + (high[0] - low[0])*random01());
    points.push_back(low[1] + (high[1] - low[1])*random01());
  }
  POLY_CHECK(points.size() == 2*n);
  clock_t t1 = clock();
  cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

  // Get the hull.
  cout << "Generating convex hull... ";
  t0 = clock();
  const double tolerance = 1.0e-9;
  polytope::PLC<2, double> hull = polytope::convexHull_2d(points, low, tolerance);
  t1 = clock();
  cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

  // Check that the hull contains all the input vertices.
  cout << "Checking hull for containment... ";
  t0 = clock();
  unsigned badFacet = 0, nFacets = hull.facets.size();
  for (unsigned i = 0; i != n; ++i) {
    badFacet = convexContains(hull, &points.front(), &points[2*i], tolerance);
    if (badFacet < nFacets) cerr << "Failed containment for facet : " << badFacet << endl;
    POLY_CHECK(badFacet == nFacets);
  }
  t1 = clock();
  cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

  // for (unsigned k = 0; k != hull.facets.size(); ++k) {
  //   cerr << "Facet " << k << " : ";
  //   for (unsigned j = 0; j != hull.facets[k].size(); ++j) cerr << " " << hull.facets[k][j];
  //   cerr << " : ";
  //   for (unsigned j = 0; j != hull.facets[k].size(); ++j) cerr << " (" << points[2*hull.facets[k][j]] << " " << points[2*hull.facets[k][j] + 1] << ")";
  //   cerr << endl;
  // }

  

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif

  return 0;
}
