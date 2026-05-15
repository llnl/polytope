// 2D nearestPoint on PLC unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "nearestPoint.hh"
#include "polytope_test_utilities.hh"

using namespace std;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Create a square PLC with one hole.
  const unsigned numVertices = 8;
  double vertices[16] = {-5.0, -5.0,
                          5.0, -5.0,
                          5.0,  5.0,
                         -5.0,  5.0,
                         -1.0, -1.0,
                          1.0, -1.0,
                          1.0,  1.0,
                         -1.0,  1.0};
  polytope::PLC<2> plc;
  plc.facets.resize(4);
  plc.holes.resize(1);
  plc.holes[0].resize(4);
  for (unsigned i = 0; i != 4; ++i) {
    plc.facets[i].resize(2);
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 4;
    plc.holes[0][i].resize(2);
    plc.holes[0][i][0] = 4 + i;
    plc.holes[0][i][1] = 4 + (i + 1) % 4;
  }

  // Test some points.
  const double tol = 1.0e-10;

  {
    double p[2] = {-10.0, -10.0}, answer[2] = {-5.0, -5.0}, result[2];
    const double dist = nearestPoint(p, numVertices, vertices, plc, result);
    POLY_CHECK2((polytope::geometry::distance<2, double>(result, answer) < tol),
                result[0] << " " << result[1] << " : " << dist);
  }

  {
    double p[2] = {-10.0, 0.0}, answer[2] = {-5.0, 0.0}, result[2];
    const double dist = nearestPoint(p, numVertices, vertices, plc, result);
    POLY_CHECK2((polytope::geometry::distance<2, double>(result, answer) < tol),
                result[0] << " " << result[1] << " : " << dist);
  }

  {
    double p[2] = {1.0, 0.25}, answer[2] = {1.0, 0.25}, result[2];
    const double dist = nearestPoint(p, numVertices, vertices, plc, result);
    POLY_CHECK2((polytope::geometry::distance<2, double>(result, answer) < tol),
                result[0] << " " << result[1] << " : " << dist);
  }

  // Make sure the new point is EXACTLY on the boundary
  { 
    double p[2] = {-4.5, -4.2127}, answer[2] = {-5.0, -4.2127}, result[2];
    const double dist = nearestPoint(p, numVertices, vertices, plc, result);
    POLY_CHECK2(result[0] == answer[0] and result[1] == answer[1],
                result[0] << " " << result[1] << " : " << dist);
  }
  
  { 
    double p[2] = {-4.2127, -4.775}, answer[2] = {-4.2127, -5.0}, result[2];
    const double dist = nearestPoint(p, numVertices, vertices, plc, result);
    POLY_CHECK2(result[0] == answer[0] and result[1] == answer[1],
                result[0] << " " << result[1] << " : " << dist);
  }
  
  { 
    double p[2]      = {4.888392983934, -4.33233323};
    double answer[2] = {5.000000000000, -4.33233323}, result[2];
    const double dist = nearestPoint(p, numVertices, vertices, plc, result);
    POLY_CHECK2(result[0] == answer[0] and result[1] == answer[1],
                result[0] << " " << result[1] << " : " << dist);
  }

  {
    double result[2];
    for (unsigned i = 0; i != 100; ++i) {
      const double xrand = random01() - 5;
      const double yrand = random01() - 4;
      double p[2] = {xrand, yrand};
      const double dist = nearestPoint(p, numVertices, vertices, plc, result);
      POLY_CHECK2(result[0] == -5.0,
                  p[0] << " " << p[1] << " --> " << result[0] << " " << result[1]);
    }
  }

  {
    double p[2] = {0.5, 0.5}, s1[2] = {1.1, 1.0}, s2[2] = {2.0, 2.0}, result[2];
    polytope::geometry::closestPointOnSegment2D(p, s1, s2, result);
    POLY_CHECK(result[0] == s1[0]);
    POLY_CHECK(result[1] == s1[1]);
  }


  
  return 0;
}
