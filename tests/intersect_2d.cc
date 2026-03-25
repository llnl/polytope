// 2D nearestPoint on PLC unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "intersect.hh"
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
  polytope::PLC<2, double> plc;
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
  
  { // 1 intersection with outside boundary
    double p1[2] = {-10.0, 0.0}, p2[2] = {-4.0, 0.0};
    vector<double> result;
    const unsigned nints = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK(nints == 1 and result.size() == 2);
    POLY_CHECK(result[0] == -5.0 and result[1] == 0.0);
  }

  { // 2 intersections with outside boundary
    double p1[2] = {-6.0, -3.0}, p2[2] = {-3.0, -6.0};
    vector<double> result;
    const unsigned nints = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK(nints == 2 and result.size() == 4);
    POLY_CHECK(result[0] == -4.0 and result[1] == -5.0 and
	       result[2] == -5.0 and result[3] == -4.0);
  }

  { // 1 intersections with outside boundary corner
    double p1[2] = {4.0, -6.0}, p2[2] = {6.0, -4.0};
    vector<double> result;
    const unsigned nints = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK(nints == 1 and result.size() == 2);
    POLY_CHECK(result[0] == 5.0 and result[1] == -5.0);
  }

  { // 2 intersections: outside boundary and hole
    double p1[2] = {0.0, -6.0}, p2[2] = {0.0, 0.0};
    vector<double> result;
    const unsigned nints = intersect(p1, p2, numVertices, vertices, plc, result);
    POLY_CHECK(nints == 2 and result.size() == 4);
    POLY_CHECK(result[0] == 0.0 and result[1] == -5.0 and
	       result[2] == 0.0 and result[3] == -1.0);
  }

  
  cout << "PASS" << endl;
  return 0;
}
