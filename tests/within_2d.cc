// 2D within on PLC unit test.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "within.hh"
#include "polytope_test_utilities.hh"

using namespace std;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Create a PLC with one hole:
  //  
  //  |\_           _/|
  //  |  \_       _/  |
  //  |    \_   _/    |
  //  |      \_/      |
  //  |               |
  //  |     |\_/|     |
  //  |     |___|     |
  //  |_______________|
  //
  const unsigned numVertices = 10;
  double vertices[20] = {-5.0, -5.0,
                          5.0, -5.0,
                          5.0,  5.0,
                          0.0,  0.0,
                         -5.0,  5.0,
                         -1.0, -4.0,
                          1.0, -4.0,
                          1.0, -1.0,
                          0.0, -2.5,
                         -1.0, -1.0};
  polytope::PLC<2, double> plc;
  plc.facets.resize(5);
  plc.holes.resize(1);
  plc.holes[0].resize(5);
  for (unsigned i = 0; i != 5; ++i) {
    plc.facets[i].resize(2);
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i + 1) % 5;
    plc.holes[0][i].resize(2);
    plc.holes[0][i][0] = 5 + i;
    plc.holes[0][i][1] = 5 + (i + 1) % 5;
  }


  {
     double p[2] = {1.0, 1.0};
     bool answer = true;
     int nvert = 16;
     double verts[32] = {0.0, 2.0, 0.0, 1.0,
                         0.0, 0.0, 1.0, 0.0,
                         1.2, 0.0, 1.2, 1.0,
                         1.2, 1.3, 1.3, 1.3,
                         1.3, 1.0, 2.0, 1.0,
                         2.0, 2.0, 1.3, 2.0,
                         1.3, 1.7, 1.2, 1.7,
                         1.2, 2.0, 1.0, 2.0};
     polytope::PLC<2,double> boundary;
     boundary.facets.resize(nvert, std::vector<int>(2));
     for (unsigned i = 0; i != nvert; ++i) {
        boundary.facets[i][0] = i;
        boundary.facets[i][1] = (i+1)%nvert;
     }
     const bool inside = within(p, nvert, verts, boundary);
     POLY_CHECK2(inside == answer, "WTF!?");
  }



  // Test some points.  
  unsigned i = 1;
  { // 1)  Inside: fully
    double p[2] = {-4.0, 1.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2( inside == answer, "Test " << i );  ++i;
  }

  { // 2) Outside: fully
    double p[2] = {-6.0, -3.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 3) Outside: within hole
    double p[2] = {0.0, -3.5};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 4) Inside: tangent to vertex
    double p[2] = {3.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }
  
  { // 5) Outside: tangent to vertex
    double p[2] = {0.0, 5.0};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 6) Inside: tangent to hole vertex
    double p[2] = {2.0, -2.5};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 7) Outside: within hole, tangent to hole vertex
    double p[2] = {-0.5, -2.5};
    bool answer = false;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 8) Inside: tangent to face of hole
    double p[2] = {-1.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 9) Inside: on left outer boundary
    double p[2] = {-5.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 10) Inside: on right outer boundary
    double p[2] = {5.0, 0.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 11) Inside: on bottom outer boundary
    double p[2] = {0.0, -5.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  { // 12) Inside: on hole boundary
    double p[2] = {-1.0, -3.0};
    bool answer = true;
    const bool inside = within(p, numVertices, vertices, plc);
    POLY_CHECK2(inside == answer, "Test " << i );  ++i;
  }

  

  cout << "PASS" << endl;
  return 0;
}
