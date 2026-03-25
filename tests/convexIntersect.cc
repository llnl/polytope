// 2D nearestPoint on PLC unit test.

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <stdlib.h>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "convexIntersect.hh"
#include "convexHull_2d.hh"

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {
   
  // Test 1
  {
    const double tol = 1.0e-10;
    const unsigned Na = 200;
    const unsigned Nb = 200;
    srand(103020);

    // Generate random points for the inner plc
    const double alow [2] = {-0.5, -0.5};
    const double ahigh[2] = { 0.5,  0.5};
    vector<double> apoints;
    for (unsigned i = 0; i < Na; ++i) {
      apoints.push_back(alow[0] + (ahigh[0]-alow[0])*random01());
      apoints.push_back(alow[1] + (ahigh[1]-alow[1])*random01());
    }
  
    // Compute its convex hull
    ReducedPLC<2, double> ahull(convexHull_2d(apoints, alow, tol), apoints);

    // Generate random points for the outer plc
    const double blow [2] = {-1.0, -1.0};
    const double bhigh[2] = { 1.0,  1.0};
    vector<double> bpoints;
    double x,y;
    unsigned i = 0;
    while (i < Nb) {
      x = blow[0] + (bhigh[0]-blow[0])*random01();
      y = blow[1] + (bhigh[1]-blow[1])*random01();
      if (not (x > alow[0] and x < ahigh[0] and
               y > alow[1] and y < ahigh[1]) ) {
        bpoints.push_back(x);
        bpoints.push_back(y);
        ++i;
      }
    }
    
    // Compute its convex hull
    ReducedPLC<2, double> bhull(convexHull_2d(bpoints, blow, tol), bpoints);    
    
    cerr << "Test 1: Hull A is inside Hull B...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(ahull, bhull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection" << endl;
  }


  // Test 2
  {
    ReducedPLC<2,double> ahull;
    double apts[6] = {-0.228435,0.0223897,-0.371005,0.335035,0.351504,0.389491};
    for (unsigned i = 0; i != 6; ++i) ahull.points.push_back(apts[i]);
    ahull.facets.resize(3, vector<int>(2));
    ahull.facets[0][0] = 1;  ahull.facets[0][1] = 0;
    ahull.facets[1][0] = 0;  ahull.facets[1][1] = 2;
    ahull.facets[2][0] = 2;  ahull.facets[2][1] = 1;
    ReducedPLC<2,double> bhull;
    double bpts[8] = {-0.58538,-0.73863,-0.950275,-0.447837,0.971177,-0.596531,-0.477071,0.927211};
    for (unsigned i = 0; i != 8; ++i) bhull.points.push_back(bpts[i]);
    bhull.facets.resize(4, vector<int>(2));
    bhull.facets[0][0] = 1;  bhull.facets[0][1] = 0;
    bhull.facets[1][0] = 0;  bhull.facets[1][1] = 2;
    bhull.facets[2][0] = 2;  bhull.facets[2][1] = 3;
    bhull.facets[3][0] = 3;  bhull.facets[3][1] = 1;
    
    cerr << "Test 2: Hull A has points inside Hull B...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(ahull, bhull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection" << endl;
  }

  // Test 3
  {
    double apts[8] = {0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0};
    double bpts[8] = {1.0, 1.0, 2.0,-0.5, 3.0, 1.0, 2.0, 2.5};
    ReducedPLC<2,double> ahull, bhull;
    for (unsigned i = 0; i != 8; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(4, vector<int>(2));
    bhull.facets.resize(4, vector<int>(2));
    for (unsigned i = 0; i != 4; ++i) {
      ahull.facets[i][0] = i;  ahull.facets[i][1] = (i+1)%4;
      bhull.facets[i][0] = i;  bhull.facets[i][1] = (i+1)%4;
    }

    cerr << "Test 3: Hull B has a point inside Hull A...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(ahull, bhull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection" << endl;
  }

  // Test 4
  {
    double apts[6] = {0.0, 1.0, 4.0, 1.0, 2.0, 4.0};
    double bpts[6] = {2.0, 0.0, 4.0, 3.0, 0.0, 3.0};
    ReducedPLC<2,double> ahull, bhull;
    for (unsigned i = 0; i != 6; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(3, vector<int>(2));
    bhull.facets.resize(3, vector<int>(2));
    for (unsigned i = 0; i != 3; ++i) {
      ahull.facets[i][0] = i;  ahull.facets[i][1] = (i+1)%3;
      bhull.facets[i][0] = i;  bhull.facets[i][1] = (i+1)%3;
    }

    cerr << "Test 4: Intersection, but neither hull has vertices inside the other...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(ahull, bhull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection" << endl;
  }
  

  // Test 5
  {
    double apts[6] = {0.0, 4.0, 4.0, 5.0, 2.0, 8.0};
    double bpts[6] = {2.0, 0.0, 4.0, 3.0, 0.0, 3.0};
    ReducedPLC<2,double> ahull, bhull;
    for (unsigned i = 0; i != 6; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(3, vector<int>(2));
    bhull.facets.resize(3, vector<int>(2));
    for (unsigned i = 0; i != 3; ++i) {
      ahull.facets[i][0] = i;  ahull.facets[i][1] = (i+1)%3;
      bhull.facets[i][0] = i;  bhull.facets[i][1] = (i+1)%3;
    }

    cerr << "Test 5: Disjoint hulls...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(ahull, bhull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == false       );
    cerr << "No Intersection" << endl;
  }

  cout << "PASS" << endl;
  return 0;
}
