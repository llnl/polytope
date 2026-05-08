// 3D convex polyhedron intersection unit test.

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <stdlib.h>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "convexIntersect.hh"
#include "convexHull_3d.hh"

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Test 1
  {
    const double tol = 1.0e-10;
    const unsigned Na = 100;
    const unsigned Nb = 100;
    srand(103020);

    // Generate random points for the inner polyhedron
    const double alow [3] = {-0.5, -0.5, -0.5};
    const double ahigh[3] = { 0.5,  0.5,  0.5};
    vector<double> apoints;
    for (unsigned i = 0; i < Na; ++i) {
      apoints.push_back(alow[0] + (ahigh[0]-alow[0])*random01());
      apoints.push_back(alow[1] + (ahigh[1]-alow[1])*random01());
      apoints.push_back(alow[2] + (ahigh[2]-alow[2])*random01());
    }

    // Compute its convex hull
    ReducedPLC<3, double> ahull(convexHull_3d(apoints, alow, tol), apoints);

    // Generate random points for the outer polyhedron
    const double blow [3] = {-1.0, -1.0, -1.0};
    const double bhigh[3] = { 1.0,  1.0,  1.0};
    vector<double> bpoints;
    double x, y, z;
    unsigned i = 0;
    while (i < Nb) {
      x = blow[0] + (bhigh[0]-blow[0])*random01();
      y = blow[1] + (bhigh[1]-blow[1])*random01();
      z = blow[2] + (bhigh[2]-blow[2])*random01();
      if (not (x > alow[0] and x < ahigh[0] and
               y > alow[1] and y < ahigh[1] and
               z > alow[2] and z < ahigh[2]) ) {
        bpoints.push_back(x);
        bpoints.push_back(y);
        bpoints.push_back(z);
        ++i;
      }
    }

    // Compute its convex hull
    ReducedPLC<3, double> bhull(convexHull_3d(bpoints, blow, tol), bpoints);

    cerr << "Test 1: Hull A is inside Hull B...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection" << endl;
  }

  // Test 2: Two tetrahedra with partial overlap
  {
    // Tetrahedron A: vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    double apts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0};
    // Tetrahedron B: shifted to overlap partially A
    // Tetrahedron C: shifted to not overlap A
    double shift = 0.5;
    ReducedPLC<3,double> ahull, bhull, chull;
    for (unsigned i = 0; i < 12; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(apts[i] + shift - 0.2);
      chull.points.push_back(apts[i] + shift);
    }
    // Tetrahedron has 4 triangular faces
    int num_tris = 4;
    ahull.facets.resize(num_tris, vector<int>(3));
    bhull.facets.resize(num_tris, vector<int>(3));
    chull.facets.resize(num_tris, vector<int>(3));
    vector<vector<int>> facets = {{0, 2, 1},
                                  {0, 1, 3},
                                  {0, 3, 2},
                                  {1, 2, 3}};
    for (auto i = 0; i < num_tris; ++i) {
      ahull.facets[i] = facets[i];
      bhull.facets[i] = facets[i];
      chull.facets[i] = facets[i];
    }
    cerr << "Test 2: Overlapping tetrahedra...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    bool aIntersectsc = convexIntersect(ahull, chull);
    bool cIntersectsa = convexIntersect(chull, ahull);
    POLY_CHECK(aIntersectsc == cIntersectsa);
    POLY_CHECK(aIntersectsc == false       );
    cerr << "Intersection" << endl;
  }

  // Test 3: Two cubes with one vertex of B inside A
  {
    // Cube A: unit cube at origin
    double apts[24] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       1.0, 1.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0,
                       1.0, 0.0, 1.0,
                       1.0, 1.0, 1.0,
                       0.0, 1.0, 1.0};
    // Cube B: overlapping corner
    double bpts[24] = {0.5, 0.5, 0.5,
                       1.5, 0.5, 0.5,
                       1.5, 1.5, 0.5,
                       0.5, 1.5, 0.5,
                       0.5, 0.5, 1.5,
                       1.5, 0.5, 1.5,
                       1.5, 1.5, 1.5,
                       0.5, 1.5, 1.5};
    ReducedPLC<3,double> ahull, bhull;
    for (unsigned i = 0; i != 24; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    // Cube has 6 square faces
    ahull.facets.resize(6, vector<int>(4));
    bhull.facets.resize(6, vector<int>(4));
    // Bottom face (z=0)
    ahull.facets[0][0] = 0; ahull.facets[0][1] = 1; ahull.facets[0][2] = 2; ahull.facets[0][3] = 3;
    // Top face (z=1)
    ahull.facets[1][0] = 4; ahull.facets[1][1] = 7; ahull.facets[1][2] = 6; ahull.facets[1][3] = 5;
    // Front face (y=0)
    ahull.facets[2][0] = 0; ahull.facets[2][1] = 4; ahull.facets[2][2] = 5; ahull.facets[2][3] = 1;
    // Back face (y=1)
    ahull.facets[3][0] = 2; ahull.facets[3][1] = 6; ahull.facets[3][2] = 7; ahull.facets[3][3] = 3;
    // Left face (x=0)
    ahull.facets[4][0] = 0; ahull.facets[4][1] = 3; ahull.facets[4][2] = 7; ahull.facets[4][3] = 4;
    // Right face (x=1)
    ahull.facets[5][0] = 1; ahull.facets[5][1] = 5; ahull.facets[5][2] = 6; ahull.facets[5][3] = 2;
    // Same for B
    bhull.facets[0][0] = 0; bhull.facets[0][1] = 1; bhull.facets[0][2] = 2; bhull.facets[0][3] = 3;
    bhull.facets[1][0] = 4; bhull.facets[1][1] = 7; bhull.facets[1][2] = 6; bhull.facets[1][3] = 5;
    bhull.facets[2][0] = 0; bhull.facets[2][1] = 4; bhull.facets[2][2] = 5; bhull.facets[2][3] = 1;
    bhull.facets[3][0] = 2; bhull.facets[3][1] = 6; bhull.facets[3][2] = 7; bhull.facets[3][3] = 3;
    bhull.facets[4][0] = 0; bhull.facets[4][1] = 3; bhull.facets[4][2] = 7; bhull.facets[4][3] = 4;
    bhull.facets[5][0] = 1; bhull.facets[5][1] = 5; bhull.facets[5][2] = 6; bhull.facets[5][3] = 2;

    cerr << "Test 3: Cube B has a vertex inside Cube A...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection" << endl;
  }

  // Test 4: Disjoint tetrahedra
  {
    // Tetrahedron A at origin
    double apts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0};
    // Tetrahedron B far away
    double bpts[12] = {5.0, 5.0, 5.0,
                       6.0, 5.0, 5.0,
                       5.0, 6.0, 5.0,
                       5.0, 5.0, 6.0};
    ReducedPLC<3,double> ahull, bhull;
    for (unsigned i = 0; i != 12; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(4, vector<int>(3));
    bhull.facets.resize(4, vector<int>(3));
    ahull.facets[0][0] = 0; ahull.facets[0][1] = 2; ahull.facets[0][2] = 1;
    ahull.facets[1][0] = 0; ahull.facets[1][1] = 1; ahull.facets[1][2] = 3;
    ahull.facets[2][0] = 0; ahull.facets[2][1] = 3; ahull.facets[2][2] = 2;
    ahull.facets[3][0] = 1; ahull.facets[3][1] = 2; ahull.facets[3][2] = 3;
    bhull.facets[0][0] = 0; bhull.facets[0][1] = 2; bhull.facets[0][2] = 1;
    bhull.facets[1][0] = 0; bhull.facets[1][1] = 1; bhull.facets[1][2] = 3;
    bhull.facets[2][0] = 0; bhull.facets[2][1] = 3; bhull.facets[2][2] = 2;
    bhull.facets[3][0] = 1; bhull.facets[3][1] = 2; bhull.facets[3][2] = 3;

    cerr << "Test 4: Disjoint tetrahedra...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == false       );
    cerr << "No Intersection" << endl;
  }

  // Test 5: Two tetrahedra touching at exactly one vertex
  {
    // Tetrahedron A at origin
    double apts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0};
    // Tetrahedron B touching at vertex (1,0,0)
    double bpts[12] = {1.0, 0.0, 0.0,
                       2.0, 0.0, 0.0,
                       1.0, 1.0, 0.0,
                       1.0, 0.0, 1.0};
    ReducedPLC<3,double> ahull, bhull;
    for (unsigned i = 0; i != 12; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(4, vector<int>(3));
    bhull.facets.resize(4, vector<int>(3));
    ahull.facets[0][0] = 0; ahull.facets[0][1] = 2; ahull.facets[0][2] = 1;
    ahull.facets[1][0] = 0; ahull.facets[1][1] = 1; ahull.facets[1][2] = 3;
    ahull.facets[2][0] = 0; ahull.facets[2][1] = 3; ahull.facets[2][2] = 2;
    ahull.facets[3][0] = 1; ahull.facets[3][1] = 2; ahull.facets[3][2] = 3;
    bhull.facets[0][0] = 0; bhull.facets[0][1] = 2; bhull.facets[0][2] = 1;
    bhull.facets[1][0] = 0; bhull.facets[1][1] = 1; bhull.facets[1][2] = 3;
    bhull.facets[2][0] = 0; bhull.facets[2][1] = 3; bhull.facets[2][2] = 2;
    bhull.facets[3][0] = 1; bhull.facets[3][1] = 2; bhull.facets[3][2] = 3;

    cerr << "Test 5: Touching at exactly one vertex...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );  // Touching counts as intersection
    cerr << "Intersection (touching)" << endl;
  }

  // Test 6: Two tetrahedra sharing a common edge
  {
    // Tetrahedron A
    double apts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.5, 1.0, 0.0,
                       0.5, 0.5, 1.0};
    // Tetrahedron B sharing edge from (0,0,0) to (1,0,0)
    double bpts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.5, -1.0, 0.0,
                       0.5, -0.5, 1.0};
    ReducedPLC<3,double> ahull, bhull;
    for (unsigned i = 0; i != 12; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(4, vector<int>(3));
    bhull.facets.resize(4, vector<int>(3));
    ahull.facets[0][0] = 0; ahull.facets[0][1] = 2; ahull.facets[0][2] = 1;
    ahull.facets[1][0] = 0; ahull.facets[1][1] = 1; ahull.facets[1][2] = 3;
    ahull.facets[2][0] = 0; ahull.facets[2][1] = 3; ahull.facets[2][2] = 2;
    ahull.facets[3][0] = 1; ahull.facets[3][1] = 2; ahull.facets[3][2] = 3;
    bhull.facets[0][0] = 0; bhull.facets[0][1] = 2; bhull.facets[0][2] = 1;
    bhull.facets[1][0] = 0; bhull.facets[1][1] = 1; bhull.facets[1][2] = 3;
    bhull.facets[2][0] = 0; bhull.facets[2][1] = 3; bhull.facets[2][2] = 2;
    bhull.facets[3][0] = 1; bhull.facets[3][1] = 2; bhull.facets[3][2] = 3;

    cerr << "Test 6: Sharing a common edge...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );  // Sharing edge counts as intersection
    cerr << "Intersection (shared edge)" << endl;
  }

  // Test 7: Two tetrahedra sharing a common face
  {
    // Tetrahedron A
    double apts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.5, 1.0, 0.0,
                       0.5, 0.5, 1.0};
    // Tetrahedron B sharing the base triangle with A
    double bpts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.5, 1.0, 0.0,
                       0.5, 0.5, -1.0};
    ReducedPLC<3,double> ahull, bhull;
    for (unsigned i = 0; i != 12; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(bpts[i]);
    }
    ahull.facets.resize(4, vector<int>(3));
    bhull.facets.resize(4, vector<int>(3));
    ahull.facets[0][0] = 0; ahull.facets[0][1] = 2; ahull.facets[0][2] = 1;
    ahull.facets[1][0] = 0; ahull.facets[1][1] = 1; ahull.facets[1][2] = 3;
    ahull.facets[2][0] = 0; ahull.facets[2][1] = 3; ahull.facets[2][2] = 2;
    ahull.facets[3][0] = 1; ahull.facets[3][1] = 2; ahull.facets[3][2] = 3;
    bhull.facets[0][0] = 0; bhull.facets[0][1] = 2; bhull.facets[0][2] = 1;
    bhull.facets[1][0] = 0; bhull.facets[1][1] = 1; bhull.facets[1][2] = 3;
    bhull.facets[2][0] = 0; bhull.facets[2][1] = 3; bhull.facets[2][2] = 2;
    bhull.facets[3][0] = 1; bhull.facets[3][1] = 2; bhull.facets[3][2] = 3;

    cerr << "Test 7: Sharing a common face...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );  // Sharing face counts as intersection
    cerr << "Intersection (shared face)" << endl;
  }

  // Test 8: Identical tetrahedra
  {
    double apts[12] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0};
    ReducedPLC<3,double> ahull, bhull;
    for (unsigned i = 0; i != 12; ++i) {
      ahull.points.push_back(apts[i]);
      bhull.points.push_back(apts[i]);
    }
    ahull.facets.resize(4, vector<int>(3));
    bhull.facets.resize(4, vector<int>(3));
    ahull.facets[0][0] = 0; ahull.facets[0][1] = 2; ahull.facets[0][2] = 1;
    ahull.facets[1][0] = 0; ahull.facets[1][1] = 1; ahull.facets[1][2] = 3;
    ahull.facets[2][0] = 0; ahull.facets[2][1] = 3; ahull.facets[2][2] = 2;
    ahull.facets[3][0] = 1; ahull.facets[3][1] = 2; ahull.facets[3][2] = 3;
    bhull.facets[0][0] = 0; bhull.facets[0][1] = 2; bhull.facets[0][2] = 1;
    bhull.facets[1][0] = 0; bhull.facets[1][1] = 1; bhull.facets[1][2] = 3;
    bhull.facets[2][0] = 0; bhull.facets[2][1] = 3; bhull.facets[2][2] = 2;
    bhull.facets[3][0] = 1; bhull.facets[3][1] = 2; bhull.facets[3][2] = 3;

    cerr << "Test 8: Identical tetrahedra...";
    bool aIntersectsb = convexIntersect(ahull, bhull);
    bool bIntersectsa = convexIntersect(bhull, ahull);
    POLY_CHECK(aIntersectsb == bIntersectsa);
    POLY_CHECK(aIntersectsb == true        );
    cerr << "Intersection (identical)" << endl;
  }

  return 0;
}
