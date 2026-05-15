// test_StarBoundary

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {
  unsigned i;
  int test = 1;
  string testName = "StarBoundary_" + tessellator.name();
  
  // Set the boundary and tessellator
  Boundary2D<double> boundary;
  boundary.setStarWithHole();
  
  // 9 input generators points
  double gens[18] = {0.11,  0.11,
                     0.10,  0.44,
                     0.05,  0.70,
                    -0.40,  0.10,
                    -0.48,  0.20,
                    -0.30, -0.32,
                     0.45, -0.55,
                     0.55,  0.05,
                     0.80,  0.22};

  // Test 1: Input generators, no points on boundary
  {
    cout << "\nTest 1: 9 input generators" << endl;
    std::vector<double> points;
    for (i = 0; i < 9; ++i){
      points.push_back(gens[2*i  ]);
      points.push_back(gens[2*i+1]);
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,testName,points,test);
    ++test;
  }
  

  // Test 2: Input generators + boundary generators
  {
    cout << "\nTest 2: 9 input generators + boundary generators" << endl;
    std::vector<double> points;
    for (i = 0; i < 9; ++i){
      points.push_back(gens[2*i  ]);
      points.push_back(gens[2*i+1]);
    }
    points.insert(points.begin(),
                  boundary.mPLCpoints.begin(),
                  boundary.mPLCpoints.end());
    Tessellation<2,double> mesh;
    //tessellator.tessellate(points,mesh);
    //tessellator.tessellate(points,points,boundary.mPLC,mesh);
    tessellator.tessellate(points,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,testName,points,test);
    ++test;
  }
  

  // Test 3: 800 random generators
  {
    cout << "\nTest 3: 800 random generators" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(800);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,testName,generators.mPoints,test);
    ++test;
  }
  

  // Test 4: 2000 random generators
  {
    cout << "\nTest 4: 2000 random generators" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(2000);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,testName,generators.mPoints,test);
    ++test;
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef POLYTOPE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#else
  POLY_CONTRACT_VAR(argc);
  POLY_CONTRACT_VAR(argv);
#endif


#ifdef POLYTOPE_ENABLE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    test(tessellator);  
  }
#endif   


#ifdef POLYTOPE_ENABLE_BOOST
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    test(tessellator);
  }
#endif
   
  

#ifdef POLYTOPE_ENABLE_MPI
  MPI_Finalize();
#endif
   return 0;
}
