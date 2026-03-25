// test_MeshEditor
//
// Unit tests for the MeshEditor class and cleanEdges routine

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "MeshEditor.hh"
#include "polytope_test_utilities.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#include "checkDistributedTessellation.hh"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef POLYTOPE_ENABLE_MPI
   MPI_Init(&argc, &argv);
#endif

   string testName = "MeshEditor";
   int itest = 1;

   // Test 1: Single cell, one small face
   {
      const double edgeTol = 0.001;
      cout << "\nTest 1: Eliminate small face on a single 2D cell" << endl;
      Tessellation<2,double> mesh;
      const unsigned ncells = 1;
      const unsigned nfaces = 6;
      const unsigned nnodes = 6;
      mesh.cells.resize(ncells, std::vector<int>(nfaces));
      mesh.faces.resize(nfaces, std::vector<unsigned>(2));
      mesh.faceCells.resize(nfaces, std::vector<int>(1));
      for (int i = 0; i != nfaces; ++i) {
         mesh.cells[0][i] = i;
         mesh.faces[i][0] = i;
         mesh.faces[i][1] = (i+1) % nnodes;
         mesh.faceCells[i][0] = 0;
      }
      double nodePoints[12] = {0.00, 0.0,
                               0.75, 0.0,
                               1.00, 0.5,
                               1.00, 0.500001,
                               0.75, 1.0,
                               0.00, 1.0};
      mesh.nodes.resize(2*nnodes);
      for (int i = 0; i != 2*nnodes; ++i) mesh.nodes[i] = nodePoints[i];
      MeshEditor<2,double> meshEditor(mesh);
      cout << mesh << endl;
      outputMesh(mesh, testName, itest);
      meshEditor.cleanEdges(edgeTol);
      cout << mesh << endl;
      outputMesh(mesh, testName, itest);
      ++itest;
   }

   // Test 2: Distributed test: edge is small on one proc, but not on the other
#ifdef POLYTOPE_ENABLE_MPI
   {
     const double edgeTol = 0.001;
     cout << "\nTest 2: Edge is small on one proc, but not on the other" << endl;
      int rank, numProcs;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
     if (numProcs == 2) {
       Tessellation<2,double> mesh;
       unsigned ncells, nfaces, nnodes;
       unsigned otherProc, sharedNode1, sharedNode2, sharedFace;
       std::vector<double> nodePoints;
       if (rank == 0) {
         ncells = 1, nfaces = 6, nnodes = 6;
         otherProc = 1, sharedNode1 = 2, sharedNode2 = 3, sharedFace = 2;       
         nodePoints.push_back(0.00);  nodePoints.push_back(0.0);
         nodePoints.push_back(0.75);  nodePoints.push_back(0.0);
         nodePoints.push_back(1.00);  nodePoints.push_back(0.5);
         nodePoints.push_back(1.00);  nodePoints.push_back(0.5000001);
         nodePoints.push_back(0.75);  nodePoints.push_back(1.0);
         nodePoints.push_back(0.00);  nodePoints.push_back(1.0);
       } else {
         ncells = 1, nfaces = 4, nnodes = 4;
         otherProc = 0, sharedNode1 = 1, sharedNode2 = 0, sharedFace = 0;
         nodePoints.push_back(1.00     );  nodePoints.push_back(0.5000001);
         nodePoints.push_back(1.00     );  nodePoints.push_back(0.5       );
         nodePoints.push_back(1.0000001);  nodePoints.push_back(0.5       );
         nodePoints.push_back(1.0000001);  nodePoints.push_back(0.5000001);
       }
       
       // Set the mesh data
       mesh.cells.resize(ncells, std::vector<int>(nfaces));
       mesh.faces.resize(nfaces, std::vector<unsigned>(2));
       mesh.faceCells.resize(nfaces, std::vector<int>(1));
       mesh.neighborDomains.resize(1);
       mesh.sharedNodes.resize(1, std::vector<unsigned>(2));
       mesh.sharedFaces.resize(1, std::vector<unsigned>(1));
       mesh.nodes.resize(2*nnodes);
       for (int i = 0; i != nfaces; ++i) {
         mesh.cells[0][i] = i;
         mesh.faces[i][0] = i;
         mesh.faces[i][1] = (i+1) % nnodes;
         mesh.faceCells[i][0] = 0;
       }
       for (int i = 0; i != 2*nnodes; ++i) mesh.nodes[i] = nodePoints[i];
       mesh.neighborDomains[0] = otherProc;
       mesh.sharedNodes[0][0]  = sharedNode1;
       mesh.sharedNodes[0][1]  = sharedNode2;
       mesh.sharedFaces[0][0]  = sharedFace;
       
       // Flipped-orientation fix-up for the shared face on proc 1
       if (rank == 1) {
         mesh.cells[0][0] = ~0;
         mesh.faces[0][0] = 1;
         mesh.faces[0][1] = 0;
         mesh.faceCells[0][0] = ~0;
       }
       MeshEditor<2,double> meshEditor(mesh);
       cout << mesh << endl;
       outputMesh(mesh, testName, itest);
       meshEditor.cleanEdges(edgeTol);
       cout << mesh << endl;
       outputMesh(mesh, testName, itest);
       ++itest;
#ifdef POLYTOPE_ENABLE_MPI
       const string parCheck = checkDistributedTessellation(mesh);
       POLY_CHECK2(parCheck == "ok", parCheck);
#endif
     }
   }
#endif

   cout << "PASS" << endl;

#ifdef POLYTOPE_ENABLE_MPI
   MPI_Finalize();
#endif
   return 0;
}
