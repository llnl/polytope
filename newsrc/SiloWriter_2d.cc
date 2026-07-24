#include "polytope.hh"
#include "SiloWriter.hh"
#include "SiloUtils.hh"
#include "Tessellation.hh"

#ifdef POLYTOPE_ENABLE_MPI
#include "Communicator.hh"
#endif

#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"

namespace polytope {

using namespace std;

namespace {

//-------------------------------------------------------------------
// Traverse the nodes of cell i within the given tessellation in
// order, writing their indices to nodes.  We rely here on two
// assumptions:
// 1.  cellFaces are given such that the faces are in counter-clockwise
//     order around the cell.
// 2.  if cellFaces[j] > 0, the nodes of the face cellFaces[j] are
//     given in counter-clockwise orientation for cell i,
//     otherwise the nodes of ~cellFaces[j] (the 1s complement)
//     are in *clockwise* order and need to be reversed.
//-------------------------------------------------------------------
template <typename RealType>
void
traverseNodes(const Tessellation<2, RealType>& mesh,
              int i,
              vector<int>& nodes) {
  const vector<int>& cellFaces = mesh.cells[i];
  for (size_t j = 0; j != cellFaces.size(); ++j) {
    int k = cellFaces[j];
    nodes.push_back(k >= 0 ? mesh.faces[ k][0] :
                    mesh.faces[~k][1]);
  }
  nodes.push_back(nodes.front());

#ifdef POLYTOPE_ENABLE_DEBUG
  // Make sure we don't have any garbage in our list of nodes.
  for (size_t n = 0; n < nodes.size(); ++n) {
    POLY_ASSERT(nodes[n] >= 0);
    POLY_ASSERT(nodes[n] < mesh.nodes.size()/2);
  }
#endif
}
//-------------------------------------------------------------------
}

//-------------------------------------------------------------------
template <typename RealType>
void
SiloWriter<2, RealType>::write(const Tessellation<2, RealType>& mesh,
                               const map<string, RealType*>& nodeFields,
                               const map<string, vector<int>*>& nodeTags,
                               const map<string, RealType*>& edgeFields,
                               const map<string, vector<int>*>& edgeTags,
                               const map<string, RealType*>& faceFields,
                               const map<string, vector<int>*>& faceTags,
                               const map<string, RealType*>& cellFields,
                               const map<string, vector<int>*>& cellTags,
                               const string& filePrefix,
                               const string& directory,
                               int cycle,
                               RealType time,
                               int numFiles,
                               int mpiTag) {
  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  int index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);
  int coord_sys = DB_CARTESIAN;
  // Open a file in Silo/HDF5 format for writing.
#ifdef POLYTOPE_ENABLE_MPI
  int nproc = Communicator::getNProcs();
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  auto& comm = Communicator::communicator();
  int localRankHasPoints = (mesh.points.size() > 0) ? 1 : 0;
  std::vector<int> ranksWithData = gatherValidRanks(localRankHasPoints);
  if (numFiles == -1) {
    int globalWriteProcs = static_cast<int>(ranksWithData.size());
    MPI_Bcast(&globalWriteProcs, 1, MPI_INT, root, comm);
    numFiles = globalWriteProcs;
  }
  POLY_ASSERT(numFiles <= nproc);

  std::string masterDirName = getMasterDirName(directory, prefix, cycle);
  if (rank == root) {
    DIR* masterDir = opendir(masterDirName.c_str());
    if (masterDir == 0) {
      mkdir(masterDirName.c_str(), S_IRWXU | S_IRWXG);
    } else {
      closedir(masterDir);
    }
  }
  Communicator::Barrier();
  bool hasPoints = bool(localRankHasPoints);

  std::string filename = getFilename(masterDirName, rank);
  std::string meshname = getLocalMeshName();

#else
  string dirname = directory;
  if (dirname.empty()) dirname = ".";
  std::string filename = getMasterFilename(prefix, cycle);

  std::string meshname = getGlobalMeshName();
  bool hasPoints = true;
#endif
  if (hasPoints) {
    DBfile* file = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, 0, DB_HDF5);
    // Add cycle/time metadata if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    DBAddOption(optlist, DBOPT_COORDSYS, &coord_sys);
    // This is optional for now, but we'll give it anyway.
    char *coordnames[2];
    coordnames[0] = (char*)"xcoords";
    coordnames[1] = (char*)"ycoords";

    // Node coordinates.
    int numNodes = mesh.nodes.size() / 2;
    vector<double> x(numNodes), y(numNodes);
    for (int i = 0; i < numNodes; ++i) {
      x[i] = mesh.nodes[2*i];
      y[i] = mesh.nodes[2*i+1];
    }
    double* coords[2];
    coords[0] = &(x[0]);
    coords[1] = &(y[0]);

    // Build the list of nodes describing the boundary faces.
    // int numBoundaryFaces = mesh.boundaryFaces.size();
    // vector<int> boundaryNodes(2*numBoundaryFaces);
    // for (int i = 0; i < numBoundaryFaces; ++i) {
    //   boundaryNodes[2*i] = mesh.faces[mesh.boundaryFaces[i]][0];
    //   boundaryNodes[2*i+1] = mesh.faces[mesh.boundaryFaces[i]][1];
    // }

    // Write the boundary face list.
    // {
    //   vector<int> shapesize(size_t(1), 2), shapecnt(size_t(1), numBoundaryFaces);
    //   DBPutFacelist(file, (char*)"boundary_faces", numBoundaryFaces,
    //                 2, &boundaryNodes[0], boundaryNodes.size(), 0,
    //                 0, &shapesize[0], &shapecnt[0], 1, 0, 0, 0);
    // }

    // All zones are polygonal.
    int numCells = mesh.cells.size();
    vector<int> shapesize(numCells, 0),
      shapetype(numCells, DB_ZONETYPE_POLYGON),
      shapecount(numCells, 1),
      nodeList;
    for (int i = 0; i < numCells; ++i) {
      // Gather the nodes from this cell in traversal order.
      vector<int> cellNodes;
      traverseNodes(mesh, i, cellNodes);
      // Insert the cell's node connectivity into the node list.
      nodeList.push_back(cellNodes.size());
      nodeList.insert(nodeList.end(), cellNodes.begin(), cellNodes.end());
    }

    // Write out the 2D polygonal mesh.
    DBPutUcdmesh(file, meshname.c_str(), 2, coordnames, coords,
                 numNodes, numCells,
                 "mesh_zonelist", NULL, DB_DOUBLE, optlist);
    DBPutZonelist2(file, "mesh_zonelist", numCells,
                   2, &nodeList[0], nodeList.size(), 0, 0, 0,
                   &shapetype[0], &shapesize[0], &shapecount[0],
                   numCells, optlist);

    // Write out the cell-face connectivity data.
    vector<int> conn(numCells);
    int elemlengths[2];
    char* elemnames[2];
    elemnames[0] = strDup("ncellfaces");
    elemlengths[0] = numCells;    
    for (int c = 0; c < numCells; ++c) {
      conn[c] = mesh.cells[c].size();
    }
    for (int c = 0; c < numCells; ++c) {
      for (size_t f = 0; f < mesh.cells[c].size(); ++f) {
        conn.push_back(mesh.cells[c][f]);
      }
    }
    // Size of conn that is the cells
    int connCellSize = static_cast<int>(conn.size()) - numCells;
    elemnames[1] = strDup("cellfaces");
    elemlengths[1] = connCellSize;
    DBPutCompoundarray(file, "conn", elemnames, elemlengths, 2,
                       (void*)&conn[0], conn.size(), DB_INT, 0);
    free(elemnames[0]);
    free(elemnames[1]);

    // Write out convex hull data.
    vector<int> hull(1+mesh.convexHull.facets.size());
    hull[0] = mesh.convexHull.facets.size();
    for (size_t f = 0; f < mesh.convexHull.facets.size(); ++f)
      hull[1+f] = mesh.convexHull.facets[f].size();
    for (size_t f = 0; f < mesh.convexHull.facets.size(); ++f)
      for (size_t n = 0; n < mesh.convexHull.facets[f].size(); ++n)
        hull.push_back(mesh.convexHull.facets[f][n]);
    elemnames[0] = strDup("nfacets");
    elemlengths[0] = 1;
    elemnames[1] = strDup("nfacetnodes");
    elemlengths[1] = mesh.convexHull.facets.size();
    elemnames[2] = strDup("facetnodes");
    elemlengths[2] = hull.size() - elemlengths[0] - elemlengths[1];
    DBPutCompoundarray(file, "convexhull", elemnames, elemlengths, 3,
                       (void*)&hull[0], hull.size(), DB_INT, 0);
    free(elemnames[0]);
    free(elemnames[1]);
    free(elemnames[2]);
    // Write out tag information.
    //writeTagsToFile(nodeTags, file, DB_NODECENT);
    writeTagsToFile(edgeTags, file, DB_EDGECENT);
    writeTagsToFile(faceTags, file, DB_FACECENT);
    writeTagsToFile(cellTags, file, DB_ZONECENT);

    // FIXME: We really should try to use the number of edges for edge fields.
    const int numFaces = mesh.faces.size();
    //writeFieldsToFile<RealType>(nodeFields, varmeshname, file, numNodes, DB_NODECENT, optlist);

    // Write cell-centered fields to CELLS directory
    writeFieldsToFile<RealType>(edgeFields, meshname, file, numFaces, DB_EDGECENT, optlist);
    writeFieldsToFile<RealType>(faceFields, meshname, file, numFaces, DB_FACECENT, optlist);
    writeFieldsToFile<RealType>(cellFields, meshname, file, numCells, DB_ZONECENT, optlist);

    int numPoints = mesh.points.size() / 2;
    vector<double> xp(numPoints), yp(numPoints);
    for (int i = 0; i < numPoints; ++i) {
      xp[i] = mesh.points[2*i];
      yp[i] = mesh.points[2*i+1];
    }
    double* pcoords[2];
    pcoords[0] = &(xp[0]);
    pcoords[1] = &(yp[0]);
    // Write point mesh
    DBPutPointmesh(file, (char*)"points", 2, pcoords, numPoints, DB_DOUBLE, optlist);
#ifdef POLYTOPE_ENABLE_DEBUG
    // Create NODES directory and write the nodes as points
    // Node coordinates
    // Write point mesh of nodes
    DBPutPointmesh(file, (char*)"nodes", 2, coords, numNodes, DB_DOUBLE, optlist);
#endif

    // Clean up.
    DBClose(file);
    DBFreeOptlist(optlist);
  }

#ifdef POLYTOPE_ENABLE_MPI

  // Finally, write the uber-master file.
  if (rank == root) {
    std::string masterFilename = getMasterFilename(prefix, cycle);
    int driver = DB_HDF5;
    DBfile* file = DBCreate(masterFilename.c_str(), DB_CLOBBER, DB_LOCAL, "Master file", driver);

    // Build mesh names
    std::vector<std::string> procPaths = getProcPaths(masterDirName, ranksWithData);
    int nblocks = static_cast<int>(ranksWithData.size());
    std::vector<char*> cellMeshNames;
    std::vector<char*> pointMeshNames;
    std::vector<int> varTypes(nblocks, DB_UCDVAR);
    std::vector<int> cellMeshTypes(nblocks, DB_UCDMESH);
    std::vector<int> pointMeshTypes(nblocks, DB_POINTMESH);
    for (const auto& p : procPaths) {
      cellMeshNames.push_back(strDup((p + meshname).c_str()));
      pointMeshNames.push_back(strDup((p + "points").c_str()));
    }
    DBoptlist* masteroptlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(masteroptlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(masteroptlist, DBOPT_DTIME, &dtime);
    std::string global_mesh_name = getGlobalMeshName();
    DBAddOption(masteroptlist, DBOPT_MMESH_NAME, global_mesh_name.data());
    DBAddOption(masteroptlist, DBOPT_COORDSYS, &coord_sys);

    DBPutMultimesh(file, global_mesh_name.c_str(), nblocks, cellMeshNames.data(), cellMeshTypes.data(), masteroptlist);
    DBPutMultimesh(file, "PPOINTS", nblocks, pointMeshNames.data(), pointMeshTypes.data(), masteroptlist);

    putCellVars<RealType>(file, edgeFields, procPaths, nblocks, varTypes, masteroptlist);
    putCellVars<RealType>(file, faceFields, procPaths, nblocks, varTypes, masteroptlist);
    putCellVars<RealType>(file, cellFields, procPaths, nblocks, varTypes, masteroptlist);
#ifdef POLYTOPE_ENABLE_DEBUG
    std::vector<char*> nodeMeshNames;
    for (const auto& p : procPaths) {
      nodeMeshNames.push_back(strDup((p + "nodes").c_str()));
    }
    DBPutMultimesh(file, "NNODES", nblocks, nodeMeshNames.data(), pointMeshTypes.data(), masteroptlist);
    for (auto f = 0; f < nodeMeshNames.size(); ++f) {
      free(nodeMeshNames[f]);
    }
#endif

    DBClose(file);

    // Clean up
    DBFreeOptlist(masteroptlist);
    for (int i = 0; i < nblocks; ++i) {
      free(cellMeshNames[i]);
      free(pointMeshNames[i]);
    }
  }
  Communicator::Barrier();
#endif
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloWriter<2, double>;

} // end namespace
