#include "polytope.hh"
#include "SiloWriter.hh"
#include "Tessellation.hh"
#include "QuantTessellation.hh"
#include "Communicator.hh"

#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>

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
template <typename TessType>
void
traverseNodes(const TessType& mesh,
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
  for (auto n = 0u; n < nodes.size(); ++n) {
    POLY_ASSERT(nodes[n] >= 0);
    POLY_ASSERT(nodes[n] < int(mesh.nodes.size()));
  }
#endif
}
//-------------------------------------------------------------------
}

//-------------------------------------------------------------------
template <typename TessType>
void
SiloWriter<2, TessType>::write(const string& filePrefix,
                               const string& directory,
                               int cycle,
                               double time,
                               int numFiles) {
  int nranks = Communicator::getNRanks();
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  auto& comm = Communicator::communicator();
  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  int index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);
  int coord_sys = DB_CARTESIAN;
  string dirname = directory;
  if (dirname.empty()) dirname = ".";
  std::string filename = getMasterFilename(prefix, cycle);

  std::string meshname = getGlobalMeshName();
  bool hasPoints = true;
  // Open a file in Silo/HDF5 format for writing.
#ifdef POLYTOPE_ENABLE_MPI
  bool doParallel = false;
  std::string masterDirname = "";
  std::vector<int> ranksWithData;
  if (nranks == 1) {
    numFiles = 1;
  }
  if (numFiles == -1 || numFiles > 1) {
    doParallel = true;
    int localRankHasPoints = (m_mesh.points.size() > 0) ? 1 : 0;
    ranksWithData = std::move(gatherValidRanks(localRankHasPoints));
    if (numFiles == -1) {
      int globalWriteProcs = static_cast<int>(ranksWithData.size());
      MPI_Bcast(&globalWriteProcs, 1, MPI_INT, root, comm);
      numFiles = globalWriteProcs;
    }
    POLY_ASSERT(numFiles <= nranks);

    masterDirname = getMasterDirname(directory, prefix, cycle);
    if (rank == root) {
      DIR* masterDir = opendir(masterDirname.c_str());
      if (masterDir == 0) {
        mkdir(masterDirname.c_str(), S_IRWXU | S_IRWXG);
      } else {
        closedir(masterDir);
      }
    }
    Communicator::Barrier();
    hasPoints = bool(localRankHasPoints);

    filename = getFilename(masterDirname, rank);
    meshname = getLocalMeshName();
  }
#endif
  if (hasPoints) {
    DBfile* file = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, 0, DB_HDF5);
    // Add cycle/time metadata if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (time >= 0.)
      DBAddOption(optlist, DBOPT_DTIME, &time);

    DBAddOption(optlist, DBOPT_COORDSYS, &coord_sys);
    // This is optional for now, but we'll give it anyway.
    char *coordnames[2];
    coordnames[0] = (char*)"xcoords";
    coordnames[1] = (char*)"ycoords";

    // Node coordinates.
    int numNodes = m_mesh.nodes.size();
    vector<double> x(numNodes), y(numNodes);
    for (int i = 0; i < numNodes; ++i) {
      x[i] = static_cast<double>(m_mesh.nodes[i].x);
      y[i] = static_cast<double>(m_mesh.nodes[i].y);
    }
    double* coords[2];
    coords[0] = &(x[0]);
    coords[1] = &(y[0]);

    int numCells = m_mesh.cells.size();
    // All zones are polygonal.
    const int numShapes = numCells;
    vector<int> shapesize(numShapes, 0),
      shapetype(numShapes, DB_ZONETYPE_POLYGON),
      shapecount(numShapes, 1),
      nodeList;
    for (int i = 0; i < numCells; ++i) {
      // Gather the nodes from this cell in traversal order.
      vector<int> cellNodes;
      traverseNodes<TessType>(m_mesh, i, cellNodes);
      // Insert the cell's node connectivity into the node list.
      nodeList.push_back(cellNodes.size());
      nodeList.insert(nodeList.end(), cellNodes.begin(), cellNodes.end());
    }

    // Write out the 2D polygonal mesh.
    DBPutUcdmesh(file, meshname.c_str(), 2, coordnames, coords,
                 numNodes, numCells,
                 "zonelist", NULL, DB_DOUBLE, optlist);
    DBPutZonelist2(file, "zonelist", numCells,
                   2, &nodeList[0], nodeList.size(), 0, 0, 0,
                   &shapetype[0], &shapesize[0], &shapecount[0],
                   numShapes, optlist);

    // Write out the cell-face connectivity data.
    // This is used when Polytope reads a Silo file it wrote.
    vector<int> conn(numCells);
    int elemlengths[2];
    char* elemnames[2];
    elemnames[0] = strDup("ncellfaces");
    elemlengths[0] = numCells;
    for (int c = 0; c < numCells; ++c) {
      conn[c] = m_mesh.cells[c].size();
    }
    for (int c = 0; c < numCells; ++c) {
      for (size_t f = 0; f < m_mesh.cells[c].size(); ++f) {
        conn.push_back(m_mesh.cells[c][f]);
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

    // FIXME: We really should try to use the number of edges for edge fields.
    writeFieldsToFile(meshname, file, optlist);

    int numPoints = m_mesh.points.size();
    vector<double> xp(numPoints), yp(numPoints);
    for (int i = 0; i < numPoints; ++i) {
      xp[i] = static_cast<double>(m_mesh.points[i].x);
      yp[i] = static_cast<double>(m_mesh.points[i].y);
    }
    double* pcoords[2];
    pcoords[0] = &(xp[0]);
    pcoords[1] = &(yp[0]);
    // Write point mesh
    DBPutPointmesh(file, (char*)"points", 2, pcoords, numPoints, DB_DOUBLE, optlist);
#ifdef POLYTOPE_ENABLE_DEBUG
    // Create NODES directory and write the nodes as points
    // Node coordinates
    DBPutPointmesh(file, (char*)"nodes", 2, coords, numNodes, DB_DOUBLE, optlist);
#endif

    // Clean up.
    DBClose(file);
    DBFreeOptlist(optlist);
  }

#ifdef POLYTOPE_ENABLE_MPI
  // Finally, write the uber-master file.
  if (rank == root && doParallel) {
    std::string masterFilename = getMasterFilename(prefix, cycle);
    writeMasterFile(masterFilename, masterDirname, ranksWithData,
                    time, cycle);
  }
  if (doParallel) {
    Communicator::Barrier();
  }
#endif
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloWriter<2, Tessellation<2, double>>;
template class SiloWriter<2, QuantTessellation<2>>;

} // end namespace
