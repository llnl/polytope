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
#include "silo.h"

namespace polytope {

using namespace std;

//-------------------------------------------------------------------
template <typename TessType>
void
SiloWriter<3, TessType>::write(const string& filePrefix,
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
    char *coordnames[3];
    coordnames[0] = (char*)"xcoords";
    coordnames[1] = (char*)"ycoords";
    coordnames[2] = (char*)"zcoords";

    // Node coordinates.
    int numNodes = m_mesh.nodes.size();
    vector<double> x(numNodes), y(numNodes), z(numNodes);
    for (int i = 0; i < numNodes; ++i) {
      x[i] = static_cast<double>(m_mesh.nodes[i].x);
      y[i] = static_cast<double>(m_mesh.nodes[i].y);
      z[i] = static_cast<double>(m_mesh.nodes[i].z);
    }
    double* coords[3];
    coords[0] = &(x[0]);
    coords[1] = &(y[0]);
    coords[2] = &(z[0]);

    const auto numFaces = m_mesh.faces.size();
    vector<int> faceNodeCounts, allFaceNodes;
    faceNodeCounts.reserve(numFaces);
    for (auto iface = 0u; iface < numFaces; ++iface) {
      faceNodeCounts.push_back(m_mesh.faces[iface].size());
      for (auto& face : m_mesh.faces[iface]) {
        allFaceNodes.push_back(static_cast<int>(face));
      }
      //std::copy(mesh.faces[iface].begin(), mesh.faces[iface].end(), std::back_inserter(allFaceNodes));
    }
    POLY_ASSERT(faceNodeCounts.size() == numFaces);

    // All zones are polygonal.
    const auto numCells = m_mesh.cells.size();
    vector<int> cellFaceCounts, allCellFaces;
    cellFaceCounts.reserve(numCells);
    for (auto i = 0u; i < numCells; ++i) {
      auto n = m_mesh.cells[i].size();
      cellFaceCounts.push_back(n);
      std::copy(m_mesh.cells[i].begin(), m_mesh.cells[i].end(), std::back_inserter(allCellFaces));
    }
    vector<char> boundaryFaceFlags(numFaces, 0x0);
    // for (vector<unsigned>::const_iterator itr = mesh.boundaryFaces.begin();
    //      itr != mesh.boundaryFaces.end();
    //      ++itr) {
    //   POLY_ASSERT(*itr < numFaces);
    //   boundaryFaceFlags[*itr] = 0x1;
    // }

    // The polyhedral zone list is referred to in the options list.
    DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"zonelist");

    // Write out the 3D polyhedral mesh.
    DBPutUcdmesh(file, meshname.c_str(), 3, coordnames, coords,
                 numNodes, numCells,
                 NULL, NULL, DB_DOUBLE, optlist);
    // Write the connectivity information.
    DBPutPHZonelist(file, (char*)"zonelist",
                    faceNodeCounts.size(), &faceNodeCounts[0],
                    allFaceNodes.size(), &allFaceNodes[0],
                    &boundaryFaceFlags[0],
                    cellFaceCounts.size(), &cellFaceCounts[0],
                    allCellFaces.size(), &allCellFaces[0],
                    0, 0, numCells-1, optlist);

    // Write out the cell-face connectivity data.
    vector<int> conn(numCells);
    int elemlengths[2];
    char* elemnames[2];
    elemnames[0] = strDup("ncellfaces");
    elemlengths[0] = numCells;
    for (auto c = 0u; c < numCells; ++c) {
      conn[c] = m_mesh.cells[c].size();
    }
    for (auto c = 0u; c < numCells; ++c) {
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

    writeFieldsToFile(meshname, file, optlist);

    int numPoints = m_mesh.points.size();
    vector<double> xp(numPoints), yp(numPoints), zp(numPoints);
    for (int i = 0; i < numPoints; ++i) {
      xp[i] = static_cast<double>(m_mesh.points[i].x);
      yp[i] = static_cast<double>(m_mesh.points[i].y);
      zp[i] = static_cast<double>(m_mesh.points[i].z);
    }
    double* pcoords[3];
    pcoords[0] = &(xp[0]);
    pcoords[1] = &(yp[0]);
    pcoords[2] = &(zp[0]);
    // Write point mesh
    DBPutPointmesh(file, (char*)"points", 3, pcoords, numPoints, DB_DOUBLE, optlist);
#ifdef POLYTOPE_ENABLE_DEBUG
    // Create NODES directory and write the nodes as points
    // Node coordinates
    // Write point mesh of nodes
    DBPutPointmesh(file, (char*)"nodes", 3, coords, numNodes, DB_DOUBLE, optlist);
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
template class SiloWriter<3, Tessellation<3, double>>;
template class SiloWriter<3, QuantTessellation<3>>;

} // end namespace
