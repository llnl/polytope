
#include "SiloReader.hh"
#include "polytope.hh"
#include "Communicator.hh"
#include "SiloUtils.hh"
#include "Serializer.hh"
#include "QuantTessellation.hh"

#include <algorithm>
#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include <type_traits>

#include "silo.h"

namespace polytope {

using namespace std;

namespace {

//-------------------------------------------------------------------
bool
tocHasName(const int n,
           char** names,
           const std::string& name) {
  for (int i = 0; i < n; ++i) {
    if (name == names[i]) return true;
  }
  return false;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
std::vector<std::string>
getSiloInputPaths(const std::string& masterFilename) {
  DBfile* file = DBOpen(masterFilename.c_str(), DB_UNKNOWN, DB_READ);
  POLY_ASSERT2(file, "Could not open file " << masterFilename);

  DBtoc* toc = DBGetToc(file);
  POLY_ASSERT2(toc, "Could not read Silo table of contents from file " << masterFilename);

  std::vector<std::string> result;
  const std::string globalMeshName = getGlobalMeshName();
  const std::string localMeshName = getLocalMeshName();

  if (tocHasName(toc->nmultimesh, toc->multimesh_names, globalMeshName)) {
    DBmultimesh* mmesh = DBGetMultimesh(file, globalMeshName.c_str());
    POLY_ASSERT2(mmesh, "Could not read multimesh " << globalMeshName << " from file " << masterFilename);
    result.resize(mmesh->nblocks);
    for (int i = 0; i < mmesh->nblocks; ++i) {
      result[i] = mmesh->meshnames[i];
    }
    DBFreeMultimesh(mmesh);
  } else if (tocHasName(toc->nucdmesh, toc->ucdmesh_names, globalMeshName)) {
    result.push_back(masterFilename + ":" + globalMeshName);
  } else if (tocHasName(toc->nucdmesh, toc->ucdmesh_names, localMeshName)) {
    result.push_back(masterFilename + ":" + localMeshName);
  } else {
    POLY_ASSERT2(false, "Could not find multimesh or mesh in file " << masterFilename);
  }

  DBClose(file);
  return result;
}
//-------------------------------------------------------------------

struct RequestedField {
  bool hasCentering;
  FieldCentering centering;
  std::string name;
};
//-------------------------------------------------------------------

}

//-------------------------------------------------------------------
template <typename TessType>
void
SiloReader<2, TessType>::read(TessType& mesh,
                              FieldTypeMap& fields,
                              const string& masterFilename) {
  using CoordType = typename std::decay<decltype(mesh.nodes[0].x)>::type;
  std::vector<std::string> inputFiles;
  int rank = 0;
  // Open a file in Silo/HDF5 format for reading.
#ifdef POLYTOPE_ENABLE_MPI
  auto& comm = Communicator::communicator();
  int nproc = Communicator::getNProcs();
  rank = Communicator::getRank();
  int root = Communicator::getRoot();
  int NBlocks = 0;
  std::vector<std::string> block_paths;
  if (rank == root) {
    block_paths = getSiloInputPaths(masterFilename);
    NBlocks = static_cast<int>(block_paths.size());
  }

  std::vector<char> buffer;
  if (rank == root) {
    serialize(block_paths, buffer);
  }
  int bsize = static_cast<int>(buffer.size());
  MPI_Bcast(&bsize, 1, MPI_INT, root, comm);
  if (rank != root) {
    buffer.resize(bsize);
  }
  if (bsize > 0) {
    MPI_Bcast(buffer.data(), bsize, MPI_CHAR, root, comm);
  }
  if (rank != root) {
    std::vector<char>::const_iterator bufIter = buffer.begin();
    deserialize<std::vector<std::string>>(block_paths, bufIter, buffer.end());
  }
  NBlocks = static_cast<int>(block_paths.size());

  for (int i = 0; i < NBlocks; ++i) {
    auto& cpath = block_paths[i];
    if (i % nproc == rank) {
      inputFiles.push_back(cpath);
    }
  }
#else
  inputFiles = getSiloInputPaths(masterFilename);
#endif

  for (const auto& cpath : inputFiles) {
    size_t colon_pos = cpath.find(":");
    std::string filename = cpath.substr(0, colon_pos);
    std::string internal_obj_path = getGlobalMeshName();
    if (colon_pos != std::string::npos) {
      internal_obj_path = cpath.substr(colon_pos + 1);
    }
    DBfile* file = DBOpen(filename.c_str(), DB_HDF5, DB_READ);
    if (!file) {
      printf("Rank %d input %s\n", rank, cpath.c_str());
      printf("Rank %d Failed to open file %s\n", rank, filename.c_str());
    }
    POLY_ASSERT2(file, "Could not open file " << filename);
    DBucdmesh* dbmesh = DBGetUcdmesh(file, internal_obj_path.c_str());
    POLY_ASSERT2(dbmesh, "Could not find mesh " << internal_obj_path << " in file " << filename);

    // Node coordinates.
    mesh.nodes.resize(dbmesh->nnodes);
    for (int i = 0; i < dbmesh->nnodes; ++i) {
      mesh.nodes[i].x = static_cast<CoordType>(((double*)(dbmesh->coords[0]))[i]);
      mesh.nodes[i].y = static_cast<CoordType>(((double*)(dbmesh->coords[1]))[i]);
    }

    // Reconstruct the cell-face and face-node connectivity.  The 2D writer
    // stores polygon zones and a separate compound array of polytope face ids.
    DBzonelist* zones = dbmesh->zones;
    POLY_ASSERT2(zones, "Could not find zonelist in file " << filename);
    const int ncells = zones->nzones;

    DBcompoundarray* conn = DBGetCompoundarray(file, "conn");
    POLY_ASSERT2(conn, "Could not find cell-face connectivity in file " << filename);
    // First element is the number of faces in each zone.
    // Second element is the list of face indices in each zone.
    // Third element is a pair of cells for each face.
    POLY_ASSERT2((conn->nelems == 2 && conn->elemlengths[0] == ncells),
                 "Found invalid cell-face connectivity in file " << filename);
    int* connData = (int*)conn->values;

    mesh.cells.resize(ncells);
    int foffset = ncells;
    int maxFace = -1;
    for (int c = 0; c < ncells; ++c) {
      int nfaces = connData[c];
      mesh.cells[c].resize(nfaces);
      copy(connData + foffset, connData + foffset + nfaces, mesh.cells[c].begin());
      for (int i = 0; i < nfaces; ++i) {
        const int iface = mesh.cells[c][i] < 0 ? ~mesh.cells[c][i] : mesh.cells[c][i];
        maxFace = std::max(maxFace, iface);
      }
      foffset += nfaces;
    }

    POLY_ASSERT2((conn->nvalues - foffset) % 2 == 0,
                 "Found invalid face-cell connectivity in file " << filename);
    const int nfaces = maxFace+1;
    mesh.faces.resize(nfaces);

    int zoffset = 0;
    int shape = 0;
    int remainingInShape = zones->nshapes > 0 ? zones->shapecnt[0] : 0;
    for (int c = 0; c < ncells; ++c) {
      while (remainingInShape == 0) {
        ++shape;
        POLY_ASSERT(shape < zones->nshapes);
        remainingInShape = zones->shapecnt[shape];
      }

      int ncellNodes = zones->shapesize[shape];
      if (ncellNodes == 0) {
        POLY_ASSERT(zoffset < zones->lnodelist);
        ncellNodes = zones->nodelist[zoffset++];
      }

      std::vector<unsigned> cellNodes(ncellNodes);
      for (int i = 0; i < ncellNodes; ++i) {
        POLY_ASSERT(zoffset < zones->lnodelist);
        cellNodes[i] = zones->nodelist[zoffset++];
      }

      const bool closed = (cellNodes.size() > 1 && cellNodes.front() == cellNodes.back());
      const int nedges = closed ? ncellNodes - 1 : ncellNodes;
      POLY_ASSERT(nedges == static_cast<int>(mesh.cells[c].size()));

      for (int i = 0; i < nedges; ++i) {
        const int signedFace = mesh.cells[c][i];
        const int iface = signedFace < 0 ? ~signedFace : signedFace;
        POLY_ASSERT(iface < static_cast<int>(mesh.faces.size()));
        const unsigned n0 = cellNodes[i];
        const unsigned n1 = cellNodes[(i + 1) % ncellNodes];
        if (mesh.faces[iface].empty()) {
          mesh.faces[iface].resize(2);
          if (signedFace >= 0) {
            mesh.faces[iface][0] = n0;
            mesh.faces[iface][1] = n1;
          } else {
            mesh.faces[iface][0] = n1;
            mesh.faces[iface][1] = n0;
          }
        }
      }

      --remainingInShape;
    }
    mesh.computeFaceCells();
    DBFreeUcdmesh(dbmesh);
    DBFreeCompoundarray(conn);
    // Get the generator points
    DBpointmesh* pmesh = DBGetPointmesh(file, "points");
    POLY_ASSERT2(pmesh, "Could not find generator points in file " << filename);
    int npts = pmesh->nels;
    mesh.points.resize(npts);
    for (int i = 0; i < npts; ++i) {
      mesh.points[i].x = static_cast<CoordType>(((double*)(pmesh->coords[0]))[i]);
      mesh.points[i].y = static_cast<CoordType>(((double*)(pmesh->coords[1]))[i]);
    }
    DBFreePointmesh(pmesh);

    // Get the cell field variables
    std::vector<RequestedField> fieldNames;
    if (fields.empty()) {
      DBtoc* contents = DBGetToc(file);
      for (int f = 0; f < contents->nucdvar; ++f) {
        fieldNames.push_back({false, FieldCentering::Cell, std::string(contents->ucdvar_names[f])});
      }
    } else {
      for (const auto& [centering, fieldmap] : fields) {
        for (const auto& field : fieldmap) {
          fieldNames.push_back({true, centering, field.first});
        }
      }
    }
    // Retrieve the fields.
    for (size_t f = 0; f < fieldNames.size(); ++f) {
      const auto& requestedField = fieldNames[f];
      const auto& fieldName = requestedField.name;
      DBucdvar* dbvar = DBGetUcdvar(file, fieldName.c_str());
      POLY_ASSERT2(dbvar, "Could not find field " << fieldName << " in file " << filename);
      FieldCentering centering = requestedField.centering;
      if (requestedField.hasCentering) {
        const int expectedCentering = siloCentering(requestedField.centering);
        POLY_ASSERT2(dbvar->centering == expectedCentering,
                     "Field " << fieldName << " has unexpected centering in file " << filename);
      } else {
        centering = fieldCenteringFromSilo(dbvar->centering);
      }
      fields[centering][fieldName].resize(dbvar->nels);
      copy((double*)(dbvar->vals[0]), (double*)(dbvar->vals[0]) + dbvar->nels,
           &(fields[centering][fieldName][0]));
      // Clean up.
      DBFreeUcdvar(dbvar);
    }
    DBClose(file);
  }
}

//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloReader<2, Tessellation<2, double>>;
template class SiloReader<2, QuantTessellation<2>>;

namespace Silo {

//-------------------------------------------------------------------
// Helper function definition
//-------------------------------------------------------------------
vector<int> findAvailableCycles(const string& prefix,
                                const string& directory) {
#ifdef POLYTOPE_ENABLE_MPI
  const int nproc = Communicator::getNProcs();
#endif

  // If the directory is not given, infer it from the prefix.
  string dir = directory;
  if (dir.empty()) {
#ifdef POLYTOPE_ENABLE_MPI
    dir = prefix + "-" + std::to_string(nproc);
#else
    dir = ".";
#endif
  }

  // Now find files that match the prefix in the directory.
  vector<int> cycles;
  DIR* d = opendir(dir.c_str());
  POLY_ASSERT2(d, "Could not find the directory " << dir);
  // if (d == 0) {
  //   error("Could not find the directory " + dir);
  // }
  dirent* entry;
  while ((entry = readdir(d)) != 0) {
    std::string path = dir + "/" + entry->d_name;
    if (entry->d_type != DT_DIR) {
      if ((path.find(prefix) != std::string::npos) &&
          (path.find(".silo") != std::string::npos)) {
        // Pull the cycle number out of the path.
        string p = string(p);
        int suffix = p.find(".silo");
        int dash = p.rfind(suffix);
        int cycle = atoi(p.substr(dash+1, suffix-dash-1).c_str());
        cycles.push_back(cycle);
      }
    }
  }

  return cycles;
}
//-------------------------------------------------------------------

}

} // end namespace
