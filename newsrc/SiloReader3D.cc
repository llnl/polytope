
#include "SiloReader.hh"
#include "polytope.hh"
#include "Communicator.hh"
#include "SiloUtils.hh"
#include "Serializer.hh"
#include "Tessellation.hh"

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
SiloReader<3, TessType>::
read(TessType& mesh,
     FieldTypeMap& fields,
     const string& masterFilename)
{
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

    // Retrieve the mesh. Note that we must deallocate the storage
    // for this object after we're through!
    DBucdmesh* dbmesh = DBGetUcdmesh(file, internal_obj_path.c_str());
    POLY_ASSERT2(dbmesh, "Could not find mesh " << internal_obj_path << " in file " << filename);

    // Node coordinates.
    mesh.nodes.resize(dbmesh->nnodes);
    for (int i = 0; i < dbmesh->nnodes; ++i)
    {
      mesh.nodes[i].x = static_cast<CoordType>(((double*)(dbmesh->coords[0]))[i]);
      mesh.nodes[i].y = static_cast<CoordType>(((double*)(dbmesh->coords[1]))[i]);
      mesh.nodes[i].z = static_cast<CoordType>(((double*)(dbmesh->coords[2]))[i]);
    }

    // Reconstruct the faces.
    mesh.faces.resize(dbmesh->faces->nfaces);
    int noffset = 0;
    for (size_t f = 0; f < mesh.faces.size(); ++f)
    {
      mesh.faces[f].resize(dbmesh->faces->shapesize[f]);
      for (size_t n = 0; n < mesh.faces[f].size(); ++n, ++noffset)
        mesh.faces[f][n] = dbmesh->faces->nodelist[noffset];
    }

    // Reconstruct the cell-face connectivity.
    mesh.cells.resize(dbmesh->zones->nzones);
    DBcompoundarray* conn = DBGetCompoundarray(file, "conn");
    POLY_ASSERT2(conn, "Could not find cell-face connectivity in file " << filename);
    // First element is the number of faces in each zone.
    // Second element is the list of face indices in each zone.
    POLY_ASSERT2((conn->nelems == 2 &&
                  conn->elemlengths[0] == dbmesh->zones->nzones),
                 "Found invalid cell-face connectivity in file " << filename);
    int* connData = (int*)conn->values;
    int foffset = dbmesh->zones->nzones;
    for (int c = 0; c < dbmesh->zones->nzones; ++c)
    {
      int nfaces = connData[c];
      mesh.cells[c].resize(nfaces);
      copy(connData + foffset, connData + foffset + nfaces, mesh.cells[c].begin());
      foffset += nfaces;
    }
    mesh.computeFaceCells();
    DBFreeUcdmesh(dbmesh);
    DBFreeCompoundarray(conn);

    // Check for convex hull data.
    // First element is the number of facets.
    // Second element is the array of numbers of nodes per facet.
    // Third element is the array of node indices for the facets.
    DBcompoundarray* hull = DBGetCompoundarray(file, "convexhull");
    if (hull != 0)
    {
      POLY_ASSERT2((hull->nelems == 3 && hull->elemlengths[0] == 1),
                   "Found invalid convex hull data in file " << filename);
      int* hullData = (int*)hull->values;
      int nfacets = hullData[0];
      mesh.convexHull.facets.resize(nfacets);
      int foffset = 1;
      for (int f = 0; f < nfacets; ++f, ++foffset)
      {
        int nnodes = hullData[foffset];
        mesh.convexHull.facets[f].resize(nnodes);
      }
      for (int f = 0; f < nfacets; ++f)
      {
        for (size_t n = 0; n < mesh.convexHull.facets[f].size(); ++n, ++foffset)
          mesh.convexHull.facets[f][n] = hullData[foffset];
      }
      DBFreeCompoundarray(hull);
    }

    DBpointmesh* pmesh = DBGetPointmesh(file, "points");
    if (pmesh != 0) {
      const int npts = pmesh->nels;
      mesh.points.resize(npts);
      for (int i = 0; i < npts; ++i) {
        mesh.points[i].x = static_cast<CoordType>(((double*)(pmesh->coords[0]))[i]);
        mesh.points[i].y = static_cast<CoordType>(((double*)(pmesh->coords[1]))[i]);
        mesh.points[i].z = static_cast<CoordType>(((double*)(pmesh->coords[2]))[i]);
      }
      DBFreePointmesh(pmesh);
    }

    // Make a list of the desired fields.
    // TODO: Retrieve all other field types as well
    vector<RequestedField> fieldNames;
    if (fields.empty()) {
      DBtoc* contents = DBGetToc(file);
      for (int f = 0; f < contents->nucdvar; ++f)
        fieldNames.push_back({false, FieldCentering::Cell, string(contents->ucdvar_names[f])});
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
      fields[centering][fieldName].resize(dbvar->nels);
      copy((double*)(dbvar->vals[0]), (double*)(dbvar->vals[0]) + dbvar->nels,
           &(fields[centering][fieldName][0]));

      // Clean up.
      DBFreeUcdvar(dbvar);
    }

    // Clean up.
    DBClose(file);
  }
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloReader<3, Tessellation<3, double>>;

} // end namespace
