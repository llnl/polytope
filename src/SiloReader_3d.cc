
#include "SiloReader.hh"
#include "polytope.hh"
#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"

#ifdef POLYTOPE_ENABLE_MPI
#include "mpi.h"
#include "pmpio.h"
#else
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif

namespace polytope {

using namespace std;

#ifdef POLYTOPE_ENABLE_MPI
namespace {
//-------------------------------------------------------------------
void*
PMPIO_createFile(const char* filename,
                 const char* dirname,
                 void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void* 
PMPIO_openFile(const char* filename, 
               const char* dirname,
               PMPIO_iomode_t iomode, 
               void* userData)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (iomode == PMPIO_WRITE)
  { 
    file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
    DBMkDir(file, dirname);
    DBSetDir(file, dirname);
  }
  else
  {
    file = DBOpen(filename, driver, DB_READ);
    DBSetDir(file, dirname);
  }
  return (void*)file;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void
PMPIO_closeFile(void* file,
                void* userData)
{
  DBClose((DBfile*)file);
}
//-------------------------------------------------------------------
} // end namespace
#endif

//-------------------------------------------------------------------
template <typename RealType>
void 
SiloReader<3, RealType>::
read(Tessellation<3, RealType>& mesh, 
     map<string, vector<RealType> >& fields,
     map<string, vector<int> >& nodeTags,
     map<string, vector<int> >& edgeTags,
     map<string, vector<int> >& faceTags,
     map<string, vector<int> >& cellTags,
     const string& filePrefix,
     const string& directory,
     int cycle,
     RealType& time,
     MPI_Comm comm,
     int numFiles,
     int mpiTag)
{
  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  int index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);

  // Open a file in Silo/HDF5 format for reading.
#ifdef POLYTOPE_ENABLE_MPI
  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (numFiles == -1)
    numFiles = nproc;
  POLY_ASSERT(numFiles <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Check for the existence of the master directory.
  string masterDirName = directory;
  if (masterDirName.empty()) {
    masterDirName = filePrefix + "-" + std::to_string(nproc);
  }
  DIR* masterDir = opendir(directory.c_str());
  if (masterDir == 0) {
    error("Could not find the directory " + masterDirName);
  }

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_baton_t* baton = PMPIO_Init(numFiles, PMPIO_READ, comm, mpiTag, 
                                    &PMPIO_createFile, 
                                    &PMPIO_openFile, 
                                    &PMPIO_closeFile,
                                    0);
  int groupRank = PMPIO_GroupRank(baton, rank);
  int rankInGroup = PMPIO_RankInGroup(baton, rank);

  // Figure out the subdirectory for this group.
  std::string groupdirname = masterDirName + "/" + std::to_string(groupRank);
  DIR* groupDir = opendir(groupdirname.c_str());
  if (groupDir == 0) {
    error("Could not find the directory " + groupdirname);
  }
  std::string filename;
  // Determine the file name.
  if (cycle >= 0) {
    filename = groupdirname + "/" + prefix + "-" + std::to_string(cycle) + ".silo";
  } else {
    filename = groupdirname + "/" + prefix + ".silo";
  }

  std::string dirname = "domain_" + std::to_string(rankInGroup);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename.c_str(), dirname.c_str());
  DBSetDir(file, dirname.c_str());
#else
  string dirname = directory;
  if (dirname.empty()) dirname = ".";
  std::string filename;
  // Determine the file name.
  if (cycle >= 0) {
    filename = dirname + "/" + prefix + "-" + std::to_string(cycle) + ".silo";
  } else {
    filename = dirname + "/" + prefix + ".silo";
  }

  int driver = DB_HDF5;
  DBfile* file = DBOpen(filename.c_str(), driver, DB_READ);
  DBSetDir(file, "/");
#endif

  // Retrieve the mesh. Note that we must deallocate the storage 
  // for this object after we're through!
  DBucdmesh* dbmesh = DBGetUcdmesh(file, "mesh");

  // Extract time.
  time = dbmesh->dtime;

  // Node coordinates.
  mesh.nodes.resize(3*dbmesh->nnodes);
  for (int i = 0; i < dbmesh->nnodes; ++i)
  {
    mesh.nodes[3*i]  = ((RealType*)(dbmesh->coords[i]))[0];
    mesh.nodes[3*i+1] = ((RealType*)(dbmesh->coords[i]))[1];
    mesh.nodes[3*i+2] = ((RealType*)(dbmesh->coords[i]))[2];
  }

  // Reconstruct the faces.
  mesh.faces.resize(dbmesh->faces->nfaces);
  int noffset = 0;
  for (int f = 0; f < mesh.faces.size(); ++f)
  {
    mesh.faces[f].resize(dbmesh->faces->shapesize[f]);
    for (int n = 0; n < mesh.faces[f].size(); ++n, ++noffset)
      mesh.faces[f][n] = dbmesh->faces->nodelist[noffset];
  }

  // Reconstruct the cell-face connectivity.
  mesh.cells.resize(dbmesh->zones->nzones);
  DBcompoundarray* conn = DBGetCompoundarray(file, "connectivity");
  if (conn == 0)
  {
    DBClose(file);
    error("Could not find cell-face connectivity in file " + filename);
  }
  // First element is the number of faces in each zone.
  // Second element is the list of face indices in each zone.
  // Third element is a pair of cells for each face.
  if ((conn->nelems != 3) or 
      (conn->elemlengths[0] != dbmesh->zones->nzones) or 
      (conn->elemlengths[2] != 2*dbmesh->faces->nfaces))
  {
    DBClose(file);
    error("Found invalid cell-face connectivity in file " + filename);
  }
  int* connData = (int*)conn->values;
  int foffset = dbmesh->zones->nzones;
  for (int c = 0; c < dbmesh->zones->nzones; ++c)
  {
    int nfaces = connData[c];
    mesh.cells[c].resize(nfaces);
    copy(connData + foffset, connData + foffset + nfaces, mesh.cells[c].begin());
    foffset += nfaces;
  }
  mesh.faceCells.resize(mesh.faces.size());
  for (size_t f = 0; f < mesh.faceCells.size(); ++f)
  {
    mesh.faceCells[f].resize(2);
    mesh.faceCells[f][0] = connData[foffset];
    mesh.faceCells[f][1] = connData[foffset+1];
    foffset += 2;
  }
  DBFreeUcdmesh(dbmesh);
  DBFreeCompoundarray(conn);
  
  // Check for convex hull data.
  // First element is the number of facets.
  // Second element is the array of numbers of nodes per facet.
  // Third element is the array of node indices for the facets.
  DBcompoundarray* hull = DBGetCompoundarray(file, "convexhull");
  if (hull != 0)
  {
    if ((hull->nelems != 3) or (hull->elemlengths[0] != 1))
    {
      DBClose(file);
      error("Found invalid convex hull data in file " + filename);
    }
    int* hullData = (int*)conn->values;
    int nfacets = hullData[0];
    mesh.convexHull.facets.resize(nfacets);
    int foffset = 1;
    for (int f = 0; f < nfacets; ++f, ++foffset)
    {
      int nnodes = hullData[foffset];
      mesh.convexHull.facets[f].resize(nnodes);
    }
    for (int f = 0; f < nfacets; ++f, ++foffset)
    {
      for (int n = 0; n < mesh.convexHull.facets[f].size(); ++n)
        mesh.convexHull.facets[f][n] = hullData[foffset];
    }
    DBFreeCompoundarray(hull);
  }

  // FIXME: Check for hole data?

  // Read any tag data.
  DBcompoundarray* tags = DBGetCompoundarray(file, "node_tags");
  if (tags != 0)
  {
    for (int i = 0; i < tags->nelems; ++i)
    {
      std::vector<int>& tag = nodeTags[tags->elemnames[i]];
      tag.resize(tags->elemlengths[i]);
      copy((int*)tags->values, (int*)((char*)tags->values + tags->elemlengths[i]), tag.begin());
    }
    DBFreeCompoundarray(tags);
  }
  tags = DBGetCompoundarray(file, "edge_tags");
  if (tags != 0)
  {
    for (int i = 0; i < tags->nelems; ++i)
    {
      std::vector<int>& tag = edgeTags[tags->elemnames[i]];
      tag.resize(tags->elemlengths[i]);
      copy((int*)tags->values, (int*)((char*)tags->values + tags->elemlengths[i]), tag.begin());
    }
    DBFreeCompoundarray(tags);
  }
  tags = DBGetCompoundarray(file, "face_tags");
  if (tags != 0)
  {
    for (int i = 0; i < tags->nelems; ++i)
    {
      std::vector<int>& tag = faceTags[tags->elemnames[i]];
      tag.resize(tags->elemlengths[i]);
      copy((int*)tags->values, (int*)((char*)tags->values + tags->elemlengths[i]), tag.begin());
    }
    DBFreeCompoundarray(tags);
  }
  tags = DBGetCompoundarray(file, "cell_tags");
  if (tags != 0)
  {
    for (int i = 0; i < tags->nelems; ++i)
    {
      std::vector<int>& tag = cellTags[tags->elemnames[i]];
      tag.resize(tags->elemlengths[i]);
      copy((int*)tags->values, (int*)((char*)tags->values + tags->elemlengths[i]), tag.begin());
    }
    DBFreeCompoundarray(tags);
  }

  // Make a list of the desired fields.
  vector<string> fieldNames;
  if (fields.empty()) {
    DBtoc* contents = DBGetToc(file);
    for (int f = 0; f < contents->nucdvar; ++f)
      fieldNames.push_back(string(contents->ucdvar_names[f]));
  } else {
    for (typename map<string, vector<RealType> >::const_iterator iter = fields.begin();
         iter != fields.end(); ++iter)
    fieldNames.push_back(iter->first);
  }

  // Retrieve the fields.
  for (int f = 0; f < fieldNames.size(); ++f) {
    DBucdvar* dbvar = DBGetUcdvar(file, fieldNames[f].c_str());
    if (dbvar == 0) {
      DBClose(file);
      error("Could not find field " + fieldNames[f] + " in file " + filename);
    }
    fields[fieldNames[f]].resize(dbvar->nels);
    copy((RealType*)(dbvar->vals[0]), (RealType*)(dbvar->vals[0]) + dbvar->nels,
         &(fields[fieldNames[f]][0]));

    // Clean up.
    DBFreeUcdvar(dbvar);
  }

  // Clean up.
  DBClose(file);
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloReader<3, double>;

} // end namespace
