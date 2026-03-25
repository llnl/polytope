
#include "polytope.hh"
#include "SiloWriter.hh"

#include <fstream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include "silo.h"

#ifdef POLYTOPE_ENABLE_MPI
#include "pmpio.h"
#endif

#include "SiloUtils.hh"

namespace polytope {

using namespace std;

namespace {

//-------------------------------------------------------------------
// Traverse the given points of a polygonal facet along their convex
// hull, writing their indices to indices.
//-------------------------------------------------------------------
template <typename RealType>
void 
traverseConvexHull(const vector<RealType>& points,
                   vector<int>& indices)
{
  // Find the "lowest" point in the set.
  RealType ymin = FLT_MAX;
  int index0 = -1;
  int numPoints = points.size() / 2;
  for (int p = 0; p < numPoints; ++p)
  {
    if (ymin > points[2*p+1])
    {
      ymin = points[2*p+1];
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  RealType thetaPrev = 0.0;
  indices.push_back(index0);

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    RealType dthetaMin = 2.0*M_PI;
    int jMin = -1;
    for (int j = 0; j < numPoints; ++j)
    {
      if (j != i)
      {
        RealType dx = points[2*j] - points[2*i],
             dy = points[2*j+1] - points[2*i+1];
        RealType theta = atan2(dy, dx);
        RealType dtheta = theta - thetaPrev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dthetaMin > dtheta)
        {
          dthetaMin = dtheta;
          jMin = j;
        }
      }
    }
    if (jMin != index0)
      indices.push_back(jMin);
    thetaPrev += dthetaMin;
    i = jMin;
  }
  while (i != index0);

  // The convex hull should be a polygon unless the input points 
  // don't form a polygon.
  POLY_ASSERT((numPoints <= 2) or 
         ((numPoints > 2) and (indices.size() > 2)));
}
//-------------------------------------------------------------------

#ifdef POLYTOPE_ENABLE_MPI

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
#endif

}

//-------------------------------------------------------------------
template <typename RealType>
void 
SiloWriter<3, RealType>::
write(const Tessellation<3, RealType>& mesh, 
      const std::map<std::string, RealType*>& nodeFields,
      const std::map<std::string, std::vector<int>*>& nodeTags,
      const std::map<std::string, RealType*>& edgeFields,
      const std::map<std::string, std::vector<int>*>& edgeTags,
      const std::map<std::string, RealType*>& faceFields,
      const std::map<std::string, std::vector<int>*>& faceTags,
      const std::map<std::string, RealType*>& cellFields,
      const std::map<std::string, std::vector<int>*>& cellTags,
      const string& filePrefix,
      const string& directory,
      int cycle,
      RealType time,
      MMPI_Comm comm,
      int numFiles,
      int mpiTag)
{
  // Strip .silo off of the prefix if it's there.
  string prefix = filePrefix;
  int index = prefix.find(".silo");
  if (index >= 0)
    prefix.erase(index);

  // Open a file in Silo/HDF5 format for writing.
  char filename[1024];
#ifdef POLYTOPE_ENABLE_MPI
  int nproc = 1, rank = 0;
  MMPI_Comm_size(comm, &nproc);
  MMPI_Comm_rank(comm, &rank);
  if (numFiles == -1)
    numFiles = nproc;
  POLY_ASSERT(numFiles <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Create the master directory if we need to.
  string masterDirName = directory;
  if (masterDirName.empty())
  {
    char dirname[1024];
    snprintf(dirname, 1024, "%s-%d", filePrefix.c_str(), nproc);
    masterDirName = dirname;
  }
  if (rank == 0)
  {
    DIR* masterDir = opendir(masterDirName.c_str());
    if (masterDir == 0)
      mkdir((char*)masterDirName.c_str(), S_IRWXU | S_IRWXG);
    else
      closedir(masterDir);
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_baton_t* baton = PMPIO_Init(numFiles, PMPIO_WRITE, comm, mpiTag, 
                                    &PMPIO_createFile, 
                                    &PMPIO_openFile, 
                                    &PMPIO_closeFile,
                                    0);
  int groupRank = PMPIO_GroupRank(baton, rank);
  int rankInGroup = PMPIO_RankInGroup(baton, rank);

  // Create a subdirectory for each group.
  char groupdirname[1024];
  snprintf(groupdirname, 1024, "%s/%d", masterDirName.c_str(), groupRank);
  if (rankInGroup == 0)
  {
    DIR* groupDir = opendir(groupdirname);
    if (groupDir == 0)
      mkdir((char*)groupdirname, S_IRWXU | S_IRWXG);
    else
      closedir(groupDir);
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  // Determine a file name.
  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", groupdirname, prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", groupdirname, prefix.c_str());

  char dirname[1024];
  snprintf(dirname, 1024, "domain_%d", rankInGroup);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dirname);
#else
  string dirname = directory;
  if (dirname.empty())
    dirname = ".";

  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", dirname.c_str(), prefix.c_str(), cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", dirname.c_str(), prefix.c_str());

  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(file, "/");
#endif

  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  double dtime = static_cast<double>(time);
  if (cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (dtime != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int numNodes = mesh.nodes.size() / 3;
  vector<double> x(numNodes), y(numNodes), z(numNodes);
  for (int i = 0; i < numNodes; ++i)
  {
    x[i] = mesh.nodes[3*i];
    y[i] = mesh.nodes[3*i+1];
    z[i] = mesh.nodes[3*i+2];
  }
  double* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Construct the silo face-node info.  We rely on the input tessellation having
  // the faces nodes arranged counter-clockwise around the face.
  const int numFaces = mesh.faces.size();
  vector<int> faceNodeCounts, allFaceNodes;
  faceNodeCounts.reserve(numFaces);
  for (int iface = 0; iface != numFaces; ++iface) 
  {
    faceNodeCounts.push_back(mesh.faces[iface].size());
    std::copy(mesh.faces[iface].begin(), mesh.faces[iface].end(), std::back_inserter(allFaceNodes));
  }
  POLY_ASSERT(faceNodeCounts.size() == numFaces);

  // Create flags indicating any exterior faces and nodes.
  vector<char> boundaryFaceFlags(numFaces, 0x0);
  for (vector<unsigned>::const_iterator itr = mesh.boundaryFaces.begin();
       itr != mesh.boundaryFaces.end();
       ++itr) 
  {
    POLY_ASSERT(*itr < numFaces);
    boundaryFaceFlags[*itr] = 0x1;
  }

  // Construct the silo cell-face info.  Silo uses the same 1's complement
  // convention polytope does for indicating face orientation, so we can
  // simply copy our faces.
  const int numCells = mesh.cells.size();
  vector<int> cellFaceCounts, allCellFaces;
  cellFaceCounts.reserve(numCells);
  int n;
  for (int icell = 0; icell != numCells; ++icell)
  {
    n = mesh.cells[icell].size();
    cellFaceCounts.push_back(n);
    std::copy(mesh.cells[icell].begin(), mesh.cells[icell].end(), std::back_inserter(allCellFaces));
  }
  POLY_ASSERT(cellFaceCounts.size() == numCells);

  // The polyhedral zone list is referred to in the options list.
  DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"mesh_zonelist");

  // Write out the 3D polyhedral mesh.
  DBPutUcdmesh(file, (char*)"mesh", 3, coordnames, coords,
               numNodes, numCells, 0, 0,
               DB_DOUBLE, optlist); 

  // Write the connectivity information.
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
                  faceNodeCounts.size(), &faceNodeCounts[0], 
                  allFaceNodes.size(), &allFaceNodes[0], 
                  &boundaryFaceFlags[0], 
                  cellFaceCounts.size(), &cellFaceCounts[0],
                  allCellFaces.size(), &allCellFaces[0], 
                  0, 0, numCells-1, optlist);

  // Write out the cell-face connectivity data.
  vector<int> conn(numCells);
  int elemlengths[3];
  char* elemnames[3];
  for (int c = 0; c < numCells; ++c)
    conn[c] = mesh.cells[c].size();
  for (int c = 0; c < numCells; ++c)
  {
    for (int f = 0; f < mesh.cells[c].size(); ++f) {
      int j = mesh.cells[c][f];
      conn.push_back(j < 0 ? ~j : j);
    }
  }
  for (int f = 0; f < mesh.faceCells.size(); ++f)
  {
    conn.push_back(mesh.faceCells[f][0]);
    conn.push_back(mesh.faceCells[f][0]);
  }
  elemnames[0] = strDup("ncellfaces");
  elemlengths[0] = numCells;
  elemnames[2] = strDup("facecells");
  elemlengths[2] = conn.size() - 2*mesh.faces.size();
  elemnames[1] = strDup("cellfaces");
  elemlengths[1] = conn.size() - elemlengths[2] - elemlengths[0];
  DBPutCompoundarray(file, "conn", elemnames, elemlengths, 3, 
                     (void*)&conn[0], conn.size(), DB_INT, 0);
  free(elemnames[0]);
  free(elemnames[1]);
  free(elemnames[2]);

  // Write out convex hull data.
  vector<int> hull(1+mesh.convexHull.facets.size());
  hull[0] = mesh.convexHull.facets.size();
  for (int f = 0; f < mesh.convexHull.facets.size(); ++f)
    hull[1+f] = mesh.convexHull.facets[f].size();
  for (int f = 0; f < mesh.convexHull.facets.size(); ++f)
    for (int n = 0; n < mesh.convexHull.facets[f].size(); ++n)
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
  writeTagsToFile(nodeTags, file, DB_NODECENT);
  writeTagsToFile(edgeTags, file, DB_EDGECENT);
  writeTagsToFile(faceTags, file, DB_FACECENT);
  writeTagsToFile(cellTags, file, DB_ZONECENT);

  // Write out the field mesh data.
  // FIXME: We really should try to use the number of edges for edge fields.
  writeFieldsToFile<RealType>(nodeFields, file, numNodes, DB_NODECENT, optlist);
  writeFieldsToFile<RealType>(edgeFields, file, numFaces, DB_EDGECENT, optlist);
  writeFieldsToFile<RealType>(faceFields, file, numFaces, DB_FACECENT, optlist);
  writeFieldsToFile<RealType>(cellFields, file, numCells, DB_ZONECENT, optlist);

#if 0
  // Vector fields.
  {
    vector<RealType> xdata(mesh.numCells()), ydata(mesh.numCells()),
                 zdata(mesh.numCells());
    char* compNames[3];
    compNames[0] = new char[1024];
    compNames[1] = new char[1024];
    compNames[2] = new char[1024];
    for (map<string, Vector<3>*>::const_iterator iter = vectorFields.begin();
         iter != vectorFields.end(); ++iter)
    {
      snprintf(compNames[0], 1024, "%s_x" , iter->first.c_str());
      snprintf(compNames[1], 1024, "%s_y" , iter->first.c_str());
      snprintf(compNames[2], 1024, "%s_z" , iter->first.c_str());
      Vector<3>* data = iter->second;
      for (int i = 0; i < mesh.numCells(); ++i)
      {
        xdata[i] = data[i][0];
        ydata[i] = data[i][1];
        zdata[i] = data[i][2];
      }
      void* vardata[3];
      vardata[0] = (void*)&xdata[0];
      vardata[1] = (void*)&ydata[0];
      vardata[2] = (void*)&zdata[0];
      DBPutUcdvar(file, (char*)iter->first.c_str(), (char*)"mesh",
                  3, compNames, vardata, mesh.numCells(), 0, 0,
                  DB_DOUBLE, DB_ZONECENT, optlist);
    }
    delete [] compNames[0];
    delete [] compNames[1];
    delete [] compNames[2];
  }
#endif

  // Clean up.
  DBFreeOptlist(optlist);

#ifdef POLYTOPE_ENABLE_MPI
  // Write the multi-block objects to the file if needed.
  int numChunks = nproc / numFiles;
  if (rankInGroup == 0)
  {
    vector<char*> meshNames(numChunks);
    vector<int> meshTypes(numChunks, DB_UCDMESH);
    vector<vector<char*> > varNames(nodeFields.size() +
                                    edgeFields.size() +
                                    faceFields.size() +
                                    cellFields.size());
    vector<int> varTypes(numChunks, DB_UCDVAR);
    for (int i = 0; i < numChunks; ++i)
    {
      // Mesh.
      char meshName[1024];
      snprintf(meshName, 1024, "domain_%d/mesh", i);
      meshNames[i] = strDup(meshName);

      // Field data.
      int fieldIndex = 0;
      appendFieldNames<RealType>(nodeFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(edgeFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(faceFields, fieldIndex, i, varNames);
      appendFieldNames<RealType>(cellFields, fieldIndex, i, varNames);
    }

    // Stick cycle and time in there if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "mesh", numChunks, &meshNames[0], 
                   &meshTypes[0], optlist);
    int fieldIndex = 0;
    putMultivarInFile<RealType>(nodeFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(edgeFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(faceFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);
    putMultivarInFile<RealType>(cellFields, fieldIndex, varNames, varTypes, file, numChunks, optlist);

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < numChunks; ++i)
      free(meshNames[i]);
    for (int f = 0; f < varNames.size(); ++f)
      for (int i = 0; i < numChunks; ++i)
        free(varNames[f][i]);
  }

  // Write the file.
  PMPIO_HandOffBaton(baton, (void*)file);
  PMPIO_Finish(baton);

  // Finally, write the uber-master file.
  if (rank == 0)
  {
    char masterFileName[1024];
    if (cycle >= 0)
      snprintf(masterFileName, 1024, "%s-%d/%s-%d.silo", prefix.c_str(), nproc, prefix.c_str(), cycle);
    else
      snprintf(masterFileName, 1024, "%s-%d/%s.silo", prefix.c_str(), nproc, prefix.c_str());
    int driver = DB_HDF5;
    DBfile* file = DBCreate(masterFileName, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    vector<char*> meshNames(numFiles*numChunks);
    vector<int> meshTypes(numFiles*numChunks, DB_UCDMESH);
    vector<vector<char*> > varNames(nodeFields.size() +
                                    edgeFields.size() +
                                    faceFields.size() +
                                    cellFields.size());
    vector<int> varTypes(numFiles*numChunks, DB_UCDVAR);
    for (int i = 0; i < numFiles; ++i)
    {
      for (int c = 0; c < numChunks; ++c)
      {
        // Mesh.
        char meshName[1024];
        if (cycle >= 0)
          snprintf(meshName, 1024, "%d/%s-%d.silo:/domain_%d/mesh", i, prefix.c_str(), cycle, c);
        else
          snprintf(meshName, 1024, "%d/%s.silo:/domain_%d/mesh", i, prefix.c_str(), c);
        meshNames[i*numChunks+c] = strDup(meshName);

        // Field data.
        int fieldIndex = 0;
        appendFieldNames<RealType>(nodeFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(edgeFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(faceFields, fieldIndex, i, c, cycle, prefix, varNames);
        appendFieldNames<RealType>(cellFields, fieldIndex, i, c, cycle, prefix, varNames);
      }
    }

    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the multimesh and variable data, and close the file.
    DBPutMultimesh(file, "mesh", numFiles*numChunks, &meshNames[0], 
                   &meshTypes[0], optlist);
    int fieldIndex = 0;
    putMultivarInFile<RealType>(nodeFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(edgeFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(faceFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    putMultivarInFile<RealType>(cellFields, fieldIndex, varNames, varTypes, file, numFiles*numChunks, optlist);
    DBClose(file);

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < numFiles*numChunks; ++i)
      free(meshNames[i]);
    for (int f = 0; f < varNames.size(); ++f)
      for (int i = 0; i < numFiles*numChunks; ++i)
        free(varNames[f][i]);
  }
#else
  // Write the file.
  DBClose(file);
#endif
}
//-------------------------------------------------------------------

// Explicit instantiation.
template class SiloWriter<3, double>;

} // end namespace

