// A collection of low-level utilities to help with silo file input/output.

#ifndef __Polytope_SiloUtils__
#define __Polytope_SiloUtils__

#include <string>
#include <vector>
#include <map>
#include <cstdio>
#include "silo.h"
#include "Communicator.hh"

namespace polytope {

// strdup isn't part of the C standard, so we can't rely on its existence.
// We keep our own handy.
char* strDup(const char* s);

// Certain names
inline std::string getLocalMeshName() {
  return "MESH";
}

inline std::string getGlobalMeshName() {
  return "MMESH";
}

void
writeTagsToFile(const std::map<std::string, std::vector<int>*>& tags,
                DBfile* file,
                int centering);

inline
void
writeFieldsToFile(const std::map<std::string, double*>& fields,
                  DBfile* file,
                  const int numElements,
                  const int centering,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, double*>::const_iterator iter = fields.begin();
       iter != fields.end();
       ++iter) {
    DBPutUcdvar1(file,
                 (char*)iter->first.c_str(),
                 (char*)"mesh",
                 (void*)iter->second,
                 numElements,
                 0,
                 0,
                 DB_DOUBLE,
                 centering,
                 optlist);
  }
}

inline
void
writeFieldsToFile(const std::map<std::string, double*>& fields,
                  const std::string& meshname,
                  DBfile* file,
                  const int numElements,
                  const int centering,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, double*>::const_iterator iter = fields.begin();
       iter != fields.end();
       ++iter) {
    DBPutUcdvar1(file,
                 (char*)iter->first.c_str(),
                 meshname.c_str(),
                 (void*)iter->second,
                 numElements,
                 0,
                 0,
                 DB_DOUBLE,
                 centering,
                 optlist);
  }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------

inline
void
putMultivarInFile(const std::map<std::string, double*>& fields,
                  int& fieldIndex,
                  std::vector<std::vector<char*> >& varNames,
                  std::vector<int>& varTypes,
                  DBfile* file,
                  const int numChunks,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, double*>::const_iterator iter = fields.begin();
       iter != fields.end();
       ++iter, ++fieldIndex)
  {
    DBPutMultivar(file, iter->first.c_str(), numChunks,
                  &varNames[fieldIndex][0], &varTypes[0], optlist);
  }
}

//-------------------------------------------------------------------
inline
std::string getMasterDirName(const std::string& directory,
                             const std::string& prefix,
                             const int cycle) {
  std::string outdir = directory;
  if (directory.empty()) {
    outdir = prefix;
  }
  if (cycle >= 0) {
    char dirchar[256];
    sprintf(dirchar, "_%04d", cycle);
    return outdir + dirchar;
  } else {
    return outdir;
  }
}

inline
std::string getFilename(const std::string& directory,
                        const int rank) {
  return directory + "/domain_" + std::to_string(rank) + ".silo";
}

inline
std::string getMasterFilename(const std::string& prefix,
                              const int& cycle = -1) {
  if (cycle >= 0) {
    char cyclechar[256];
    sprintf(cyclechar, "_%04d", cycle);
    return prefix + std::string(cyclechar) + ".silo";
  } else {
    return prefix + ".silo";
  }
}

#ifdef POLYTOPE_ENABLE_MPI
// Gather all ranks that have valid tessellation data on them
inline
std::vector<int> gatherValidRanks(int hasData) {
  int nproc = Communicator::getNProcs();
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  auto& comm = Communicator::communicator();
  std::vector<int> flags;
  if (rank == root) {
    flags.resize(nproc);
  }
  MPI_Gather(&hasData, 1, MPI_INT,
             rank == root ? flags.data() : nullptr,
             1, MPI_INT, root, comm);
  std::vector<int> ranksWithData;
  if (rank == root) {
    int k = 0;
    for (const auto& f : flags) {
      if (f) {
        ranksWithData.push_back(k);
      }
      k++;
    }
  }
  return ranksWithData;
}
#endif

inline
std::vector<std::string> getProcPaths(const std::string& directory,
                                      const std::vector<int> ranksWithData) {
  std::vector<std::string> out;
  for (const auto& p : ranksWithData) {
    auto filename = getFilename(directory, p) + ":";
    out.push_back(filename);
  }
  return out;
}

inline
void putCellVars(DBfile* file,
                 const std::map<std::string, double*>& fields,
                 const std::vector<std::string>& procPaths,
                 const size_t nblocks,
                 const std::vector<int>& varTypes,
                 DBoptlist* optlist) {
  for (auto iter = fields.begin(); iter != fields.end(); ++iter) {
    std::vector<char*> varNames;
    for (const auto& procPath : procPaths) {
      auto varName = procPath + iter->first;
      varNames.push_back(strDup(varName.c_str()));
    }
    auto gvarName = iter->first;
    DBPutMultivar(file, gvarName.c_str(), nblocks, varNames.data(), varTypes.data(), optlist);
    for (auto f = 0u; f < varNames.size(); ++f) {
      free(varNames[f]);
    }
  }
}

//-------------------------------------------------------------------

}
#endif
