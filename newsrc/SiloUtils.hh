// A collection of low-level utilities to help with silo file input/output.

#ifndef POLYTOPE_SILOUTILS_HH
#define POLYTOPE_SILOUTILS_HH

#include <string>
#include <vector>
#include <map>
#include <cstdio>
#include "silo.h"

namespace polytope {

// strdup isn't part of the C standard, so we can't rely on its existence.
// We keep our own handy.
char* strDup(const char* s);

void
writeTagsToFile(const std::map<std::string, std::vector<int>*>& tags,
                DBfile* file,
                int centering);

template <typename RealType>
void
writeFieldsToFile(const std::map<std::string, RealType*>& fields,
                  DBfile* file,
                  const int numElements,
                  const int centering,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
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
template <typename RealType>
void
writeFieldsToFile(const std::map<std::string, RealType*>& fields,
                  const std::string& meshname,
                  DBfile* file,
                  const int numElements,
                  const int centering,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
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
template <typename RealType>
void
appendFieldNames(const std::map<std::string, RealType*>& fields,
                 int& fieldIndex,
                 const int ichunk,
                 std::vector<std::vector<char*> >& varNames) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
       iter != fields.end();
       ++iter, ++fieldIndex)
  {
    char varName[1024];
    snprintf(varName, 1024, "domain_%d/%s", ichunk, iter->first.c_str());
    varNames[fieldIndex].push_back(strDup(varName));
  }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename RealType>
void
appendFieldNames(const std::map<std::string, RealType*>& fields,
                 int& fieldIndex,
                 const int ifile,
                 const int ichunk,
                 const int cycle,
                 const std::string& prefix,
                 std::vector<std::vector<char*> >& varNames) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
       iter != fields.end(); ++iter, ++fieldIndex)
  {
    char varName[1024];
    if (cycle >= 0)
      snprintf(varName, 1024, "%d/%s-%d.silo:/domain_%d/%s", ifile, prefix.c_str(), cycle, ichunk, iter->first.c_str());
    else
      snprintf(varName, 1024, "%d/%s.silo:/domain_%d/%s", ifile, prefix.c_str(), ichunk, iter->first.c_str());
    varNames[fieldIndex].push_back(strDup(varName));
  }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename RealType>
void
putMultivarInFile(const std::map<std::string, RealType*>& fields,
                  int& fieldIndex,
                  std::vector<std::vector<char*> >& varNames,
                  std::vector<int>& varTypes,
                  DBfile* file,
                  const int numChunks,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
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
std::string getFileName(const std::string& directory,
                        const std::string& prefix,
                        const int groupRank) {
  char groupchar[256];
  sprintf(groupchar, "_%04d", groupRank);
  return directory + "/" + prefix + std::string(groupchar) + ".silo";
}

inline
std::string getRankDir(const int rankInGroup) {
  char rgchar[256];
  sprintf(rgchar, "_%04d", rankInGroup);
  return "domain" + std::string(rgchar);
}

inline
std::string getMasterFileName(const std::string& prefix,
                              const int& cycle = -1) {
  if (cycle >= 0) {
    char cyclechar[256];
    sprintf(cyclechar, "_%04d", cycle);
    return prefix + std::string(cyclechar) + ".silo";
  } else {
    return prefix + ".silo";
  }
}

inline
std::vector<std::string> getProcPaths(const std::string& directory,
                                      const std::string& prefix,
                                      const size_t nproc,
                                      const std::vector<std::pair<int, int>>& pgr) {
  std::vector<std::string> out;
  for (int p = 0; p < nproc; ++p) {
    int groupRank = pgr[p].first;
    int rankInGroup = pgr[p].second;
    auto filename = getFileName(directory, prefix, groupRank);
    auto rankdir = getRankDir(rankInGroup);
    out.push_back(filename + ":" + rankdir + "/");
  }
  return out;
}

template<typename RealType>
void putCellVars(DBfile* file,
                 const std::map<std::string, RealType*>& fields,
                 const std::vector<std::string>& procPaths,
                 const size_t nproc,
                 const std::vector<int>& varTypes,
                 DBoptlist* optlist) {
  for (auto iter = fields.begin(); iter != fields.end(); ++iter) {
    std::vector<char*> varNames;
    for (const auto& procPath : procPaths) {
      auto varName = procPath + iter->first;
      varNames.push_back(strDup(varName.c_str()));
    }
    auto gvarName = iter->first;
    DBPutMultivar(file, gvarName.c_str(), nproc, varNames.data(), varTypes.data(), optlist);
    for (auto f = 0; f < varNames.size(); ++f) {
      free(varNames[f]);
    }
  }
}

//-------------------------------------------------------------------
// Build mesh name paths for multi-block aggregation
//-------------------------------------------------------------------
void
buildMeshNames(std::vector<char*>& meshNames,
               const std::string& subdir,
               const std::string& meshName,
               const int numChunks,
               const int startChunk = 0);

//-------------------------------------------------------------------
// Build mesh name paths for master file aggregation
//-------------------------------------------------------------------
void
buildMasterMeshNames(std::vector<char*>& meshNames,
                     const std::string& subdir,
                     const std::string& meshName,
                     const int numFiles,
                     const int numChunks,
                     const int cycle,
                     const std::string& prefix);

//-------------------------------------------------------------------
// Extended version of appendFieldNames with subdirectory support
//-------------------------------------------------------------------
template <typename RealType>
void
appendFieldNamesWithSubdir(const std::map<std::string, RealType*>& fields,
                           int& fieldIndex,
                           const int ichunk,
                           const std::string& subdir,
                           std::vector<std::vector<char*> >& varNames) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
       iter != fields.end();
       ++iter, ++fieldIndex)
  {
    char varName[1024];
    snprintf(varName, 1024, "domain_%d/%s/%s", ichunk, subdir.c_str(), iter->first.c_str());
    varNames[fieldIndex].push_back(strDup(varName));
  }
}

//-------------------------------------------------------------------
// Extended version for master file with subdirectory support
//-------------------------------------------------------------------
template <typename RealType>
void
appendFieldNamesWithSubdir(const std::map<std::string, RealType*>& fields,
                           int& fieldIndex,
                           const int ifile,
                           const int ichunk,
                           const int cycle,
                           const std::string& prefix,
                           const std::string& subdir,
                           std::vector<std::vector<char*> >& varNames) {
  for (typename std::map<std::string, RealType*>::const_iterator iter = fields.begin();
       iter != fields.end(); ++iter, ++fieldIndex)
  {
    char varName[1024];
    if (cycle >= 0)
      snprintf(varName, 1024, "%d/%s-%d.silo:/domain_%d/%s/%s",
               ifile, prefix.c_str(), cycle, ichunk, subdir.c_str(), iter->first.c_str());
    else
      snprintf(varName, 1024, "%d/%s.silo:/domain_%d/%s/%s",
               ifile, prefix.c_str(), ichunk, subdir.c_str(), iter->first.c_str());
    varNames[fieldIndex].push_back(strDup(varName));
  }
}
//-------------------------------------------------------------------

}
#endif
