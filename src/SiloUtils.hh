// A collection of low-level utilities to help with silo file input/output.

#ifndef POLYTOPE_SILOUTILS_HH
#define POLYTOPE_SILOUTILS_HH

#include <string>
#include <vector>
#include <map>
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
       ++iter)
  {
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

}
#endif
