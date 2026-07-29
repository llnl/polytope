// A collection of low-level utilities to help with silo file input/output.

#include <cstdlib>
#include <cstring>
#include <iostream>
#include "polytope_internal.hh"
#include "SiloUtils.hh"

using namespace std;

namespace polytope {

char* strDup(const char* s)
{
  if (s == NULL)
    return NULL;
  char* dup = (char*)malloc(sizeof(char) * (strlen(s) + 1));
  strcpy(dup, s);
  return dup;
}

void
writeTagsToFile(const map<string, vector<int>*>& tags,
                DBfile* file,
                int centering)
{
  if (tags.empty()) 
    return;

  // Pack the tags into a compound array.
  vector<int> elemLengths;
  vector<char*> elemNames;
  vector<int> tagData;

  for (map<string, vector<int>*>::const_iterator 
       iter = tags.begin(); iter != tags.end(); ++iter)
  {
    vector<int>& tag = *iter->second;
    elemLengths.push_back(static_cast<int>(tag.size()));
    string tagName;
    elemNames.push_back(strDup(tagName.c_str()));
    for (size_t i = 0; i < tag.size(); ++i)
      tagData.push_back(tag[i]);
  }

  // Write the compound array.
  string tagListName;
  if (centering == DB_NODECENT)
    tagListName = "node_tags";
  else if (centering == DB_EDGECENT)
    tagListName = "edge_tags";
  else if (centering == DB_FACECENT)
    tagListName = "face_tags";
  else
  {
    POLY_ASSERT(centering == DB_ZONECENT);
    tagListName = "cell_tags";
  }
  DBPutCompoundarray(file, tagListName.c_str(), &elemNames[0], &elemLengths[0], 
                     static_cast<int>(tags.size()), (void*)&tagData[0], 
                     static_cast<int>(tagData.size()), DB_INT, 0);

  // Clean up.
  for (size_t i = 0; i < elemNames.size(); ++i)
    free(elemNames[i]);
}

}
