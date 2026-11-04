// A collection of low-level utilities to help with silo file input/output.

#ifndef __Polytope_SiloUtils__
#define __Polytope_SiloUtils__

#include <string>
#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include "silo.h"
#include "Communicator.hh"
#include "FieldData.hh"
#include "polytope_internal.hh"

namespace polytope {

//! Translate mesh field centering to Silo's centering constants.
inline int siloFieldCentering(const FieldCentering centering) {
  switch (centering) {
    case FieldCentering::Node: return DB_NODECENT;
    case FieldCentering::Edge: return DB_EDGECENT;
    case FieldCentering::Face: return DB_FACECENT;
    case FieldCentering::Cell: return DB_ZONECENT;
  }
  throw std::invalid_argument("Invalid field centering");
}

inline FieldCentering fieldCenteringFromSilo(const int centering) {
  switch (centering) {
    case DB_NODECENT: return FieldCentering::Node;
    case DB_EDGECENT: return FieldCentering::Edge;
    case DB_FACECENT: return FieldCentering::Face;
    case DB_ZONECENT: return FieldCentering::Cell;
  }
  throw std::invalid_argument("Invalid Silo field centering");
}

inline void storeReadField(FieldData& fieldData,
                           const FieldCentering centering,
                           const std::string& name,
                           const DBucdvar& dbvar) {
  if (dbvar.datatype == DB_INT) {
    std::vector<int> values;
    if (dbvar.nels > 0) {
      const auto* data = static_cast<const int*>(dbvar.vals[0]);
      values.assign(data, data + dbvar.nels);
    }
    fieldData.fields[centering].insert_or_assign(name, std::move(values));
  } else if (dbvar.datatype == DB_DOUBLE) {
    std::vector<double> values;
    if (dbvar.nels > 0) {
      const auto* data = static_cast<const double*>(dbvar.vals[0]);
      values.assign(data, data + dbvar.nels);
    }
    fieldData.fields[centering].insert_or_assign(name, std::move(values));
  } else {
    throw std::invalid_argument("Unsupported Silo field data type");
  }
}

//! Silo data types
enum class FieldDataType {
  Double = DB_DOUBLE,
  Int = DB_INT
};

//! Overlink attributes
static std::map<std::string, int> OverlinkAttr = {
    {"ATTR_NODAL",        0},
    {"ATTR_ZONAL",        1},
    {"ATTR_FACE",         2},
    {"ATTR_EDGE",         3},
    {"ATTR_INTENSIVE",    0},
    {"ATTR_EXTENSIVE",    1},
    {"ATTR_FIRST_ORDER",  0},
    {"ATTR_SECOND_ORDER", 1},
    {"ATTR_INTEGER",      0},
    {"ATTR_FLOAT",        1}
};

// strdup isn't part of the C standard, so we can't rely on its existence.
// We keep our own handy.
inline char* strDup(const char* s) {
  if (s == NULL)
    return NULL;
  char* dup = (char*)malloc(sizeof(char) * (strlen(s) + 1));
  strcpy(dup, s);
  return dup;
}

// Certain names
inline std::string getLocalMeshName() {
  return "MESH";
}

inline std::string getGlobalMeshName() {
  return "MMESH";
}

inline std::string getLocalMatName() {
  return "MATERIAL";
}

inline std::string getGlobalMatName() {
  return "MMATERIAL";
}

// Name of custom variables, tessellator name, coord bit count
inline std::string getTessNameVar() {
  return "tessellatorName";
}

inline std::string getCoordBitVar() {
  return "coord_bits";
}

inline std::string getEscapePodVar() {
  return "escape_pod";
}

template<typename FieldType>
void
writeField(const std::string& name,
           const std::vector<FieldType>& values,
           const std::string& meshname,
           DBfile* file,
           const int numElements,
           const int centering,
           DBoptlist* optlist) {
  int siloDataType = static_cast<int>(FieldDataType::Int);
  if constexpr (std::is_floating_point_v<FieldType>) {
    siloDataType = static_cast<int>(FieldDataType::Double);
  }
  POLY_CHECK(DBPutUcdvar1(file,
                          const_cast<char*>(name.c_str()),
                          meshname.c_str(),
                          const_cast<FieldType*>(values.data()),
                          numElements,
                          0,
                          0,
                          siloDataType,
                          centering,
                          optlist) != 1);
}

template<typename FieldType>
void writeFields(const std::map<std::string, std::vector<FieldType>>& fields,
                 const std::string& meshname,
                 DBfile* file,
                 const int numElements,
                 const int centering,
                 DBoptlist* optlist) {
  for (const auto& [name, values] : fields) {
    writeField(name, values, meshname, file, numElements, centering, optlist);
  }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------

inline
void
putMultivarInFile(const std::map<std::string, std::vector<double>>& fields,
                  int& fieldIndex,
                  std::vector<std::vector<char*> >& varNames,
                  std::vector<int>& varTypes,
                  DBfile* file,
                  const int numChunks,
                  DBoptlist* optlist) {
  for (typename std::map<std::string, std::vector<double>>::const_iterator iter = fields.begin();
       iter != fields.end();
       ++iter, ++fieldIndex) {
    POLY_CHECK(DBPutMultivar(file, iter->first.c_str(), numChunks,
                             &varNames[fieldIndex][0], &varTypes[0], optlist) != 1);
  }
}

//-------------------------------------------------------------------
inline
std::string getMasterDirname(const std::string& directory,
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
  return directory + "/domain" + std::to_string(rank) + ".silo";
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
  int nranks = Communicator::getNRanks();
  int rank = Communicator::getRank();
  int root = Communicator::getRoot();
  auto& comm = Communicator::communicator();
  std::vector<int> flags;
  if (rank == root) {
    flags.resize(nranks);
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

template<typename FieldMapType>
void putCellVars(DBfile* file,
                 const FieldMapType& fields,
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
    POLY_CHECK(DBPutMultivar(file, gvarName.c_str(), nblocks,
                             varNames.data(), varTypes.data(), optlist) != 1);
    for (auto f = 0u; f < varNames.size(); ++f) {
      free(varNames[f]);
    }
  }
}

//-------------------------------------------------------------------

}
#endif
