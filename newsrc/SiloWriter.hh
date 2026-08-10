#ifndef __Polytope_SiloWriter__
#define __Polytope_SiloWriter__

#include <string>
#include <float.h>
#include <map>
#include <vector>
#include "Communicator.hh"
#include "silo.h"

namespace polytope {

template<int Dimension, typename RealType> class Tessellation;

//! \class SiloWriter
//! This class provides a static interface for writing Silo files
//! containing tessellations made by polytope.
template <int Dimension, typename TessType>
class SiloWriter {
  // No general recipe
};

//! Partial specialization for 2D tessellations.
template <typename TessType>
class SiloWriter<2, TessType> {
public:
  // Map of a field name and it's values
  using FieldMap = std::map<std::string, std::vector<double>>;
  using TagMap = std::map<std::string, std::vector<int>*>;
  // Map of the centering type (DB_ZONECENT, DB_EDGECENT, ect) and it's FieldMap
  using FieldTypeMap = std::map<int, FieldMap>;
  using TagTypeMap = std::map<int, TagMap>;

  //! Write an arbitrary polygonal mesh, an associated set of
  //! (node, edge, face, cell)-centered fields, and a corresponding set of
  //! tags, to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const TagTypeMap& tags,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
                    int numFiles = -1);

  //! Write an arbitrary polygonal mesh and an associated set of
  //! (node, edge, face, cell)-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
                    int numFiles = -1) {
    // Just call the general function with no tags.
    TagTypeMap tags;
    write(mesh, fields, tags, filePrefix, directory, cycle, time, numFiles);
  }

  //! Write an arbitrary polygonal mesh and an associated set of
  //! (node, edge, face, cell)-centered fields to a SILO file. This version generates a
  //! directory name automatically. For parallel runs, the directory
  //! name is filePrefix-nproc. For serial runs, the directory is
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    int cycle,
                    double time,
                    int numFiles = -1) {
    write(mesh, fields, filePrefix, "", cycle, time, numFiles);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int numFiles = -1) {
    write(mesh, fields,filePrefix, directory, -1, -1., numFiles);
  }

  //! This version of write omits the cycle and time arguments and
  //! automatically generates the directory name from the file prefix.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    int numFiles = -1) {
    write(mesh, fields, filePrefix, "", -1, -1., numFiles);
  }

  static void write(const TessType& mesh,
                    const std::string& filePrefix,
                    int numFiles = -1,
                    int cycle = 0,
                    double time = 0) {
    size_t meshSize = mesh.cells.size();
    std::vector<double> index(meshSize);
    std::vector<double> genx (meshSize);
    std::vector<double> geny (meshSize);
    FieldMap cellFields;
#ifdef POLYTOPE_ENABLE_MPI
    std::vector<double> rankField;
#endif
    for (auto i = 0u; i < meshSize; ++i) {
      index[i] = double(i);
      genx[i] = static_cast<double>(mesh.points[i].x);
      geny[i] = static_cast<double>(mesh.points[i].y);
    }
    cellFields["cell_index"] = index;
    cellFields["gen_x"     ] = genx;
    cellFields["gen_y"     ] = geny;
#ifdef POLYTOPE_ENABLE_MPI
    if (mesh.cellRank.size() > 0) {
      rankField.resize(meshSize);
      for (auto i = 0u; i < meshSize; ++i) {
        rankField[i] = static_cast<double>(mesh.cellRank[i]);
      }
      cellFields["rank"    ] = rankField;
    } else {
      int rank = Communicator::getRank();
      rankField.assign(meshSize, static_cast<double>(rank));
      cellFields["rank"    ] = rankField;
    }
#endif
    FieldTypeMap fields;
    fields[DB_ZONECENT] = cellFields;
    write(mesh, fields, filePrefix, "", cycle, time, numFiles);
  }

};

//! Partial specialization for 3D tessellations.
template <typename TessType>
class SiloWriter<3, TessType> {
public:
  // Map of a field name and it's values
  using FieldMap = std::map<std::string, std::vector<double>>;
  using TagMap = std::map<std::string, std::vector<int>*>;
  // Map of the centering type (DB_ZONECENT, DB_EDGECENT, ect) and it's FieldMap
  using FieldTypeMap = std::map<int, FieldMap>;
  using TagTypeMap = std::map<int, TagMap>;

  //! Write an arbitrary polygonal mesh, an associated set of
  //! (node, edge, face, cell)-centered fields, and a corresponding set of
  //! tags, to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const TagTypeMap& tags,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
                    int numFiles = -1);

  //! Write an arbitrary polyhedral mesh and an associated set of
  //! (node, edge, face, cell)-centered fields to a SILO file in the given directory.
  //! \param numFiles The number of files that will be written. If this
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int cycle,
                    double time,
                    int numFiles = -1) {
    // Just call the general function with no tags.
    TagTypeMap tags;
    write(mesh, fields, tags, filePrefix, directory, cycle, time, numFiles);
  }

  //! Write an arbitrary polyhedral mesh and an associated set of
  //! (node, edge, face, cell)-centered fields to a SILO file. This version generates a
  //! directory name automatically. For parallel runs, the directory
  //! name is filePrefix-nproc. For serial runs, the directory is
  //! the current working directory.
  //! \param numFiles The number of files that will be written. If this
  //!                 is set to -1, one file will be written for each process.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    int cycle,
                    double time,
                    int numFiles = -1)
  {
    write(mesh, fields, filePrefix, "", cycle, time, numFiles);
  }

  //! This version of write omits the cycle and time arguments.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    const std::string& directory,
                    int numFiles = -1)
  {
    write(mesh, fields, filePrefix, directory, -1, -1., numFiles);
  }

  //! This version of write omits the cycle and time arguments and
  //! automatically generates the directory name from the file prefix.
  static void write(const TessType& mesh,
                    const FieldTypeMap& fields,
                    const std::string& filePrefix,
                    int numFiles = -1)
  {
    write(mesh, fields, filePrefix, "", -1, -1., numFiles);
  }

  static void write(const TessType& mesh,
                    const std::string& filePrefix,
                    int numFiles = -1) {
    size_t meshSize = mesh.cells.size();
    std::vector<double> index(meshSize);
    std::vector<double> genx (meshSize);
    std::vector<double> geny (meshSize);
    std::vector<double> genz (meshSize);
    FieldMap cellFields;
#ifdef POLYTOPE_ENABLE_MPI
    std::vector<double> rankField;
#endif
    for (auto i = 0u; i < meshSize; ++i) {
      index[i] = double(i);
      genx[i] = static_cast<double>(mesh.points[i].x);
      geny[i] = static_cast<double>(mesh.points[i].y);
      genz[i] = static_cast<double>(mesh.points[i].z);
    }
    cellFields["cell_index"] = index;
    cellFields["gen_x"     ] = genx;
    cellFields["gen_y"     ] = geny;
    cellFields["gen_z"     ] = genz;
#ifdef POLYTOPE_ENABLE_MPI
    if (mesh.cellRank.size() > 0) {
      rankField.resize(meshSize);
      for (auto i = 0u; i < meshSize; ++i) {
        rankField[i] = static_cast<double>(mesh.cellRank[i]);
      }
      cellFields["rank"    ] = rankField;
    } else {
      int rank = Communicator::getRank();
      rankField.assign(meshSize, static_cast<double>(rank));
      cellFields["rank"      ] = rankField;
    }
#endif
    FieldTypeMap fields;
    fields[DB_ZONECENT] = cellFields;
    write(mesh, fields, filePrefix, "", -1, -1., numFiles);
  }

};

} // end namespace

#endif
