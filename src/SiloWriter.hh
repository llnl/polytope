#ifndef __Polytope_SiloWriter__
#define __Polytope_SiloWriter__

#include <string>
#include <float.h>
#include <array>
#include <map>
#include <type_traits>
#include <vector>
#include "Communicator.hh"
#include "SiloUtils.hh"
#include "Tessellation.hh"

namespace polytope {

template <int Dimension, typename TessType>
class SiloWriter;

//! \class SiloWriter
//! This class writes a tessellation and its associated fields to Silo files.
template <int Dimension, typename TessType>
class SiloWriterBase {
public:
  // Map of a field name and it's values
  template<typename FieldType>
  using FieldMap = std::map<std::string, std::vector<FieldType>>;
  // Map of a Polytope field centering and its FieldMap
  template<typename FieldType>
  using FieldTypeMap = std::map<FieldCentering, FieldMap<FieldType>>;

  SiloWriterBase(const TessType& mesh) :
    m_mesh(mesh) { }
  virtual ~SiloWriterBase() = default;

  // Routines for filling field data
  template<typename FieldType,
           std::enable_if_t<std::is_same_v<FieldType, double> ||
                            std::is_same_v<FieldType, int>, int> = 0>
  void addField(const FieldCentering& centering,
                const FieldMap<FieldType>& field) {
    if constexpr (std::is_floating_point_v<FieldType>) {
      m_doubleFields[centering].insert(field.begin(), field.end());
    } else {
      m_intFields[centering].insert(field.begin(), field.end());
    }
  }

  template<typename FieldType,
           std::enable_if_t<std::is_same_v<FieldType, double> ||
                            std::is_same_v<FieldType, int>, int> = 0>
  void addField(const FieldCentering& centering,
                const std::string& name,
                const std::vector<FieldType>& vals) {
    addField<FieldType>(centering, FieldMap<FieldType>{{name, vals}});
  }

  // template<typename FieldType,
  //          std::enable_if_t<std::is_same_v<FieldType, double> ||
  //                           std::is_same_v<FieldType, int>, int> = 0>
  // void addField(const FieldCentering& centering,
  //               const std::string& name,
  //               const std::vector<FieldType>&& vals) {
  //   addField<FieldType>(centering, FieldMap<FieldType>{{name, vals}});
  // }

  void addDoubleCellField(const std::string& name,
                          const std::vector<double>& vals) {
    addField<double>(FieldCentering::Cell, name, vals);
  }

  virtual void write(const std::string& filePrefix,
                     const std::string& directory,
                     int cycle,
                     double time,
                     int numFiles = -1) = 0;

  // This version omits the cycle and time arguments
  void write(const std::string& filePrefix,
             const std::string& directory,
             int numFiles = -1) {
    write(filePrefix, directory, -1, -1., numFiles);
  }

  // This version omits the directory
  void write(const std::string& filePrefix,
             int cycle,
             double time,
             int numFiles = -1) {
    write(filePrefix, "", cycle, time, numFiles);
  }

  // This version omits the cycle, time, and directory arguments
  void write(const std::string& filePrefix,
             int numFiles = -1) {
    write(filePrefix, "", numFiles);
  }

  // Generate variables used in testing
  void generateTestVars() {
    const auto numCells = m_mesh.cells.size();
    std::vector<int> indices(numCells);
    std::array<std::vector<double>, Dimension> pos;
    for (auto& component : pos) component.resize(numCells);
    for (auto i = 0u; i < numCells; ++i) {
      for (int d = 0; d < Dimension; ++d) {
        pos[d][i] = m_mesh.points[i][d];
      }
      indices[i] = int(i);
    }
    std::array<std::string, 3> pos_names = {"x", "y", "z"};
    FieldMap<double> cellDoubleFields;
    for (int d = 0; d < Dimension; ++d) {
      cellDoubleFields[pos_names[d]] = pos[d];
    }
    FieldMap<int> cellIntFields;
    cellIntFields["cell_index"] = indices;
#ifdef POLYTOPE_ENABLE_MPI
    if (!m_mesh.cellRank.empty()) {
      std::vector<int> ranks = m_mesh.cellRank;
      cellIntFields["rank"] = ranks;
    } else {
      int rank = Communicator::getRank();
      std::vector<int> ranks(numCells, rank);
      cellIntFields["rank"] = ranks;
    }
#endif
    addField<double>(FieldCentering::Cell, cellDoubleFields);
    addField<int>(FieldCentering::Cell, cellIntFields);
  }

  // Whether to write Overlink file type
  bool m_overlinkType = false;
  void writeOvlType(const bool inBool) { m_overlinkType = inBool; }
  bool writeType() { return m_overlinkType; }
  FieldTypeMap<double> m_doubleFields;
  FieldTypeMap<int> m_intFields;
  const TessType& m_mesh;
protected:
  void writeFieldsToFile(const std::string& meshname,
                         DBfile* file,
                         DBoptlist* optlist) {
    fieldWrite<double>(m_doubleFields, meshname, file, optlist);
    fieldWrite<int>(m_intFields, meshname, file, optlist);
  }

  template<typename FieldType>
  void fieldWrite(const FieldTypeMap<FieldType>& fields,
                  const std::string& meshname,
                  DBfile* file,
                  DBoptlist* optlist) {
    const auto numFaces = m_mesh.faces.size();
    const auto numNodes = m_mesh.nodes.size();
    const auto numCells = m_mesh.cells.size();
    for (const auto& [centering, fieldmap] : fields) {
      const int siloType = static_cast<int>(centering);
      if (centering == FieldCentering::Face ||
          centering == FieldCentering::Edge) {
        writeFields<FieldType>(fieldmap, meshname, file, numFaces, siloType, optlist);
      } else if (centering == FieldCentering::Node) {
        writeFields<FieldType>(fieldmap, meshname, file, numNodes, siloType, optlist);
      } else {
        writeFields<FieldType>(fieldmap, meshname, file, numCells, siloType, optlist);
      }
    }
  }
};

//! Instance writer for two-dimensional tessellations.
template <typename TessType>
class SiloWriter<2, TessType>: public SiloWriterBase<2, TessType> {
  using Base = SiloWriterBase<2, TessType>;
protected:
  using Base::writeFieldsToFile;
public:
  using Base::Base;
  using Base::write;
  using Base::m_doubleFields;
  using Base::m_intFields;
  using Base::m_mesh;

  void write(const std::string& filePrefix,
             const std::string& directory,
             int cycle,
             double time,
             int numFiles = -1) override;
};

//! Instance writer for three-dimensional tessellations.
template <typename TessType>
class SiloWriter<3, TessType>: public SiloWriterBase<3, TessType> {
  using Base = SiloWriterBase<3, TessType>;
protected:
  using Base::writeFieldsToFile;
public:
  using Base::Base;
  using Base::write;
  using Base::m_doubleFields;
  using Base::m_intFields;
  using Base::m_mesh;

  void write(const std::string& filePrefix,
             const std::string& directory,
             int cycle,
             double time,
             int numFiles = -1) override;
};

} // end namespace

#endif
