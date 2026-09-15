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

template<int Dimension, typename TessType>
class SiloWriter;

//! \class SiloWriter
//! This class writes a tessellation and its associated fields to Silo files.
template<int Dimension, typename TessType>
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
  void generateTestVars();

  // Whether to write Overlink file type
  bool m_overlinkType = false;
  void setOvlType(const bool inBool) { m_overlinkType = inBool; }
  bool getOvlType() { return m_overlinkType; }
  FieldTypeMap<double> m_doubleFields;
  FieldTypeMap<int> m_intFields;
  const TessType& m_mesh;

protected:

  // Write the double and int fields
  void writeFieldsToFile(const std::string& meshname,
                         DBfile* file,
                         DBoptlist* optlist) {
    fieldWrite<double>(m_doubleFields, meshname, file, optlist);
    fieldWrite<int>(m_intFields, meshname, file, optlist);
  }

  // Determine the centering and write the appropriate field
  template<typename FieldType>
  void fieldWrite(const FieldTypeMap<FieldType>& fields,
                  const std::string& meshname,
                  DBfile* file,
                  DBoptlist* optlist);

  void writeMasterFile(const std::string& masterFilename,
                       const std::string& masterDirname,
                       const std::vector<int>& ranksWithData,
                       double& time,
                       int cycle);
};

//! Instance writer for two-dimensional tessellations.
template<typename TessType>
class SiloWriter<2, TessType>: public SiloWriterBase<2, TessType> {
  using Base = SiloWriterBase<2, TessType>;
protected:
  using Base::writeFieldsToFile;
  using Base::writeMasterFile;
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
template<typename TessType>
class SiloWriter<3, TessType>: public SiloWriterBase<3, TessType> {
  using Base = SiloWriterBase<3, TessType>;
protected:
  using Base::writeFieldsToFile;
  using Base::writeMasterFile;
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

#include "SiloWriterInline.hh"

#endif
