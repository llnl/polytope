#ifndef __Polytope_FieldData__
#define __Polytope_FieldData__

#include <map>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

namespace polytope {

//! Centering locations for mesh fields.
enum class FieldCentering {
  Node,
  Edge,
  Face,
  Cell
};

inline std::map<std::string, FieldCentering> FieldCenteringMap = {
  {"Node", FieldCentering::Node},
  {"Edge", FieldCentering::Edge},
  {"Face", FieldCentering::Face},
  {"Cell", FieldCentering::Cell}
};

//! Numeric fields associated with a tessellation.
class FieldData {
public:
  using FieldValues = std::variant<std::vector<int>, std::vector<double>>;
  using FieldMap = std::map<std::string, FieldValues>;
  using FieldTypeMap = std::map<FieldCentering, FieldMap>;

  template<typename FieldType>
  using TypedFieldMap = std::map<std::string, std::vector<FieldType>>;

  template<typename FieldType>
  void addField(const FieldCentering centering,
                const TypedFieldMap<FieldType>& newFields) {
    for (const auto& [name, values] : newFields) {
      addField(centering, name, values);
    }
  }

  template<typename FieldType>
  void addField(const FieldCentering centering,
                const std::string& name,
                const std::vector<FieldType>& values) {
    if constexpr (std::is_same_v<FieldType, int> ||
                  std::is_same_v<FieldType, double>) {
      fields[centering].emplace(name, values);
    } else {
      static_assert(std::is_same_v<FieldType, int> ||
                    std::is_same_v<FieldType, double>,
                    "FieldData supports only int and double fields");
    }
  }

  void addDoubleCellField(const std::string& name,
                          const std::vector<double>& values) {
    addField<double>(FieldCentering::Cell, name, values);
  }

  void clear() {
    fields.clear();
  }

  FieldTypeMap fields;
};

} // end namespace polytope

#endif
