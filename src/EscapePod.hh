//-----------------------------------------------------------------------------//
// EscapePod
//
// Helpers for QuantTessellation escape pod debug files.
//-----------------------------------------------------------------------------//

#ifndef __Polytope_EscapePod__
#define __Polytope_EscapePod__

#include "Communicator.hh"
#include "Quantizer.hh"

#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace polytope {
namespace escape_pod {

inline
std::string
filename(const std::string& baseName) {
  std::ostringstream os;
  os << baseName << "." << Communicator::getRank();
  return os.str();
}

inline
void
expectToken(std::istream& in,
            const std::string& expected) {
  std::string token;
  in >> token;
  POLY_VERIFY2(in and token == expected,
               "Malformed escape pod file: expected token '" << expected
               << "' but found '" << token << "'");
}

template<typename TessellationType>
inline
void
rebuildTessellationPointMetadata(TessellationType& tessellation) {
  const auto& Q = Quantizer<2>::instance();
  tessellation.hashes.clear();
  tessellation.hashes.reserve(tessellation.points.size());
  for (unsigned i = 0; i < tessellation.points.size(); ++i) {
    auto& point = tessellation.points[i];
    point.index = i;
    tessellation.hashes.push_back(Q.encode(point));
  }
}

template<typename QPLCType>
inline
void
rebuildQPLCPointMetadata(QPLCType& qplc) {
  const auto& Q = Quantizer<2>::instance();
  qplc.m_loBounds = Q.maxBound;
  qplc.m_hiBounds = -qplc.m_loBounds;
  for (unsigned i = 0; i < qplc.points.size(); ++i) {
    auto& point = qplc.points[i];
    point.index = i;
    qplc.m_loBounds = qplc.m_loBounds.minElements(point);
    qplc.m_hiBounds = qplc.m_hiBounds.maxElements(point);
  }
}

}
}

#endif
