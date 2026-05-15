//------------------------------------------------------------------------------
// Take a tessellation and a PLC.  Make sure all points within the given
// tolerance are on the boundary, and that all PLC boundary points are accounted
// for.
//------------------------------------------------------------------------------
#ifndef __polytope_snapToBoundary__
#define __polytope_snapToBoundary__

#include <vector>
#include "Tessellation.hh"
#include "PLC.hh"

namespace polytope {

// 2D
template<typename RealType>
void snapToBoundary(Tessellation<2, RealType>& mesh,
                    const std::vector<RealType>& points,
                    const PLC<2>& geometry,
                    const RealType degeneracy);

// 3D
template<typename RealType>
void snapToBoundary(Tessellation<3, RealType>& mesh,
                    const std::vector<RealType>& points,
                    const PLC<3>& geometry,
                    const RealType degeneracy);

}

#endif
