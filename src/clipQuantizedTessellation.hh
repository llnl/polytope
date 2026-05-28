//------------------------------------------------------------------------------
// Clip a QuantizedTessellation against a PLC boundary.
//------------------------------------------------------------------------------
#ifndef __polytope_clipQuantizedTessellation__
#define __polytope_clipQuantizedTessellation__

#include "PLC.hh"

namespace polytope {

template<int Dimension, typename RealType> class Tessellator;
template<typename IntType, typename RealType> class QuantizedTessellation2d;
template<typename IntType, typename RealType> class QuantizedTessellation3d;

// 2D
template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation2d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<2>& geometry,
                               const Tessellator<2, RealType>& tessellator);

// 3D
template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation3d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<3>& geometry,
                               const Tessellator<3, RealType>& tessellator);

}

#endif
