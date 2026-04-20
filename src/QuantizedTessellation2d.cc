//------------------------------------------------------------------------
// QuantizedTessellation2d
//
// An intermediate representation for 2D tessellations in integer
// coordinates.
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "QuantizedTessellation2d.hh"
#include "polytope.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {
  
//------------------------------------------------------------------------------
// Construct with the given generators.  Finds the bounding limits and sets the
// quantized generators.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation2d<IntType, RealType>::
QuantizedTessellation2d(const std::vector<RealType>& points,
                        const std::vector<RealType>& boundaryPoints) {
  geometry::computeBoundingBox<2, RealType>(&boundaryPoints[0], boundaryPoints.size(), true, xmin, xmax);
  length = std::max(xmax[0] - xmin[0], xmax[1] - xmin[1]);
  xmin[0] -= 2.0*length;
  xmin[1] -= 2.0*length;
  xmax[0] += 2.0*length;
  xmax[1] += 2.0*length;
  length *= 5.0;
  this->construct(points);
}

//------------------------------------------------------------------------------
// Construct with the given generators using the specified bounds.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation2d<IntType, RealType>::
QuantizedTessellation2d(const std::vector<RealType>& points,
                        const RealType xmin_in[2],
                        const RealType xmax_in[2]) {
  xmin[0] = xmin_in[0];
  xmin[1] = xmin_in[1];
  xmax[0] = xmax_in[0];
  xmax[1] = xmax_in[1];
  length = std::max(xmax[0] - xmin[0], xmax[1] - xmin[1]);
  this->construct(points);
}

//------------------------------------------------------------------------------
// Special copy constructor -- take all the state of other QuantizedTessellation
// except use the given integer generators.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation2d<IntType, RealType>::
QuantizedTessellation2d(const std::vector<IntPoint>& generators,
                        const QuantizedTessellation2d<IntType, RealType>& otherqmesh):
  length(otherqmesh.length),
  generators(generators),
  guardGenerators(otherqmesh.guardGenerators),
  nodes(),
  edges(),
  cellEdges() {
  xmin[0] = otherqmesh.xmin[0];
  xmin[1] = otherqmesh.xmin[1];
  xmax[0] = otherqmesh.xmax[0];
  xmax[1] = otherqmesh.xmax[1];
}

//------------------------------------------------------------------------------
// Convert real coordinates to integers.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation2d<IntType, RealType>::
quantize(const RealType* realcoords, IntType* intcoords) const {
  const RealType dx = length/(coordMax - coordMin);
  intcoords[0] = coordMin + IntType((realcoords[0] - xmin[0])/dx);
  intcoords[1] = coordMin + IntType((realcoords[1] - xmin[1])/dx);
}

//------------------------------------------------------------------------------
// Convert int coordinates to reals.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation2d<IntType, RealType>::
dequantize(const IntType* intcoords, RealType* realcoords) const {
  const RealType dx = length/(coordMax - coordMin);
  realcoords[0] = xmin[0] + (intcoords[0] - coordMin)*dx;
  realcoords[1] = xmin[1] + (intcoords[1] - coordMin)*dx;
}

//------------------------------------------------------------------------------
// Read out the current quantized tessellation state to a regular tessellation.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation2d<IntType, RealType>::
fillTessellation(Tessellation<2, RealType>& mesh) const {
  const unsigned numNodes = nodes.size();
  const unsigned numFaces = edges.size();
  const unsigned numGenerators = generators.size();
  mesh.nodes.resize(2*numNodes);
  mesh.faces.resize(numFaces, std::vector<unsigned>(2));
  mesh.faceCells.resize(numFaces);
  mesh.cells = cellEdges;
  Point2<RealType> p;
  for (unsigned i = 0; i != numNodes; ++i) {
    this->dequantize(&nodes[i].x, &mesh.nodes[2*i]);
  }
  for (unsigned i = 0; i != numFaces; ++i) {
    POLY_ASSERT(mesh.faces[i].size() == 2);
    mesh.faces[i][0] = edges[i].first;
    mesh.faces[i][1] = edges[i].second;
  }
  for (unsigned i = 0; i != numGenerators; ++i) {
    const unsigned nf = mesh.cells[i].size();
    for (unsigned j = 0; j != nf; ++j) {
      int k = mesh.cells[i][j];
      if (k < 0) {
        POLY_ASSERT2(~k < numFaces, k << " " << ~k << " " << numFaces);
        mesh.faceCells[~k].push_back(~i);
      } else {
        POLY_ASSERT2(k < numFaces, k << " " << numFaces);
        mesh.faceCells[k].push_back(i);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Internal method to construct once we have the bounds set.
// Assumes xmin, xmax, and length are already set.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation2d<IntType, RealType>::
construct(const std::vector<RealType>& points) {
  POLY_ASSERT(points.size() % 2 == 0);
  // POLY_ASSERT(std::abs(xmax[0] - xmin[0] - length) < 1e-10*length);
  // POLY_ASSERT(std::abs(xmax[1] - xmin[1] - length) < 1e-10*length);
  const int numGenerators = points.size()/2;
  generators.resize(numGenerators);
  for (unsigned i = 0; i < numGenerators; ++i) {
    this->quantize(&points[2*i], &generators[i].x);
    generators[i].index = i;
  }
  const unsigned nx = 2;
  IntType dx = (coordMax - coordMin)/nx;
  unsigned k = numGenerators;
  for (unsigned ix = 0; ix <= nx; ++ix) {
    guardGenerators.push_back(IntPoint(coordMin + ix*dx, coordMin, k++));
    guardGenerators.push_back(IntPoint(coordMin + ix*dx, coordMax, k++));
    guardGenerators.push_back(IntPoint(coordMin, coordMin + ix*dx, k++));
    guardGenerators.push_back(IntPoint(coordMax, coordMin + ix*dx, k++));
  }
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType> 
IntType
QuantizedTessellation2d<IntType, RealType>::coordMin = std::numeric_limits<IntType>::min()/2;

template<typename IntType, typename RealType> 
IntType
QuantizedTessellation2d<IntType, RealType>::coordMax = std::numeric_limits<IntType>::max()/2;

template struct QuantizedTessellation2d<int, double>;
template struct QuantizedTessellation2d<uint64_t, double>;
}
