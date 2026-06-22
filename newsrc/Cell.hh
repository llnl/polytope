//-----------------------------------------------------------------------------//
// Cell
//
// Generalized class for handling cells in 2D and 3D.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_Cell__
#define __Polytope_Cell__

#include "Point.hh"

namespace polytope {

template<int Dimension, typename CoordType> struct Cell;

template<typename CoordType> struct Cell<2, CoordType> {
  using PointType = Point<2, CoordType>;
  using CellType = std::vector<PointType>;
  // Extract using the points and a vector of vector of indices
  static CellType extractCell(const std::vector<PointType>& points,
                              const std::vector<std::vector<int>>& facets) {
    CellType facePoints;
    facePoints.reserve(facets.size());
    for (const auto& f : facets) {
      facePoints.push_back(points[f[0]]);
    }
    return facePoints;
  }
  // Extract with a layer of indirection
  static CellType extractCell(const std::vector<PointType>& points,
                              const std::vector<int>& faceIndices,
                              const std::vector<std::vector<int>>& facets) {
    CellType facePoints;
    facePoints.reserve(faceIndices.size());
    for (const auto& f : faceIndices) {
      if (f < 0) {
        facePoints.push_back(points[facets[~f][1]]);
      } else {
        facePoints.push_back(points[facets[f][0]]);
      }
    }
    return facePoints;
  }
};

template<typename CoordType> struct Cell<3, CoordType> {
  using PointType = Point<3, CoordType>;
  using CellType = std::vector<std::vector<PointType>>;
  // Extract using the points and a vector of vector of indices
  static CellType extractCell(const std::vector<PointType>& points,
                              const std::vector<std::vector<int>>& facets) {
    CellType facePoints;
    facePoints.reserve(facets.size());
    for (const auto& face : facets) {
      facePoints.push_back(std::vector<PointType>());
      for (const auto& f : face) {
        facePoints.back().push_back(points[f]);
      }
    }
    return facePoints;
  }
  // Extract with a layer of indirection
  static CellType extractCell(const std::vector<PointType>& points,
                              const std::vector<int>& faceIndices,
                              const std::vector<std::vector<int>>& facets) {
    CellType facePoints;
    facePoints.reserve(faceIndices.size());
    for (const auto& fi : faceIndices) {
      facePoints.push_back(std::vector<PointType>());
      for (const auto& f : facets[fi]) {
        if (f < 0) {
          facePoints.back().push_back(points[~f]);
        } else {
          facePoints.back().push_back(points[f]);
        }
      }
    }
    return facePoints;
  }
};
}

#endif
