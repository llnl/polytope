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

  Cell() = default;

  Cell(const std::vector<PointType>& points,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, facets);
  }

  Cell(const std::vector<PointType>& points,
       const std::vector<int>& faceIndices,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, faceIndices, facets);
  }

  Cell(const std::vector<CoordType>& points,
       const std::vector<int>& faceIndices,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, faceIndices, facets);
  }

  const auto& points() const noexcept {
    return m_points;
  }

  auto& points() noexcept {
    return m_points;
  }

  size_t size() const noexcept {
    return m_points.size();
  }

  bool empty() const noexcept {
    return m_points.empty();
  }

  PointType& operator[](std::size_t index) noexcept {
    return m_points[index];
  }

  const PointType& operator[](std::size_t index) const noexcept {
    return m_points[index];
  }

  // Extract using the points and a vector of vector of indices
  void init(const std::vector<PointType>& points,
            const std::vector<std::vector<unsigned>>& facets) {
    m_points.reserve(facets.size());
    for (const auto& f : facets) {
      m_points.push_back(points[f[0]]);
    }
  }

  // Extract with a layer of indirection
  void init(const std::vector<PointType>& points,
            const std::vector<int>& faceIndices,
            const std::vector<std::vector<unsigned>>& facets) {
    m_points.reserve(faceIndices.size());
    for (const auto& f : faceIndices) {
      if (f < 0) {
        m_points.push_back(points[facets[~f][1]]);
      } else {
        m_points.push_back(points[facets[f][0]]);
      }
    }
  }

  // Extract with a layer of indirection
  void init(const std::vector<CoordType>& points,
            const std::vector<int>& faceIndices,
            const std::vector<std::vector<unsigned>>& facets) {
    m_points.reserve(faceIndices.size());
    for (const auto& f : faceIndices) {
      if (f < 0) {
        int findx = facets[~f][1];
        CoordType f0 = points[2*findx];
        CoordType f1 = points[2*findx+1];
        m_points.push_back(Point2<CoordType>(f0, f1));
      } else {
        int findx = facets[f][0];
        CoordType f0 = points[2*findx];
        CoordType f1 = points[2*findx+1];
        m_points.push_back(Point2<CoordType>(f0, f1));
      }
    }
  }

  bool operator==(const Cell& other) const {
    auto lhs = m_points;
    auto rhs = other.m_points;
    std::sort(lhs.begin(), lhs.end());
    std::sort(rhs.begin(), rhs.end());
    return lhs == rhs;
  }

private:
  CellType m_points;
};

template<typename CoordType> struct Cell<3, CoordType> {
  using PointType = Point<3, CoordType>;
  using CellType = std::vector<std::vector<PointType>>;

  Cell() = default;

  Cell(const std::vector<PointType>& points,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, facets);
  }

  Cell(const std::vector<PointType>& points,
       const std::vector<int>& faceIndices,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, faceIndices, facets);
  }

  Cell(const std::vector<CoordType>& points,
       const std::vector<int>& faceIndices,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, faceIndices, facets);
  }
  // Extract using the points and a vector of vector of indices
  void init(const std::vector<PointType>& points,
            const std::vector<std::vector<unsigned>>& facets) {
    m_points.reserve(facets.size());
    for (const auto& face : facets) {
      m_points.push_back(std::vector<PointType>());
      for (const auto& f : face) {
        m_points.back().push_back(points[f]);
      }
    }
  }

  // Extract with a layer of indirection
  void init(const std::vector<PointType>& points,
            const std::vector<int>& faceIndices,
            const std::vector<std::vector<unsigned>>& facets) {
    m_points.reserve(faceIndices.size());
    for (const auto& fi : faceIndices) {
      m_points.push_back(std::vector<PointType>());
      for (const auto& f : facets[fi]) {
        if (f < 0) {
          m_points.back().push_back(points[~f]);
        } else {
          m_points.back().push_back(points[f]);
        }
      }
    }
  }

  // Extract with a layer of indirection
  void init(const std::vector<CoordType>& points,
            const std::vector<int>& faceIndices,
            const std::vector<std::vector<unsigned>>& facets) {
    m_points.reserve(faceIndices.size());
    for (const auto& fi : faceIndices) {
      m_points.push_back(std::vector<PointType>());
      for (const auto& f : facets[fi]) {
        if (f < 0) {
          int findx = facets[~f];
          CoordType f0 = points[3*findx];
          CoordType f1 = points[3*findx+1];
          CoordType f2 = points[3*findx+2];
          m_points.back().push_back(Point3<CoordType>(f0, f1, f2));
        } else {
          int findx = facets[f];
          CoordType f0 = points[3*findx];
          CoordType f1 = points[3*findx+1];
          CoordType f2 = points[3*findx+2];
          m_points.back().push_back(Point3<CoordType>(f0, f1, f2));
        }
      }
    }
  }

  const auto& points() const noexcept {
    return m_points;
  }

  auto& points() noexcept {
    return m_points;
  }

  size_t size() const noexcept {
    return m_points.size();
  }

  bool empty() const noexcept {
    return m_points.empty();
  }

  std::vector<PointType>& operator[](std::size_t index) noexcept {
    return m_points[index];
  }

  const std::vector<PointType>& operator[](std::size_t index) const noexcept {
    return m_points[index];
  }

  bool operator==(const Cell& other) const {
    auto lhs = m_points;
    auto rhs = other.m_points;
    std::sort(lhs.begin(), lhs.end());
    std::sort(rhs.begin(), rhs.end());
    return lhs == rhs;
  }

private:
  CellType m_points;
};

template<typename CoordType>
std::ostream&
operator<<(std::ostream& s, const Cell<2, CoordType>& cell) {
  s << "v = [";
  for (const auto& p : cell.points()) {
    s << p << ", ";
  }
  s << "]\n";
  return s;
}

template<typename CoordType>
std::ostream&
operator<<(std::ostream& s, const Cell<3, CoordType>& cell) {
  s << "v = [";
  for (const auto& face : cell.points()) {
    s << "[";
    for (const auto& p : face) {
      s << p << ", ";
    }
    s << "]";
  }
  s << "]\n";
  return s;
}
}

#endif
