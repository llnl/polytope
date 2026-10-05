//-----------------------------------------------------------------------------//
// Cell
//
// Generalized class for handling cells in 2D and 3D.
// FIXME: Consolidate repeated routines
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

  Cell(const CellType& points) :
    m_points(points) {
  }

  Cell(const std::vector<PointType>& points,
       const std::vector<std::vector<unsigned>>& facets) {
    init(points, facets);
  }

  Cell(const std::vector<PointType>& points,
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

  bool operator==(const Cell& other) const {
    auto lhs = m_points;
    auto rhs = other.m_points;
    std::sort(lhs.begin(), lhs.end());
    std::sort(rhs.begin(), rhs.end());
    return lhs == rhs;
  }

  template<typename OtherType>
  Cell<2, OtherType> type_cast() const {
    std::vector<Point<2, OtherType>> rpoints;
    rpoints.reserve(size());
    for (const auto& point : m_points) {
      rpoints.push_back(point.template type_cast<OtherType>());
    }
    return Cell<2, OtherType>(rpoints);
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

template<int Dimension, typename CoordType>
bool isApprox(const std::vector<Point<Dimension, CoordType>>& points0,
              const std::vector<Point<Dimension, CoordType>>& points1,
              const CoordType& relTol) {
  for (const auto& p0 : points0) {
    bool found = false;
    for (const auto& p1 : points1) {
      if (isApprox(p0, p1, relTol)) {
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

template<typename CoordType>
bool isApprox(const Cell<2, CoordType>& cell0,
              const Cell<2, CoordType>& cell1,
              const CoordType& relTol) {
  return isApprox<2, CoordType>(cell0.points(), cell1.points(), relTol);
}

template<typename CoordType>
bool isApprox(const Cell<3, CoordType>& cell0,
              const Cell<3, CoordType>& cell1,
              const CoordType& relTol) {
  for (const auto& f0 : cell0.points()) {
    bool found = false;
    for (const auto& f1 : cell1.points()) {
      if (isApprox(f0, f1, relTol)) {
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

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
