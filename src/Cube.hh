#ifndef __Polytope_Cube__
#define __Polytope_Cube__

//------------------------------------------------------------------------------
// Simple cube routines  3D.
//------------------------------------------------------------------------------

#include "Point.hh"

namespace polytope {

template<typename CoordType>
class Cube {
public:
  using PointType = Point<3, CoordType>;
  std::vector<PointType> nodes;
  std::vector<std::vector<unsigned>> faces = {{0, 1, 2, 3},  // bottom (-z)
                                              {4, 7, 6, 5},  // top (+z)
                                              {0, 4, 5, 1},  // front (-y)
                                              {2, 6, 7, 3},  // back (+y)
                                              {0, 3, 7, 4},  // left (-x)
                                              {1, 5, 6, 2}   // right (+x);
  };

  // Default constructor
  Cube() = default;

  // 3D version using bit pattern (no specific ordering required yet)
  Cube(const Point3<CoordType>& min,
       const Point3<CoordType>& max) {
    init(min, max);
  }

  inline std::vector<std::vector<unsigned>> createCubeFaces() {
    return faces;
  }

  inline std::vector<CoordType> flatNodes() {
    return flattenCoords(nodes);
  }

  bool within(const PointType& pos) {
    return (pos.x >= nodes[0].x && pos.x <= nodes[1].x) &&
      (pos.y >= nodes[0].y && pos.y <= nodes[2].y) &&
      (pos.z >= nodes[0].z && pos.z <= nodes[4].z);
  }

    // 3D version using bit pattern (no specific ordering required yet)
  void init(const Point3<CoordType>& min,
            const Point3<CoordType>& max) {
    int count = 8;
    unsigned kk = 0;
    for (int i = 0; i < count; ++i) {
      Point3<CoordType> corner;
      for (int d = 0; d < 3; ++d) {
        corner[d] = min[d] + max[d]*((i >> d) & 1);
      }
      corner.index = kk++;
      nodes.push_back(corner);
    }
  }
};
}
#endif
