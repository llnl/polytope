
#include "TetgenTessellator.hh"
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

#define TETLIBRARY
#include "tetgen.h"

namespace polytope {

void TetgenTessellator::tessellateQuantizedImpl(QT& result) const {
  tetgenio in = createTetgenPoints(result);
  tetgenio out;
  // Create the Delaunay
  tetrahedralize((char*)"qQd", &in, &out);

  for (int i = 0; i < out.numberoftetrahedra; ++i) {
    const int* tet =
      &out.tetrahedronlist[i * out.numberofcorners];
    std::array<std::array<double, 3>, 4> tetra;
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 3; ++k) {
        tetra[j][k] = out.pointlist[3*tet[j]+k];
      }
    }
  }
}

// Create Tetgen class
tetgenio TetgenTessellator::createTetgenPoints(const QT& result) const {
  tetgenio in;
  auto generators = flattenCoords(result.getRealPoints());
  const auto N = result.m_points.size();
  in.pointlist = new REAL[N*3];
  in.numberofpoints = N;
  std::copy(generators.begin(), generators.end(), in.pointlist);
  return in;
}

}
