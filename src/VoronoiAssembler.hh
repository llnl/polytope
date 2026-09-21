//-----------------------------------------------------------------------------//
// VoronoiAssembler
//
// Convert backend-produced Voronoi primitives into a bounded quantized
// tessellation.  Specializations define the primitive representation and the
// dimension-specific clipping/topology construction.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_VoronoiAssembler__
#define __Polytope_VoronoiAssembler__

#include "QuantTessellation.hh"
#include "EdgeUtils.hh"
#include "Intersections.hh"
#include "Shapes.hh"
#include "Clipping2D.hh"

#include <vector>

namespace polytope {

template<int Dimension>
struct VoronoiPrimitive {};

template<int Dimension>
using VoronoiPrimitiveCells =
  std::vector<std::vector<VoronoiPrimitive<Dimension>>>;

template<int Dimension>
class VoronoiAssembler;

template<>
struct VoronoiPrimitive<2> {
  GenPair gp;
  Point2<double> rp0;
  Point2<double> rp1;
  int gp3 = -1;

  bool inf0 = true;
  bool inf1 = true;

  VoronoiPrimitive(const int gp0,
                   const int gp1):
    gp(orderGenPair(gp0, gp1)) { }
  void setV0(const Point2<double>& v0) {
    inf0 = false;
    rp0 = v0;
  }
  void setV1(const Point2<double>& v1) {
    inf1 = false;
    rp1 = v1;
  }
  void setThirdPoint(const int thirdIndex) {
    gp3 = thirdIndex;
  }
};

template<>
class VoronoiAssembler<2> {
public:
  VoronoiAssembler(QuantTessellation<2>& input):
    result(input) {
    cornerIndices = addBoxPoints(node2id, result.nodes);
  }

  //! Assemble all cells using one shared edge/node cache.
  void assemble(const VoronoiPrimitiveCells<2>& cellPrimitives) {
    POLY_ASSERT(cellPrimitives.size() == result.points.size());
    result.cells.resize(cellPrimitives.size());
    for (auto cellIndex = 0u; cellIndex < cellPrimitives.size(); ++cellIndex) {
      constructEdges(cellPrimitives[cellIndex], cellIndex);
    }
  }

private:
  std::map<QuantizedPoint<2>, int> node2id;
  EdgeToFaceMap edgeToFace;
  GenPairToClippedEdgeMap genPairToEdge;
  std::map<BoxSide, unsigned> cornerIndices;
  QuantTessellation<2>& result;

  void constructEdges(const std::vector<VoronoiPrimitive<2>>& vps,
                      const int cellIndex) {
    std::vector<ClippedEdge> clippedEdges;
    for (const auto& vp : vps) {
      const auto& gp = vp.gp;
      auto cacheIt = genPairToEdge.find(gp);
      if (cacheIt != genPairToEdge.end()) {
        clippedEdges.push_back(flipEdge(cacheIt->second));
        continue;
      }
      ClippedEdge clippedEdge;
      Clip2D<QuantizedCoordinate<2>> clipper;
      const int gindx1 = cellIndex;
      const int gindx2 = (gp.first != cellIndex) ? gp.first : gp.second;
      clipper.gen0 = result.points[gindx1];
      clipper.gen1 = result.points[gindx2];
      if (vp.inf0) {
        clipper.inf0 = true;
        clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1);
      } else {
        clipper.rp0 = vp.rp0;
      }
      if (vp.inf1) {
        clipper.inf1 = true;
        clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1);
      } else {
        clipper.rp1 = vp.rp1;
      }
      if (vp.gp3 >= 0) {
        clipper.normalRay = outwardRay(clipper.gen0, clipper.gen1, result.points[vp.gp3]);
      }
      if (!vp.inf0 && !vp.inf1) {
        clipper.normalRay = pointDirection<QuantizedCoordinate<2>>(clipper.rp0, clipper.rp1);
      }
      if (clipper.doClipping()) {
        continue;
      }
      clippedEdge.clippedSides = std::make_pair(clipper.ifirstSide, clipper.isecondSide);
      clippedEdge.curEdge = updateNodeMap(clipper.p0, clipper.p1, node2id, result.nodes);
      genPairToEdge[gp] = clippedEdge;
      clippedEdges.push_back(clippedEdge);
    }
    result.cells[cellIndex] =
      makeCCWCellFromClippedEdges<QuantizedCoordinate<2>>(clippedEdges, result.points[cellIndex],
                                                          cornerIndices, result.nodes, result.faces, edgeToFace);
  }
};

} // namespace polytope

#endif
