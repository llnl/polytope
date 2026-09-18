#ifndef __Polytope_VoronoiConstructor__
#define __Polytope_VoronoiConstructor__

#include "QuantTessellation.hh"
#include "EdgeUtils.hh"
#include "Intersections.hh"
#include "Shapes.hh"
#include "Clipping2D.hh"

namespace polytope {

struct VoronoiPrimitive {
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

class VoronoiConstructor {
public:
  // Map QuantizedPoint coordinates to our node indices (for deduplication)
  std::map<QuantizedPoint<2>, int> node2id;
  // Map canonical edges to face indices for oriented edge tracking
  EdgeToFaceMap edgeToFace;
  // Map generator pairs to ClippedEdges
  // ClippedEdge contains an edge and two ints for the box sides that
  // possibly clipped the start and the end of the segment
  GenPairToClippedEdgeMap genPairToEdge;
  // Corner indices
  std::map<BoxSide, unsigned> cornerIndices;
  // Reference to the quantized tessellation
  QuantTessellation<2>& result;

  VoronoiConstructor(QuantTessellation<2>& input) :
    result(input) {
    cornerIndices = addBoxPoints(node2id, result.nodes);
  }

  // For a set of VoronoiPrimitives, create an unordered vector of ClippedEdges
  // and update the nodes, faces, and cells in the result QuantTessellation.
  void constructEdges(std::vector<VoronoiPrimitive>& vps,
                      const int cellIndex) {
    std::vector<ClippedEdge> clippedEdges;
    for (auto& vp : vps) {
      auto& gp = vp.gp;
      auto cacheIt = genPairToEdge.find(gp);
      // Check if we have already solved for this Voronoi edge
      if (cacheIt != genPairToEdge.end()) {
        clippedEdges.push_back(flipEdge(cacheIt->second));
        continue;
      }
      ClippedEdge clippedEdge;
      // Check if edge must be clipped by quantized space
      // TODO: Simplify the clipping logic
      Clip2D<QuantizedCoordinate<2>> clipper;
      int gindx1 = cellIndex;
      int gindx2 = (gp.first != cellIndex) ? gp.first : gp.second;
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
}

#endif
