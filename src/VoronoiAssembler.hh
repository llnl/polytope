//-----------------------------------------------------------------------------//
// VoronoiAssembler
//
// Convert backend-produced ClippedEdges into a bounded quantized tessellation.
// Specializations define the dimension-specific clipping/topology construction.
//-----------------------------------------------------------------------------//
#ifndef __Polytope_VoronoiAssembler__
#define __Polytope_VoronoiAssembler__

#include "QuantTessellation.hh"
#include "EdgeUtils.hh"
#include "Intersections.hh"
#include "Shapes.hh"
#include "Clipping2D.hh"

#include <array>
#include <vector>

namespace polytope {

template<int Dimension>
class VoronoiAssembler;

template<>
class VoronoiAssembler<2> {
public:
  using RawEdge = DirectedVoronoiEdge<QuantizedCoordinate<2>>;
  std::vector<ClippedEdge> vps;
  VoronoiAssembler(QuantTessellation<2>& input):
    result(input) {
    cornerIndices = addBoxPoints(node2id, result.nodes);
    genVPS.resize(input.points.size());
  }

  //-----------------------------------------------------------------------------//
  // Routines used by the specific tessellators
  //-----------------------------------------------------------------------------//
  std::set<GenPair> genPairs; // Track which gen pairs have been found
  // List of indices into vps for each generator index
  std::vector<std::vector<unsigned>> genVPS;
  std::vector<Point<2, double>> tessNodes; // Nodes provided by tessellator

  //-----------------------------------------------------------------------------//
  // Handle special cases like 2 generators or collinear generators
  //-----------------------------------------------------------------------------//
  void degenerateAssembly() {
    const auto N = result.points.size();
    // Implies collinear or 2 generator case
    for (auto cellIndex = 0u; cellIndex < N-1; ++cellIndex) {
      addInfLines(cellIndex, cellIndex+1);
    }
  }

  //-----------------------------------------------------------------------------//
  // Assemble Voronoi primitives from a Delaunay triangulation. For every
  // triangle i, triangles[i][v] is its generator at local vertex v.
  // neighbors[i][side] is the adjacent triangle across the side opposite
  // triangles[i][side], or -1 when that side belongs to the convex hull.
  // Triangle indices are also the indices of the generated circumcenters.
  //-----------------------------------------------------------------------------//
  void assembleDelaunay(const std::vector<std::array<int, 3>>& triangles,
                        const std::vector<std::array<int, 3>>& neighbors) {
    POLY_ASSERT(triangles.size() == neighbors.size());
    const auto ntri = triangles.size();
    tessNodes.reserve(ntri);

    // Extract the circumcenters of triangles (these become Voronoi vertices)
    for (auto i = 0u; i < ntri; ++i) {
      int ia = triangles[i][0];
      int ib = triangles[i][1];
      int ic = triangles[i][2];
      auto a = result.points[ia].template type_cast<double>();
      auto b = result.points[ib].template type_cast<double>();
      auto c = result.points[ic].template type_cast<double>();
      tessNodes.push_back(circumcenter(a, b, c));
    }

    // Each Delaunay side has one dual Voronoi primitive: a segment for an
    // interior side and a ray for a hull side.  Emit interior sides once.
    for (auto triIndex = 0u; triIndex < ntri; ++triIndex) {
      const auto& tri = triangles[triIndex];
      const auto& adjacent = neighbors[triIndex];
      for (auto side = 0u; side < 3; ++side) {
        const int neighbor = adjacent[side];
        const int gen0 = tri[(side + 1) % 3];
        const int gen1 = tri[(side + 2) % 3];

        if (neighbor == -1) {
          addRay(gen0, gen1, triIndex, tri[side]);
        } else {
          POLY_ASSERT(neighbor >= 0 &&
                      static_cast<size_t>(neighbor) < ntri);
          if (triIndex < static_cast<size_t>(neighbor)) {
            addFiniteEdge(gen0, gen1, triIndex, neighbor);
          }
        }
      }
    }
  }

  // Insert edge that is infinite in both directions
  void addInfLines(const int gen0,
                   const int gen1) {
    RawEdge edge = makeRawEdge(gen0, gen1,
                               outwardRay(result.points[gen0], result.points[gen1]));
    addVP(gen0, gen1, edge);
  }

  //-----------------------------------------------------------------------------//
  // These routines are for tessellators that can provide the nodes list
  // explicitly and can set edges based on the indices to that list
  //-----------------------------------------------------------------------------//
  void fillTessNodes(const std::vector<Point<2, double>>& nodes) {
    tessNodes = nodes;
  }

  // Add finite edge segment
  void addFiniteEdge(const int gen0,
                     const int gen1,
                     const int vertex0,
                     const int vertex1) {
    addFiniteEdge(gen0, gen1, tessNodes[vertex0], tessNodes[vertex1]);
  }

  // Insert ray with a vertex and direction
  void addRay(const int gen0,
              const int gen1,
              const int vertex0,
              const Point<2, double>& dir) {
    addRay(gen0, gen1, tessNodes[vertex0], dir);
  }

  // Insert ray with third point for determining ray direction
  void addRay(const int gen0,
              const int gen1,
              const int vertex,
              const int gen2) {
    addRay(gen0, gen1, tessNodes[vertex], gen2);
  }

  // Insert ray that goes to infinity
  void addStartRay(const int gen0,
                   const int gen1,
                   const int vertex0) {
    addStartRay(gen0, gen1, tessNodes[vertex0]);
  }

  // Insert ray that starts at infinity
  void addEndRay(const int gen0,
                 const int gen1,
                 const int vertex1) {
    addEndRay(gen0, gen1, tessNodes[vertex1]);
  }

  //-----------------------------------------------------------------------------//
  // These routines are for tessellators that provide edge endpoints directly.
  //-----------------------------------------------------------------------------//
  void addFiniteEdge(const int gen0,
                     const int gen1,
                     const Point<2, double>& node0,
                     const Point<2, double>& node1) {
    RawEdge edge = makeRawEdge(gen0, gen1,
                               pointDirection<QuantizedCoordinate<2>>(node0, node1));
    edge.hasStart = true;
    edge.start = node0;
    edge.hasEnd = true;
    edge.end = node1;
    addVP(gen0, gen1, edge);
  }

  // Insert ray with a vertex and a direction
  void addRay(const int gen0,
              const int gen1,
              const Point<2, double>& node0,
              const Point<2, double>& dir) {
    // For consistency, recompute the direction using the generators
    // Only use the provided direction to determine generator order
    auto ray = rayDirection(gen0, gen1);
    if (dot(dir, ray.template type_cast<double>()) < 0.) {
      ray = -ray;
    }
    RawEdge edge = makeRawEdge(gen0, gen1, ray);
    edge.hasStart = true;
    edge.start = node0;
    addVP(gen0, gen1, edge);
  }

  // Insert edge that goes to or starts at infinity with third point
  // for determining ray direction
  void addRay(const int gen0,
              const int gen1,
              const Point<2, double>& node0,
              const int gen2) {
    RawEdge edge = makeRawEdge(gen0, gen1,
                               rayDirection(gen0, gen1, gen2));
    edge.hasStart = true;
    edge.start = node0;
    addVP(gen0, gen1, edge);
  }

  // Insert ray that goes to infinity
  void addStartRay(const int gen0,
                   const int gen1,
                   const Point<2, double>& node0) {
    RawEdge edge = makeRawEdge(gen0, gen1, rayDirection(gen0, gen1));
    edge.hasStart = true;
    edge.start = node0;
    addVP(gen0, gen1, edge);
  }

  // Insert ray that starts at infinity
  void addEndRay(const int gen0,
                 const int gen1,
                 const Point<2, double>& node1) {
    RawEdge edge = makeRawEdge(gen0, gen1, rayDirection(gen0, gen1));
    edge.hasEnd = true;
    edge.end = node1;
    addVP(gen0, gen1, edge);
  }

  //-----------------------------------------------------------------------------//
  // Routines used by the SerialTessellator
  // Assemble all cells using one shared edge/node cache.
  //-----------------------------------------------------------------------------//
  void assemble() {
    const auto N = result.points.size();
    result.cells.resize(N);
    for (auto cellIndex = 0u; cellIndex < N; ++cellIndex) {
      constructEdges(cellIndex);
    }
  }

private:
  std::map<QuantizedPoint<2>, int> node2id;
  EdgeToFaceMap edgeToFace;
  std::map<BoxSide, unsigned> cornerIndices;
  QuantTessellation<2>& result;

  RawEdge makeRawEdge(const int gen0,
                      const int gen1,
                      const QuantizedPoint<2>& direction) const {
    RawEdge edge;
    edge.bisectorPoint = midPoint(result.points[gen0], result.points[gen1]);
    edge.direction = direction;
    return edge;
  }

  QuantizedPoint<2> rayDirection(const int gen0,
                                 const int gen1,
                                 const int gen2 = -1) const {
    if (gen2 >= 0) {
      return outwardRay(result.points[gen0], result.points[gen1],
                        result.points[gen2]);
    }
    return outwardRay(result.points[gen0], result.points[gen1]);
  }

  // Clip once, then store a canonical edge directed with gp.first on its left.
  void addVP(const int gen0,
             const int gen1,
             const RawEdge& rawEdge) {
    ClippedEdge vp(gen0, gen1);
    if (genPairs.count(vp.gp) != 0) return;

    Clip2D<QuantizedCoordinate<2>> clipper(rawEdge);
    if (!clipper.clip()) return;

    vp.csides = std::make_pair(clipper.startSide, clipper.endSide);
    vp.edge = updateNodeMap(clipper.p0, clipper.p1, node2id, result.nodes);
    orientClippedEdge(vp, clipper.p0, clipper.p1, result.points[vp.gp.first]);

    const auto edgeIndex = vps.size();
    genPairs.insert(vp.gp);
    genVPS[vp.gp.first].push_back(edgeIndex);
    genVPS[vp.gp.second].push_back(edgeIndex);
    vps.push_back(vp);
  }

  void constructEdges(const int cellIndex) {
    // Make local clipped edges for this particular generator point
    std::vector<ClippedEdge> clippedEdges;
    // Retrieve all clipped edges for this generator point
    const std::vector<unsigned> curVPS = genVPS[cellIndex];
    for (const auto i : curVPS) {
      const auto& vp = vps[i];
      const auto& gp = vp.gp;
      // Edges are oriented CCW relative to the first generator in the pair
      // Flip the edges if this isn't the first generator
      bool flipEdges = (gp.first == cellIndex) ? false : true;
      if (flipEdges) {
        clippedEdges.push_back(flipEdge(vp));
      } else {
        clippedEdges.push_back(vp);
      }
    }
    result.cells[cellIndex] =
      makeCCWCellFromClippedEdges<QuantizedCoordinate<2>>(clippedEdges, cornerIndices, result.faces, edgeToFace);
  }
};

} // namespace polytope

#endif
