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
#include "SiloWriter.hh"

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
    cornerIndices = addBoxPoints(node2id, result.nodes, result.quantizer());
    genVPS.resize(input.points.size());
    auto N = static_cast<long double>(input.points.size());
    // Determine the mean location of the generator points
    for (auto& p : input.points) {
      auto rp = p.template type_cast<long double>();
      centerPoint += rp/N;
    }
  }

  //-----------------------------------------------------------------------------//
  // Routines used by the specific tessellators
  //-----------------------------------------------------------------------------//
  std::set<GenPair> genPairs; // Track which gen pairs have been found
  // List of indices into vps for each generator index
  std::vector<std::vector<unsigned>> genVPS;
  std::vector<Point<2, double>> tessNodes; // Nodes provided by tessellator
  Point<2, long double> centerPoint; // Mean location of all generators

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
    std::vector<bool> invalidTris;
    invalidTris.reserve(ntri);

    // Extract the circumcenters of triangles (these become Voronoi vertices)
    for (auto i = 0u; i < ntri; ++i) {
      int ia = triangles[i][0];
      int ib = triangles[i][1];
      int ic = triangles[i][2];
      auto a = result.points[ia].template type_cast<long double>();
      auto b = result.points[ib].template type_cast<long double>();
      auto c = result.points[ic].template type_cast<long double>();
      Point2<double> out;
      invalidTris.push_back(circumcenter(a, b, c, out));
      tessNodes.push_back(out);
    }

    // Each Delaunay side has one dual Voronoi primitive: a segment for an
    // interior side and a ray for a hull side.  Emit interior sides once.
    for (auto triIndex = 0u; triIndex < ntri; ++triIndex) {
      // FIXME: Handle invalid triangles better
      //if (invalidTris[triIndex]) continue;
      const auto& tri = triangles[triIndex];
      const auto& adjacent = neighbors[triIndex];
      for (auto side = 0u; side < 3; ++side) {
        const int neighbor = adjacent[side];
        const int gen0 = tri[(side + 1) % 3];
        const int gen1 = tri[(side + 2) % 3];

        if (neighbor == -1) {// || invalidTris[neighbor]) {
          addRay(gen0, gen1, triIndex, tri[side]);
        } else {
          POLY_ASSERT(neighbor >= 0 &&
                      static_cast<size_t>(neighbor) < ntri);
          if (triIndex < static_cast<size_t>(neighbor)) {
            addFiniteEdge(gen0, gen1, triIndex, neighbor, tri[side]);
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
                     const int vertex1,
                     const int gen2 = -1) {
    addFiniteEdge(gen0, gen1, tessNodes[vertex0], tessNodes[vertex1], gen2);
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

  // Insert ray and use mean generator location to determine direction
  void addRay(const int gen0,
              const int gen1,
              const int vertex) {
    addRay(gen0, gen1, tessNodes[vertex]);
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
                     const Point<2, double>& node1,
                     const int gen2 = -1) {
    RawEdge edge;
    if (gen2 >= 0) {
      edge = makeRawEdge(gen0, gen1, gen2);
    } else {
      auto dir = node1 - node0;
      edge = makeRawEdge(gen0, gen1, rayDir(gen0, gen1, dir));
    }
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
    RawEdge edge = makeRawEdge(gen0, gen1, rayDir(gen0, gen1, dir));
    edge.hasStart = true;
    edge.start = node0;
    addVP(gen0, gen1, edge);
  }

  // Insert ray with a vertex and use the generator mean location to
  // determine the direction
  void addRay(const int gen0,
              const int gen1,
              const Point<2, double>& node0) {
    RawEdge edge = makeRawEdge(gen0, gen1, rayDirMean(gen0, gen1));
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
    RawEdge edge = makeRawEdge(gen0, gen1, gen2);
    edge.hasStart = true;
    edge.start = node0;
    addVP(gen0, gen1, edge);
  }

  // Insert ray that goes to infinity
  void addStartRay(const int gen0,
                   const int gen1,
                   const Point<2, double>& node0) {
    RawEdge edge = makeRawEdge(gen0, gen1, rayDir(gen0, gen1));
    edge.hasStart = true;
    edge.start = node0;
    addVP(gen0, gen1, edge);
  }

  // Insert ray that starts at infinity
  void addEndRay(const int gen0,
                 const int gen1,
                 const Point<2, double>& node1) {
    RawEdge edge = makeRawEdge(gen0, gen1, rayDir(gen0, gen1));
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

  // Make raw edge using 2 generators and a bisector direction
  RawEdge makeRawEdge(const int gen0,
                      const int gen1,
                      const QuantizedPoint<2>& direction) const {
    RawEdge edge;
    edge.bisectorPoint = midPoint(result.points[gen0], result.points[gen1]);
    edge.direction = direction;
    return edge;
  }

  // Make raw edge using 3 generator points
  RawEdge makeRawEdge(const int gen0,
                      const int gen1,
                      const int gen2) const {
    RawEdge edge;
    const auto& a = result.points[gen0];
    const auto& b = result.points[gen1];
    const auto& c = result.points[gen2];
    edge.bisectorPoint = midPoint(a, b);
    const auto winding = qcross(b - a, c - a);
    POLY_CHECK(winding != 0);
    if (winding > 0) {
      edge.leftGen = gen1;
      edge.direction = outwardRay(b, a);
    } else if (winding < 0) {
      edge.leftGen = gen0;
      edge.direction = outwardRay(a, b);
    }
    return edge;
  }

  // Determine ray direction using the mean point
  QuantizedPoint<2> rayDirMean(const int gen0,
                               const int gen1) const {
    const auto& qa = result.points[gen0];
    const auto& qb = result.points[gen1];

    const auto a = qa.template type_cast<long double>();
    const auto b = qb.template type_cast<long double>();
    const auto midpoint = (a + b) * 0.5;

    auto direction = outwardRay(qa, qb);
    const auto directionD = direction.template type_cast<long double>();
    const long double side = dot(directionD, midpoint - centerPoint);
    POLY_CHECK2(std::abs(side) > 1.E-12,
                "Cannot determine hull-ray direction from center");
    if (side < 0.0) {
      direction = -direction;
    }
    return direction;
  }

  // Determine the get outward direction
  QuantizedPoint<2> rayDir(const int gen0,
                           const int gen1) {
    return outwardRay(result.points[gen0], result.points[gen1]);
  }

  // Determine the ray direction using floating point direction
  QuantizedPoint<2> rayDir(const int gen0,
                           const int gen1,
                           const Point<2, double>& dir) {
    auto ray = outwardRay(result.points[gen0], result.points[gen1]);
    if (dot(dir, ray.template type_cast<double>()) < 0.) {
      ray = -ray;
    }
    return ray;
  }

  // Clip once, then store a canonical edge directed with gp.first on its left.
  void addVP(const int gen0,
             const int gen1,
             const RawEdge& rawEdge) {
    ClippedEdge vp(gen0, gen1);
    if (genPairs.count(vp.gp) != 0) return;

    Clip2D<QuantizedCoordinate<2>> clipper(rawEdge, result.quantizer());
    if (!clipper.clip()) return;

    vp.csides = std::make_pair(clipper.startSide, clipper.endSide);
    vp.edge = updateNodeMap(clipper.p0, clipper.p1, node2id, result.nodes);
    if (rawEdge.leftGen >= 0) {
      if (rawEdge.leftGen != vp.gp.first) {
        vp.flipEdge();
      }
    } else {
      orientClippedEdge(vp, clipper.p0, clipper.p1, result.points[vp.gp.first]);
    }

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
    std::string errorStr =
      makeCCWCellFromClippedEdges<QuantizedCoordinate<2>>(clippedEdges,
                                                          cornerIndices,
                                                          result.faces,
                                                          edgeToFace,
                                                          result.cells[cellIndex]);
    // If error occurs above, generate an escape pod
    if (!errorStr.empty() && !result.escapePod) {
      // FIXME: Move functionality elsewhere
      using PodTess = Tessellation<2, QuantizedCoordinate<2>>;
      PodTess escapePod;
      escapePod.escapePod = true;
      escapePod.tessellatorName = result.tessellatorName;
      escapePod.points.push_back(result.points[cellIndex]);
      std::vector<int> cell;
      int curCell = 0;
      unsigned faceIndex = 0;
      // Loop over current clipped edges and fill the escape pod tessellation
      std::set<int> foundGens;
      foundGens.insert(cellIndex);
      for (const auto i : clippedEdges) {
        int gen = (i.gp.first == cellIndex) ? i.gp.second : i.gp.first;
        if (foundGens.count(gen) == 0) {
          foundGens.insert(gen);
          escapePod.points.push_back(result.points[gen]);
        }
        // Loop over any generators paired with paired generators
        for (const auto ig : genVPS[gen]) {
          auto& gpi = vps[ig].gp;
          int geni0 = gpi.first;
          int geni1 = gpi.second;
          if (foundGens.count(geni0) == 0) {
            foundGens.insert(geni0);
            escapePod.points.push_back(result.points[geni0]);
          }
          if (foundGens.count(geni1) == 0) {
            foundGens.insert(geni1);
            escapePod.points.push_back(result.points[geni1]);
          }
        }
        int e0 = i.edge.first;
        int e1 = i.edge.second;
        if (e0 >= 0) {
          escapePod.nodes.push_back(result.nodes[e0]);
          faceIndex++;
        }
        if (e1 >= 0) {
          escapePod.nodes.push_back(result.nodes[e1]);
          faceIndex++;
        }
        if (e0 >= 0 && e1 >= 0) {
          escapePod.faces.push_back({faceIndex-2, faceIndex-1});
          cell.push_back(curCell++);
        }
      }
      escapePod.cells.push_back(cell);
      SiloWriter<2, PodTess> writer(escapePod);
      writer.setOvlType(false);
      std::string filename = "narrowEscapePod" + std::to_string(Communicator::getRank());
      writer.write(filename, "", 1);
      // Also write the entire tessellation but remove unfilled cells
      result.cells.resize(cellIndex);
      result.escapePod = true;
      SiloWriter<2, QuantTessellation<2>> qwriter(result);
      qwriter.setOvlType(false);
      filename = "wideEscapePod" + std::to_string(Communicator::getRank());
      qwriter.write(filename, "", 1);
      POLY_CHECK2(false, errorStr << " wrote escape pod file " << filename);
    } else if (!errorStr.empty()) {
      POLY_CHECK2(false, errorStr);
    }
  }
};

} // namespace polytope

#endif
