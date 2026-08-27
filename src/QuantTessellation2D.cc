//-----------------------------------------------------------------------------//
// QuantTessellation<2> specializations
//-----------------------------------------------------------------------------//

#include "QuantTessellation.hh"
#include "EscapePod.hh"
#include "Intersections.hh"
#include "Tessellator.hh"
#include "EdgeUtils.hh"
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <algorithm>
#ifdef POLYTOPE_ENABLE_BOOST
#include "BoostPolygonIntersections.hh"
#include "boost/polygon/polygon_set_data.hpp"
#endif

namespace polytope {

//------------------------------------------------------------------------------
// Escape pod methods
//------------------------------------------------------------------------------
template<>
void
QuantTessellation<2>::ejectEscapePod(std::string filename,
                                     const std::vector<unsigned>& genPoints,
                                     const QuantPLC<2>& QPLC,
                                     const std::string& tessellatorName) {
  if (m_isEscapePod) return;
  const auto podName = escape_pod::filename(filename);
  std::ofstream out(podName.c_str());
  POLY_VERIFY2(out, "Unable to open escape pod file for writing: " << podName);

  out << "PolytopeQuantTessellation2DEscapePod 2\n";
  out << "rank " << Communicator::getRank() << "\n";
  out << "tessellator " << (tessellatorName.empty() ? "unknown" : tessellatorName) << "\n";

  out << "generators " << genPoints.size() << "\n";
  for (const auto genIndex : genPoints) {
    POLY_VERIFY2(genIndex < points.size(),
                 "Escape pod generator index " << genIndex
                 << " exceeds generator count " << points.size());
    out << points[genIndex] << "\n";
  }

  const bool haveQPLC = not QPLC.empty();
  out << "qplc " << haveQPLC << "\n";
  if (haveQPLC) {
    out << "qplc_points " << QPLC.points.size() << "\n";
    for (const auto& point : QPLC.points) out << point << "\n";

    out << "qplc_topology\n";
    QPLC.write(out);
  }

  POLY_VERIFY2(out, "Failed while writing escape pod file: " << podName);
}

template<>
void
QuantTessellation<2>::loadEscapePod(std::string filename,
                                    QuantPLC<2>& QPLC,
                                    std::string& tessellatorName) {
  auto podName = filename;
  std::ifstream in(podName.c_str());
  if (!in) {
    podName = escape_pod::filename(filename);
    in.clear();
    in.open(podName.c_str());
  }
  POLY_VERIFY2(in, "Unable to open escape pod file for reading: " << podName);

  escape_pod::expectToken(in, "PolytopeQuantTessellation2DEscapePod");
  unsigned version = 0;
  in >> version;
  POLY_VERIFY2(in and (version == 1 or version == 2),
               "Unsupported escape pod version " << version
               << " in file " << podName);

  escape_pod::expectToken(in, "rank");
  int fileRank = -1;
  in >> fileRank;
  POLY_VERIFY2(in, "Malformed escape pod file while reading rank");

  tessellatorName.clear();
  if (version >= 2) {
    escape_pod::expectToken(in, "tessellator");
    in >> tessellatorName;
    POLY_VERIFY2(in, "Malformed escape pod file while reading tessellator name");
  }

  clear();

  escape_pod::expectToken(in, "generators");
  unsigned ngenerators = 0;
  in >> ngenerators;
  POLY_VERIFY2(in, "Malformed escape pod file while reading generator count");
  points.resize(ngenerators);
  for (auto& point : points) {
    in >> point;
    POLY_VERIFY2(in, "Malformed escape pod file while reading a generator point");
  }
  escape_pod::rebuildTessellationPointMetadata(*this);

  escape_pod::expectToken(in, "qplc");
  bool haveQPLC = false;
  in >> haveQPLC;
  POLY_VERIFY2(in, "Malformed escape pod file while reading QPLC presence");

  QuantPLC<2> localQPLC;
  if (haveQPLC) {
    escape_pod::expectToken(in, "qplc_points");
    unsigned nqplcPoints = 0;
    in >> nqplcPoints;
    POLY_VERIFY2(in, "Malformed escape pod file while reading QPLC point count");
    localQPLC.points.resize(nqplcPoints);
    for (auto& point : localQPLC.points) {
      in >> point;
      POLY_VERIFY2(in, "Malformed escape pod file while reading a QPLC point");
    }

    escape_pod::expectToken(in, "qplc_topology");
    localQPLC.read(in);
    POLY_VERIFY2(in, "Malformed escape pod file while reading QPLC topology");
    escape_pod::rebuildQPLCPointMetadata(localQPLC);
  }
  QPLC = std::move(localQPLC);

  std::string trailing;
  POLY_VERIFY2(!(in >> trailing),
               "Unexpected trailing token in escape pod file "
               << podName << ": " << trailing);
}

template<>
void
QuantTessellation<2>::loadEscapePod(std::string filename,
                                    QuantPLC<2>& QPLC) {
  std::string tessellatorName;
  loadEscapePod(filename, QPLC, tessellatorName);
}

template<>
void
QuantTessellation<2>::cullExternalPoints(const QuantPLC<2>& QPLC) {
  const auto& Q = Quantizer<2>::instance();
  auto N = points.size();
  std::vector<QuantizedPoint<2>> newPoints;
  newPoints.reserve(N);
  std::vector<QuantizedKey<2>> newHashes;
  newHashes.reserve(N);
  unsigned indx = 0;
  for (auto i = 0u; i < N; ++i) {
    auto point = points[i];
    if (!QPLC.within(point)) continue;
    point.index = indx++;
    newPoints.push_back(point);
    newHashes.push_back(Q.encode(point));
  }
  points = std::move(newPoints);
  hashes = std::move(newHashes);
}

// All clipping functionality relies on Boost
#if defined(POLYTOPE_ENABLE_BOOST) && !defined(POLYTOPE_ENABLE_HIBIT2D)
namespace bp = boost::polygon;
// Need this to use the -=, +=, etc operators
using namespace boost::polygon::operators;

// Intersect the current cell with a boundary and extend the orphans vector if necessary.
// Return the polygon that contains the generator point
PolygonWithHoles
cellPolygonIntersect(const QuantTessellation<2>::QuantizedCell& currentCell,
                     const QuantizedPoint<2>& genPointP,
                     const PolygonWithHoles boundary,
                     std::vector<PolygonWithHoles>& orphans) {
  bp::point_data<QuantizedCoordinate<2>> genPoint =
    bp::construct<QuantizedPoint<2>>(genPointP.x, genPointP.y);
  std::vector<PolygonWithHoles> cellSet = boostIntersect(currentCell, boundary);
  unsigned fragIndex = 0;
  auto NFrag = cellSet.size();
  if (NFrag == 0) return PolygonWithHoles(); // Cell was completely outside boundary
  if (NFrag > 1) {
    // Find which part owns the generator
    while (fragIndex < NFrag and
           not bp::contains(cellSet[fragIndex], genPoint)) ++fragIndex;
    for (unsigned iPoly = 0; iPoly < NFrag; ++iPoly) {
      if (iPoly != fragIndex) {
        // Check if orphan can be added to other orphans
        // We should check all existing orphans for robustness
        bool foundUnion = false;
        for (auto& orphan : orphans) {
          if (validUnion(cellSet[iPoly], orphan)) {
            foundUnion = true;
          }
        }
        if (!foundUnion) {
          orphans.push_back(cellSet[iPoly]);
        }
      }
    }
  } // End of multiple fragment check
  return cellSet[fragIndex];
}

template<>
void
QuantTessellation<2>::clipTessellation(const QuantPLC<2>& QPLC,
                                       Tessellator<2, double>& tessellator) {
  const auto& Q = Quantizer<2>::instance();
  auto boundaryPoints = QPLC.getCell().points();
  PolygonWithHoles boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());
  auto holePoints = QPLC.getHolePoints();
  std::vector<Polygon> holes_vector;
  for (const auto& hole : holePoints) {
    Polygon holepoly;
    bp::set_points(holepoly, hole.points().begin(), hole.points().end());
    holes_vector.push_back(holepoly);
  }
  if (holePoints.size() > 0) {
    bp::set_holes(boundary, holes_vector.begin(), holes_vector.end());
  }

  // Keep track of any orphans
  std::vector<PolygonWithHoles> orphans;
  std::vector<PolygonWithHoles> cellPolygons;
  std::vector<QuantizedPoint<2>> localGenPoints;
  std::vector<int> polyIndex;
  // Map between hashed vertices to associated polygons in cellPolygons
  std::unordered_map<QuantizedKey<2>, std::set<unsigned>, QuantizedKeyHasher<2>> vertexMap;

  // Loop over cells and intersect them with the boundary
  for (auto i = 0u; i < cells.size(); ++i) {
    auto clippedPolygon = cellPolygonIntersect(getCell(i), points[i], boundary, orphans);
    if (clippedPolygon.size() == 0) continue;
    const auto curP = cellPolygons.size();
    // Hash and map the vertices for this polygon
    for (const auto& v : clippedPolygon) {
      auto vertexHash = Q.encode(BoostToPolytope(v));
      vertexMap[vertexHash].insert(curP);
    }
    polyIndex.push_back(i);
    localGenPoints.push_back(points[i]);
    cellPolygons.push_back(clippedPolygon);
  }

  // Loop over each orphan, keep track of any newly created orphans
  std::vector<PolygonWithHoles> remainingOrphans;
  for (auto& orphan : orphans) {
    std::vector<QuantizedPoint<2>> genPoints;
    std::set<unsigned> genIndex;
    // For converting between this smaller subset of generators to the larger one
    std::unordered_map<unsigned, unsigned> smallToLarge;
    PolygonWithHoles orphanBound = orphan;
    // Grab any neighboring cells by looping over the vertices
    for (const auto& v : orphan) {
      auto vertexHash = Q.encode(BoostToPolytope(v));
      for (const auto& pi : vertexMap[vertexHash]) {
        auto [it, added] = genIndex.insert(pi);
        if (added && validUnion(cellPolygons[pi], orphanBound)) {
          smallToLarge[genPoints.size()] = pi;
          genPoints.push_back(localGenPoints[pi]);
        }
      }
    }

    // Skip if no valid neighboring cells found
    if (genPoints.empty()) continue;

    // If only 1 cell is nearby, simply add orphan to it
    if (genPoints.size() == 1) {
      cellPolygons[smallToLarge[0]] = orphanBound;
    } else {
      // Retessellate with these select generators, convert final product into boost polygons for simplicity
      QuantTessellation<2> newQT(genPoints);
      newQT.m_isEscapePod = m_isEscapePod;
      tessellator.tessellateQuantized(newQT);
      POLY_ASSERT2(newQT.cells.size() == genPoints.size(), "Number of gen points should not change");
      for (auto i = 0u; i < newQT.cells.size(); ++i) {
        auto newPolygon = cellPolygonIntersect(newQT.getCell(i), genPoints[i], orphanBound, remainingOrphans);
        cellPolygons[smallToLarge[i]] = newPolygon;
      }
    }
  }
  std::vector<std::vector<int>> newCells;
  std::vector<QuantizedPoint<2>> newNodes;
  std::vector<std::vector<unsigned>> newFaces;

  // Storage for generator points corresponding to surviving cells
  std::vector<QuantizedPoint<2>> newPoints;
  std::vector<QuantizedKey<2>> newHashes;
  std::vector<int> newCellRank;
  bool updateCellRank = cellRank.size() == 0 ? false : true;

  // Map from QuantizedPoint to index in newNodes (for vertex deduplication)
  std::map<QuantizedPoint<2>, int> node2id;

  // Map from canonical edge to face index (for oriented edge tracking)
  edge::EdgeToFaceMap edgeToFace;
  unsigned i = 0;
  for (auto& cellPoly : cellPolygons) {
    // Check if any remaining orphans can be added to the current cell
    // If they can, remove them from the remainingOrphans vector
    if (remainingOrphans.size() > 0u) {
      remainingOrphans.erase(
        std::remove_if(remainingOrphans.begin(), remainingOrphans.end(),
                       [&cellPoly](const PolygonWithHoles& pcell) {
                         if (validUnion(pcell, cellPoly)) {
                           return true;
                         }
                         return false;
                       }),
        remainingOrphans.end());
    }
    std::vector<QuantizedPoint<2>> keptVertices = bp::BoostToPolytope(cellPoly);
    removeCollinear(keptVertices);
    auto nv = keptVertices.size();
    std::vector<int> localCellIndex;
    localCellIndex.reserve(nv);
    // Gather cell indices from the node2id and update newNodes
    for (const auto& p : keptVertices) {
      auto it = node2id.find(p);
      if (it != node2id.end()) {
        localCellIndex.push_back(it->second);
      } else {
        node2id[p] = newNodes.size();
        newNodes.push_back(p);
        newNodes.back().index = node2id[p];
        localCellIndex.push_back(node2id[p]);
      }
    }
    // Build edges with oriented indexing
    std::vector<int> cellEdgeIndices;
    for (auto ic = 0u; ic < nv; ++ic) {
      auto j = (ic+1) % nv;
      int n0 = localCellIndex[ic];
      int n1 = localCellIndex[j];
      // Add oriented edge - returns signed index (negative if reversed)
      int signedFaceIndex = edge::addOrientedEdge(n0, n1, newFaces, edgeToFace);
      cellEdgeIndices.push_back(signedFaceIndex);
    }
    int curIndex = polyIndex[i++];
    newCells.push_back(cellEdgeIndices);
    newPoints.push_back(points[curIndex]);
    newHashes.push_back(hashes[curIndex]);
    if (updateCellRank) {
      newCellRank.push_back(cellRank[curIndex]);
    }
  }
  if (remainingOrphans.size() != 0) {
    ejectEscapePod("remainingorphan", QPLC, tessellator.name());
    std::cerr << "Rank " << Communicator::getRank() << ": Outstanding orphan remains" << std::endl;
    Communicator::abort();
  }
  nodes = std::move(newNodes);
  faces = std::move(newFaces);
  cells = std::move(newCells);
  points = std::move(newPoints);
  hashes = std::move(newHashes);
  if (updateCellRank) {
    cellRank = std::move(newCellRank);
  }
}
#else
// Must enable Boost to use clipping methods
template<>
void
QuantTessellation<2>::clipTessellation(const QuantPLC<2>& /*QPLC*/,
                                       Tessellator<2, double>& /*tessellator*/) {
  std::cerr << "Must compile with POLYTOPE_ENABLE_BOOST=ON and POLYTOPE_ENABLE_HIBIT2D=OFF to use clipping\n";
  Communicator::abort();
}
#endif // POLYTOPE_ENABLE_BOOST

//------------------------------------------------------------------------------
// Fill 2D tessellation mesh
//------------------------------------------------------------------------------
template<>
void
QuantTessellation<2>::fillTessellation(TessellationType& mesh) {
  auto& Q = Quantizer<2>::instance();
  compactUnusedNodesAndFaces();
  const unsigned numNodes = nodes.size();
  const unsigned numFaces = faces.size();
  const unsigned numCells = points.size();  // Number of generators

  // Allocate space for mesh data
  mesh.nodes.resize(numNodes);
  mesh.faces.resize(numFaces, std::vector<unsigned>(2));
  mesh.points.resize(numCells);
  mesh.cells = cells;
  POLY_ASSERT2(cells.size() == numCells, "Differing number of cells and generator points");

  for (unsigned i = 0; i < numCells; ++i) {
    mesh.points[i] = Q.dequantize(points[i]);
  }

  // Dequantize nodes from integer coordinates to real coordinates
  for (unsigned i = 0; i < numNodes; ++i) {
    mesh.nodes[i] = Q.dequantize(nodes[i]);
  }

  // Copy face topology (each face has 2 nodes in 2D)
  for (unsigned i = 0; i < numFaces; ++i) {
    POLY_ASSERT(faces[i].size() == 2);
    POLY_ASSERT(mesh.faces[i].size() == 2);
    mesh.faces[i][0] = faces[i][0];
    mesh.faces[i][1] = faces[i][1];
  }
  mesh.computeFaceCells();
  if (cellRank.size() == numCells) {
    mesh.cellRank = cellRank;
  }
}

//------------------------------------------------------------------------------
// Return the visible keys or points depending on the ExchangeType
// TODO: Put this in header file since it should be Dimension agnostic
//------------------------------------------------------------------------------
template<>
template<typename ExchangeType>
std::vector<ExchangeType>
QuantTessellation<2>::visibleGenerators(const std::vector<ExchangeType>& et) {
  auto N = points.size();
  std::vector<ExchangeType> result;
  if (!convexHull.m_convex) {
    makeConvexHull();
  }
  if (N == 0) return result;
  if (convexHull.isValid()) {
    result.reserve(N);
    auto plc_cell = convexHull.getCell();
    for (auto i = 0u; i < N; ++i) {
      auto qcell = getCell(i);
      if (convexBoundaryIntersect<QuantizedCoordinate<2>>(qcell, plc_cell)) {
        result.push_back(et[i]);
      }
    }
  } else {
    result = et;
  }
  return result;
}

//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class QuantTessellation<2>;

}
