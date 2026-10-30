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

template<>
void
QuantTessellation<2>::cullExternalPoints(const QuantPLC<2>& QPLC) {
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
    newHashes.push_back(m_Q.encode(point));
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
  bp::point_data<bp::BoostCoordinate2D> genPoint = bp::polytopeToBoost(genPointP);
  std::vector<PolygonWithHoles> cellSet = boostIntersect(currentCell, boundary);
  unsigned fragIndex = 0;
  auto NFrag = cellSet.size();
  if (NFrag == 0) return PolygonWithHoles(); // Cell was completely outside boundary
  // Find the fragment that owns the generator.
  while (fragIndex < NFrag and
         not bp::contains(cellSet[fragIndex], genPoint)) ++fragIndex;
  POLY_ASSERT2(fragIndex < NFrag,
               "No clipped fragment contains generator " << genPointP
               << " (" << NFrag << " fragments)");
  if (NFrag > 1) {
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
  const auto boundaryPoints = bp::polytopeToBoostPoints(QPLC.getCell());
  PolygonWithHoles boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());
  auto holePoints = QPLC.getHolePoints();
  std::vector<Polygon> holes_vector;
  for (const auto& hole : holePoints) {
    Polygon holepoly;
    const auto holeBoostPoints = bp::polytopeToBoostPoints(hole);
    bp::set_points(holepoly, holeBoostPoints.begin(), holeBoostPoints.end());
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
      auto vertexHash = m_Q.encode(BoostToPolytope(v));
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
      auto vertexHash = m_Q.encode(BoostToPolytope(v));
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
      QuantTessellation<2> newQT(m_Q, genPoints);
      newQT.escapePod = escapePod;
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
  EdgeToFaceMap edgeToFace;
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
    POLY_ASSERT2(nv >= 3,
                 "Clipped cell " << i << " has only " << nv
                 << " vertices after conversion to quantized coordinates");
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
      int signedFaceIndex = addOrientedEdge(n0, n1, newFaces, edgeToFace);
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
  // FIXME: Implement escape pod with QPLC
  // if (remainingOrphans.size() != 0) {
  //   ejectEscapePod("remainingorphan", QPLC, tessellator.name());
  //   std::cerr << "Rank " << Communicator::getRank() << ": Outstanding orphan remains" << std::endl;
  //   Communicator::abort();
  // }
  POLY_ASSERT2(newCells.size() == newPoints.size(),
               "Clipped cell and generator counts differ: "
               << newCells.size() << " cells, " << newPoints.size() << " generators");
  POLY_ASSERT2(newCells.empty() or !newNodes.empty(),
               "Clipped tessellation has cells but no nodes");
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
  compactUnusedNodesAndFaces();
  const unsigned numNodes = nodes.size();
  const unsigned numFaces = faces.size();
  const unsigned numCells = points.size();  // Number of generators

  POLY_ASSERT2(numCells == 0 or numNodes > 0,
               "Tessellation has " << numCells << " cells but no nodes");
  // Allocate space for mesh data
  mesh.nodes.resize(numNodes);
  mesh.faces.resize(numFaces, std::vector<unsigned>(2));
  mesh.points.resize(numCells);
  mesh.cells = cells;
  POLY_ASSERT2(cells.size() == numCells, "Differing number of cells and generator points");
  for (auto cellIndex = 0u; cellIndex < cells.size(); ++cellIndex) {
    POLY_ASSERT2(cells[cellIndex].size() >= 3,
                 "Cell " << cellIndex << " has fewer than three faces");
  }

  for (unsigned i = 0; i < numCells; ++i) {
    mesh.points[i] = m_Q.dequantize(points[i]);
  }

  // Dequantize nodes from integer coordinates to real coordinates
  for (unsigned i = 0; i < numNodes; ++i) {
    mesh.nodes[i] = m_Q.dequantize(nodes[i]);
  }

  // Copy face topology (each face has 2 nodes in 2D)
  for (unsigned i = 0; i < numFaces; ++i) {
    POLY_ASSERT(faces[i].size() == 2);
    POLY_ASSERT(mesh.faces[i].size() == 2);
    POLY_ASSERT2(faces[i][0] < numNodes and faces[i][1] < numNodes and
                 faces[i][0] != faces[i][1],
                 "Face " << i << " has invalid node indices");
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
      if (convexBoundaryIntersect<QCoord>(qcell, plc_cell)) {
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
