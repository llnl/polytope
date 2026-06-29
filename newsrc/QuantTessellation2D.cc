//-----------------------------------------------------------------------------//
// QuantTessellation<2> specializations
//-----------------------------------------------------------------------------//

#include "QuantTessellation.hh"
#include "Intersections.hh"
#include "Tessellator.hh"
#include "EdgeUtils.hh"
#include <map>
#include <set>
#include <algorithm>
#include "RegisterBoostPolygonTypes.hh"
#include "boost/polygon/polygon_set_data.hpp"

namespace polytope {

namespace bp = boost::polygon;
// Need this to use the -=, +=, etc operators
using namespace boost::polygon::operators;

namespace {
using IntType = HashKey<2>::IntType;
using PolygonWithHoles = bp::polygon_with_holes_data<IntType>;
using PolygonSet = bp::polygon_set_data<IntType>;
void printpoint(const Quantizer<2>& Q,
                const Point2<HashKey<2>::IntType>& point) {
  auto qp = Q.dequantize(point);
  std::cout << qp << "," << std::endl;
}

void printPolygon(const Quantizer<2>& Q,
                  const PolygonWithHoles& polygon) {
  auto hole_points = bp::innerPoints(polygon);
  auto points = bp::outerPoints(polygon);
  std::cout << "v = [";
  for (const auto& point : points) {
    printpoint(Q, point);
  }
  std::cout << "vertices.append(v)" << std::endl;
  for (const auto& hole : hole_points) {
    std::cout << "hole " << std::endl;
    for (const auto& p : hole) {
      printpoint(Q, p);
    }
  }
}
}

// Remove any generator points that are outside our clipping region
template<>
void
QuantTessellation<2>::cullExternalPoints(const QuantPLC<2>& QPLC) {
  auto boundaryPoints = QPLC.getFacetPoints();
  PolygonWithHoles boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());
  auto holePoints = QPLC.getHolePoints();
  std::vector<Polygon> holes_vector;
  for (const auto& hole : holePoints) {
    Polygon holepoly;
    bp::set_points(holepoly, hole.begin(), hole.end());
    holes_vector.push_back(holepoly);
  }
  if (holePoints.size() > 0) {
    bp::set_holes(boundary, holes_vector.begin(), holes_vector.end());
  }
  auto N = m_points.size();
  std::vector<IntPoint> newPoints;
  newPoints.reserve(N);
  std::vector<CoordHash> newHashes;
  newHashes.reserve(N);
  unsigned indx = 0;
  for (int i = 0; i < N; ++i) {
    auto point = m_points[i];
    bp::point_data<IntType> genPoint = bp::construct<IntPoint>(point.x, point.y);
    if (!bp::contains(boundary, genPoint)) continue;
    point.index = indx++;
    newPoints.push_back(point);
    newHashes.push_back(m_Q.hash(point));
  }
  m_points = std::move(newPoints);
  m_hashes = std::move(newHashes);
}

template<>
void
QuantTessellation<2>::clipTessellation(const QuantPLC<2>& QPLC,
                                       const Tessellator<2, double>& tessellator) {
  auto boundaryPoints = QPLC.getFacetPoints();
  PolygonWithHoles boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());
  auto holePoints = QPLC.getHolePoints();
  std::vector<Polygon> holes_vector;
  for (const auto& hole : holePoints) {
    Polygon holepoly;
    bp::set_points(holepoly, hole.begin(), hole.end());
    holes_vector.push_back(holepoly);
  }
  if (holePoints.size() > 0) {
    bp::set_holes(boundary, holes_vector.begin(), holes_vector.end());
  }

  // Keep track of any orphans
  std::vector<PolygonWithHoles> orphans;
  std::vector<PolygonWithHoles> cellPolygons;
  std::vector<Point2<IntType>> localGenPoints;
  std::vector<int> polyIndex;
  // Map between hashed vertices to associated polygons in cellPolygons
  std::unordered_map<CoordHash, std::set<unsigned>> vertexMap;

  // Loop over cells and intersect them with the boundary
  for (auto i = 0; i < m_cells.size(); ++i) {
    bp::point_data<IntType> genPoint = bp::construct<IntPoint>(m_points[i].x, m_points[i].y);
    std::vector<PolygonWithHoles> cellSet = boostIntersect(getCell(i), boundary);
    auto NFrag = cellSet.size();
    if (NFrag == 0) continue; // Cell was completely outside the boundary
    // Check if orphans were generated
    unsigned fragIndex = 0;
    if (NFrag > 1) {
      // Find which part owns the generator
      while (fragIndex < NFrag and not bp::contains(cellSet[fragIndex], genPoint)) ++fragIndex;
      for (unsigned iPoly = 0; iPoly < NFrag; ++iPoly) {
        if (iPoly != fragIndex) {
          // Check if orphan can be added to other orphans
          bool foundUnion = false;
          for (auto& orphan : orphans) {
            auto trialUnion = boostUnion(orphan, cellSet[iPoly]);
            if (trialUnion.size() == 1) {
              orphan = trialUnion[0];
              foundUnion = true;
              break;
            }
          }
          if (!foundUnion) {
            orphans.push_back(cellSet[iPoly]);
          }
        }
      }
    } // End of multiple fragment check
    const auto curP = cellPolygons.size();
    // Hash and map the vertices for this polygon
    for (const auto& v : cellSet[fragIndex]) {
      auto vertexHash = m_Q.hash(BoostToPolytope(v));
      vertexMap[vertexHash].insert(curP);
    }
    polyIndex.push_back(i);
    localGenPoints.push_back(m_points[i]);
    cellPolygons.push_back(cellSet[fragIndex]);
  }

  // Loop over each orphan
  for (auto& orphan : orphans) {
    std::vector<Point2<IntType>> genPoints;
    std::set<unsigned> genIndex;
    // For converting between this smaller subset of generators to the larger one
    std::unordered_map<unsigned, unsigned> smallToLarge;
    PolygonSet orphanBounds;
    orphanBounds += orphan;
    // Grab any neighboring cells by looping over the vertices
    for (const auto& v : orphan) {
      auto vertexHash = m_Q.hash(BoostToPolytope(v));
      for (const auto& pi : vertexMap[vertexHash]) {
        auto [it, added] = genIndex.insert(pi);
        if (added) {
          orphanBounds |= cellPolygons[pi];
          smallToLarge[genPoints.size()] = pi;
          genPoints.push_back(localGenPoints[pi]);
        }
      }
    }

    // Skip if no valid neighboring cells found
    if (genPoints.empty()) continue;

    // If only 1 cell is nearby, simply add orphan to it
    std::vector<PolygonWithHoles> orphanBound;
    orphanBounds.get(orphanBound);
    POLY_ASSERT2(orphanBound.size() == 1, "Should only have 1 polygon as the boundary for the orphans");
    if (genPoints.size() == 1) {
      cellPolygons[smallToLarge[0]] = orphanBound[0];
    } else {
      // Retessellate with these select generators, convert final product into boost polygons for simplicity
      QuantTessellation<2> newQT(genPoints, *this);
      tessellator.tessellateQuantized(QPLC, newQT);
      POLY_ASSERT2(newQT.m_cells.size() == genIndex.size(), "Number of gen points should not change");
      for (auto i = 0; i < newQT.m_cells.size(); ++i) {
        std::vector<PolygonWithHoles> newPolygon = boostIntersect(newQT.getCell(i), orphanBound[0]);
        POLY_ASSERT2(newPolygon.size() > 0, "Final clipped polygon must exist");
        POLY_ASSERT2(newPolygon.size() == 1, "Only one polygon per generator after merging and clipping");
        cellPolygons[smallToLarge[i]] = newPolygon[0];
      }
    }
  }
  std::vector<std::vector<int>> newCells;
  std::vector<IntPoint> newNodes;
  std::vector<std::vector<int>> newFaces;

  // Storage for generator points corresponding to surviving cells
  std::vector<IntPoint> newPoints;
  std::vector<CoordHash> newHashes;

  // Map from IntPoint to index in newNodes (for vertex deduplication)
  std::map<IntPoint, int> nodeMap;

  // Map from canonical edge to face index (for oriented edge tracking)
  edge::EdgeToFaceMap edgeToFace;
  unsigned i = 0;
  for (auto& cellPoly : cellPolygons) {
    int curIndex = polyIndex[i];
    std::vector<IntPoint> keptVertices = bp::BoostToPolytope(cellPoly);
    //removeCollinear(keptVertices);
    auto nv = keptVertices.size();
    std::vector<int> localCellIndex;
    localCellIndex.reserve(nv);
    // Gather cell indices from the nodeMap and update newNodes
    for (const auto& p : keptVertices) {
      auto it = nodeMap.find(p);
      if (it != nodeMap.end()) {
        localCellIndex.push_back(it->second);
      } else {
        nodeMap[p] = newNodes.size();
        newNodes.push_back(p);
        newNodes.back().index = nodeMap[p];
        localCellIndex.push_back(nodeMap[p]);
      }
    }
    // Build edges with oriented indexing
    std::vector<int> cellEdgeIndices;
    for (auto ic = 0; ic < nv; ++ic) {
      auto j = (ic+1) % nv;
      int n0 = localCellIndex[ic];
      int n1 = localCellIndex[j];
      // Add oriented edge - returns signed index (negative if reversed)
      int signedFaceIndex = edge::addOrientedEdge(n0, n1, newFaces, edgeToFace);
      cellEdgeIndices.push_back(signedFaceIndex);
    }
    newCells.push_back(cellEdgeIndices);
    newPoints.push_back(m_points[curIndex]);
    newHashes.push_back(m_hashes[curIndex]);
    i++;
  }
  m_nodes = std::move(newNodes);
  m_faces = std::move(newFaces);
  m_cells = std::move(newCells);
  m_points = std::move(newPoints);
  m_hashes = std::move(newHashes);
}

//------------------------------------------------------------------------------
// Fill 2D tessellation mesh
//------------------------------------------------------------------------------
template<>
void
QuantTessellation<2>::fillTessellation(TessellationType& mesh) const {
  const unsigned numNodes = m_nodes.size();
  const unsigned numFaces = m_faces.size();
  const unsigned numCells = m_points.size();  // Number of generators

  // Allocate space for mesh data
  // In 2D: nodes are stored as [x0, y0, x1, y1, ...]
  mesh.nodes.resize(2 * numNodes);
  mesh.faces.resize(numFaces, std::vector<unsigned>(2));
  mesh.faceCells.resize(numFaces);
  mesh.points.resize(2 * numCells);
  mesh.cells = m_cells;
  POLY_ASSERT2(m_cells.size() == numCells, "Differing number of cells and generator points");

  for (unsigned i = 0; i < numCells; ++i) {
    RealPoint rp = m_Q.dequantize(m_points[i]);
    mesh.points[2*i]   = rp.x;
    mesh.points[2*i+1] = rp.y;
  }

  // Dequantize nodes from integer coordinates to real coordinates
  for (unsigned i = 0; i < numNodes; ++i) {
    RealPoint rp = m_Q.dequantize(m_nodes[i]);
    mesh.nodes[2*i]     = rp.x;
    mesh.nodes[2*i + 1] = rp.y;
  }

  // Copy face topology (each face has 2 nodes in 2D)
  for (unsigned i = 0; i < numFaces; ++i) {
    POLY_ASSERT(m_faces[i].size() == 2);
    mesh.faces[i][0] = m_faces[i][0];
    mesh.faces[i][1] = m_faces[i][1];
  }

  // Build faceCells connectivity: for each cell, mark which faces touch it
  // Cells stored signed face indices where negative means inverted orientation.
  for (unsigned i = 0; i < numCells; ++i) {
    const unsigned nf = mesh.cells[i].size();
    for (unsigned j = 0; j != nf; ++j) {
      auto k = mesh.cells[i][j];
      if (k < 0) {
        // Negative index: inverted face orientation
        POLY_ASSERT2(~k < numFaces, k << " " << ~k << " " << numFaces);
        mesh.faceCells[~k].push_back(~i);
      } else {
        // Positive index: normal face orientation
        POLY_ASSERT2(k < numFaces, k << " " << numFaces);
        mesh.faceCells[k].push_back(i);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class QuantTessellation<2>;

}
