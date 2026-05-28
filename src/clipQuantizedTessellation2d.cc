//------------------------------------------------------------------------------
// 2D implementation of clipQuantizedTessellation.
//
// This method relies upon the Boost.Polygon library for geometric operations.
//------------------------------------------------------------------------------

#include <vector>
#include <map>
#include <set>

#include "clipQuantizedTessellation.hh"
#include "RegisterBoostPolygonTypes.hh"
#include "removeElements.hh"
#include "IntPointMap.hh"
#include "Tessellator.hh"
#include "QuantizedTessellation2d.hh"

using namespace std;
namespace bp = boost::polygon;
using namespace boost::polygon::operators;

namespace polytope {

namespace {
//------------------------------------------------------------------------------
// Remove collinear points from a polygon.
//------------------------------------------------------------------------------
template<typename IntType>
void
removeCollinearPoints(bp::polygon_data<IntType>& poly) {
  using namespace bp;
  typedef polygon_data<IntType> Polygon;
  typedef segment_data<IntType> Segment;
  typedef typename polygon_traits<Polygon>::point_type Point;
  typedef std::vector<polygon_data<IntType> > PolygonSet;

  // Find the unique set of points bounding the union polygon.
  const int n = poly.size();
  std::vector<Point> points;
  Segment segment;
  for (int i = 0; i != n-1; ++i) {
    segment = bp::construct<Segment>(*(poly.begin()+(i-1)%n), *(poly.begin()+(i+1)%n));
    if (!bp::contains(segment, *(poly.begin() + i))) points.push_back(*(poly.begin() + i));
  }
  POLY_ASSERT2(points.size() >= 3, points.size());

  // Construct the return value.
  poly.set(points.begin(), points.end());
}

}

//------------------------------------------------------------------------------
// The public method.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation2d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<2>& geometry,
                               const Tessellator<2, RealType>& tessellator) {

  using namespace bp;
  typedef polygon_data<IntType> Polygon;
  typedef polygon_with_holes_data<IntType> PolygonWithHoles;
  typedef typename polygon_traits<Polygon>::point_type Point;
  typedef std::vector<polygon_data<IntType> > PolygonSet;

  // { // BLAGO!
  //   // Check for any "boundary" edges on the interior.
  //   cerr << "Checking for internal boundary edges BEFORE clipping." << endl;
  //   vector<unsigned> edgeCount(qmesh.edges.size());
  //   for (unsigned i = 0; i != qmesh.cellEdges.size(); ++i) {
  //     for (unsigned j = 0; j != qmesh.cellEdges[i].size(); ++j) {
  //       ++edgeCount[internal::positiveID(qmesh.cellEdges[i][j])];
  //     }
  //   }
  //   typename QuantizedTessellation2d<IntType, RealType>::RealPoint rp1, rp2;
  //   for (unsigned i = 0; i != edgeCount.size(); ++i) {
  //     if (edgeCount[i] < 2) {
  //       qmesh.dequantize(&qmesh.nodes[qmesh.edges[i].first].x, &rp1.x);
  //       qmesh.dequantize(&qmesh.nodes[qmesh.edges[i].second].x, &rp2.x);
  //       if ((rp1.x > -0.49 and rp1.x < 0.49 and rp1.y > -0.49 and rp1.y < 0.49) or
  //           (rp2.x > -0.49 and rp2.x < 0.49 and rp2.y > -0.49 and rp2.y < 0.49)) {
  //         cerr << "Internal edge : " << i << " " << rp1 << " " << rp2 << " "
  //              << qmesh.nodes[qmesh.edges[i].first] << " " << qmesh.nodes[qmesh.edges[i].second] << endl;
  //       }
  //     }
  //   }
  // } // BLAGO!

  // Copy the input PLC geometry to a polygon Boost.Polygon will recognize.
  // We have registered a std::vector<IntPoint> as a polygon.
  const unsigned numBoundaryPoints = geometry.facets.size();
  std::vector<Point> boundaryPoints(numBoundaryPoints);
  IntPoint p;
  for (unsigned i = 0; i != geometry.facets.size(); ++i) {
    POLY_ASSERT(geometry.facets[i].size() == 2);
    const int k = geometry.facets[i][0];
    qmesh.quantize(&PLCpoints[2*k], &p.x);
    boundaryPoints[i] = bp::construct<Point>(p.x, p.y);
  }
  PolygonWithHoles boundary;
  bp::set_points(boundary, boundaryPoints.begin(), boundaryPoints.end());

  // Add the holes.
  std::vector<std::vector<Point> > holePoints(geometry.holes.size());
  for (unsigned ihole = 0; ihole != geometry.holes.size(); ++ihole) {
    const unsigned n = geometry.holes[ihole].size();
    holePoints[ihole].resize(n);
    for (unsigned i = 0; i != n; ++i) {
      POLY_ASSERT(geometry.holes[ihole][i].size() == 2);
      const int k = geometry.holes[ihole][i][0];
      qmesh.quantize(&PLCpoints[2*k], &p.x);
      holePoints[ihole][i] = bp::construct<Point>(p.x, p.y);
    }
  }
  boundary.set_holes(holePoints.begin(), holePoints.end());
  
  // Make the polygon_set boundary.
  PolygonSet boundarySet;
  boundarySet += boundary;

  // Prepare new mesh data.
  const unsigned ncells = qmesh.cellEdges.size();
  std::vector<bp::IntPoint> newNodes;
  std::vector<std::pair<int, int> > newEdges;
  std::vector<std::vector<int> > newCellEdges(ncells);
  internal::IntPointMap<2, IntType> node2id(2);
  typedef std::map<std::pair<int, int>, int> EdgeIDMap;
  EdgeIDMap edge2id;
  std::vector<std::set<int> > newNodeCells;
  std::vector<Polygon> cellPolygons;
  
  // We also prepare to keep track of the orphans for each cell.
  typedef PolygonSet Orphanage;
  Orphanage orphans;

  // Now walk each cell and clip it with the boundary.
  for (unsigned i = 0; i != ncells; ++i) {
    unsigned nverts = qmesh.cellEdges[i].size();

    // Copy the cell geometry to the quantized polygon type.
    Polygon cell;
    std::vector<Point> cellPoints(nverts);
    for (unsigned j = 0; j != nverts; ++j) {
      int k = qmesh.cellEdges[i][j];
      if (k < 0) {
        k = ~k;
        p = qmesh.nodes[qmesh.edges[k].second];
      } else {
        p = qmesh.nodes[qmesh.edges[k].first];
      }
      cellPoints[j] = bp::construct<Point>(p.x, p.y);
    }
    bp::set_points(cell, cellPoints.begin(), cellPoints.end());

    // Clip the cell against the boundary.
    PolygonSet cellSet;
    cellSet += cell;
    cellSet &= boundarySet;

    // Check if we generated more than one polygon.  If so, some are orphans.
    Point gen = bp::construct<Point>(qmesh.generators[i].x, qmesh.generators[i].y);
    unsigned polygonIndex = 0;
    if (cellSet.size() > 1) {

      // This intersection generated orphans.  Find the polygon that contains the generator -- that becomes this new cell.
      // The other fragments we insert into the orphanage for adoption.
      while (polygonIndex < cellSet.size() and not bp::contains(cellSet[polygonIndex], gen)) ++polygonIndex;
      POLY_ASSERT(polygonIndex < cellSet.size());

      // For each new orphan check if it can be unioned with an existing orphan.  Otherwise add it to the orphan set.
      unsigned k;
      PolygonSet trialunion;
      typename QuantizedTessellation2d<IntType, RealType>::RealPoint rp;
      IntType pp[2] = {gen.x(), gen.y()};
      qmesh.dequantize(pp, &rp.x);
      for (unsigned ipoly = 0; ipoly != cellSet.size(); ++ipoly) {
        if (ipoly != polygonIndex) {
          for (k = 0; k != orphans.size(); ++k) {
            PolygonSet set1, set2, trialunion;
            set1 += cellSet[ipoly];
            set2 += orphans[k];
            bp::assign(trialunion, set1 | set2);
            if (trialunion.size() == 1) {
              // Looks like these orphans can be unioned.  Replace the exising orphan with the union.
              orphans[k] = trialunion[0];
              removeCollinearPoints(orphans[k]);
              break;
            }
          }

          // If we didn't find another orphan to union with, add this one to the orphanage.
          if (k == orphans.size()) orphans += cellSet[ipoly];
        }
      }
    }

    // Read out the final cell geometry to the new QuantizedTessellation.
    // It appears Boost.Polygon gives us back the same node for the beginning and ending nodes of the
    // polygon.
    POLY_ASSERT(bp::contains(cellSet[polygonIndex], gen));
    cellPolygons.push_back(cellSet[polygonIndex]);
    nverts = cellSet[polygonIndex].size();
    POLY_ASSERT(nverts > 3);
    POLY_ASSERT(*(cellSet[polygonIndex].begin()) == *(cellSet[polygonIndex].begin() + nverts - 1));
    for (unsigned j = 0; j != nverts-1; ++j) {
      const Point& v0 = *(cellSet[polygonIndex].begin() + j);
      const Point& v1 = *(cellSet[polygonIndex].begin() + j + 1);
   
      // Insert vertex 0.
      p.x = v0.x();
      p.y = v0.y();
      const int j0 = node2id.index(p);
      if (j0 == newNodes.size()) {
        newNodes.push_back(p);
        newNodeCells.push_back(std::set<int>());
      }
      newNodeCells[j0].insert(i);

      // Insert vertex 1.
      p.x = v1.x();
      p.y = v1.y();
      const int j1 = node2id.index(p);
      if (j1 == newNodes.size()) {
        newNodes.push_back(p);
        newNodeCells.push_back(std::set<int>());
      }
      newNodeCells[j1].insert(i);
    
      // Now insert the edge if non-degenerate.
      if (j0 != j1) {
        const std::pair<int, int> edge = internal::hashEdge(j0, j1);
        POLY_ASSERT((edge.first == j0 and edge.second == j1) or
                    (edge.first == j1 and edge.second == j0));
        const int old_size = edge2id.size();
        const int e1 = internal::addKeyToMap(edge, edge2id);
        if (e1 == old_size) {
          POLY_ASSERT(e1 == newEdges.size());
          newEdges.push_back(edge);
        }
        if (edge.first == j0) {
          newCellEdges[i].push_back(e1);
        } else {
          newCellEdges[i].push_back(~e1);
        }
      }
    }
  }
  POLY_ASSERT(cellPolygons.size() == ncells);

  // Now deal with any orphans.
  for (typename Orphanage::iterator orphanItr = orphans.begin();
       orphanItr != orphans.end();
       ++orphanItr) {
    const Polygon& orphan = *orphanItr;

    // Start the local boundary with just the orphan.
    PolygonSet localBoundary;
    localBoundary += orphan;

    // Walk the edges of this orphan and find the neighboring cells.
    // We use those neigbors to augment the local boundary by union.
    std::set<unsigned> neighbors;
    for (typename Polygon::iterator_type pItr = orphan.begin();
         pItr != orphan.end() - 1;
         ++pItr) {

      // Check the generators touching this node.
      const IntPoint p0(pItr->x(), pItr->y());
      if (node2id.have(p0)) {
        const int nodeID = node2id.index(p0);
        POLY_ASSERT(nodeID < newNodeCells.size());
        neighbors.insert(newNodeCells[nodeID].begin(), newNodeCells[nodeID].end());
      }
    }
    POLY_ASSERT(neighbors.size() > 0);

    // Find the overall boundary of the orphan unioned with all the neighbors, and accumulate the
    // local generator positions.
    vector<IntPoint> localgenerators;
    std::set<unsigned>::const_iterator gitr = neighbors.begin();
    for (unsigned i = 0; i != neighbors.size(); ++i, ++gitr) {
      const unsigned igen = *gitr;
      PolygonSet oldcell;
      oldcell += cellPolygons[igen];
      bp::assign(localBoundary, localBoundary | oldcell);
      localgenerators.push_back(qmesh.generators[igen]);
      localgenerators.back().index = i;
    }
    POLY_ASSERT(localBoundary.size() == 1);

    // Construct the tessellation of the neighbor generators.
    QuantizedTessellation2d<IntType, RealType> localqmesh(localgenerators, qmesh);
    tessellator.tessellateQuantized(localqmesh);

    // Clip each new cell geometry against the local boundary to obtain the final cell geometries.
    gitr = neighbors.begin();
    for (unsigned i = 0; i != localgenerators.size(); ++i, ++gitr) {
      const unsigned igen = *gitr; // localgenerators[i].index;
      Polygon cell;
      unsigned nverts = localqmesh.cellEdges[i].size();
      std::vector<Point> cellPoints(nverts);
      for (unsigned j = 0; j != nverts; ++j) {
        int k = localqmesh.cellEdges[i][j];
        if (k < 0) {
          k = ~k;
          p = localqmesh.nodes[localqmesh.edges[k].second];
        } else {
          p = localqmesh.nodes[localqmesh.edges[k].first];
        }
        typename QuantizedTessellation2d<IntType, RealType>::RealPoint rp;
        qmesh.dequantize(&p.x, &rp.x);
        cellPoints[j] = bp::construct<Point>(p.x, p.y);
      }
      bp::set_points(cell, cellPoints.begin(), cellPoints.end());

      // Clip the cell against the boundary.
      PolygonSet cellSet;
      cellSet += cell;
      cellSet &= localBoundary;
      POLY_ASSERT2(cellSet.size() == 1, cellSet.size());  // Hopefully no more new orphans!

      removeCollinearPoints(cellSet[0]); // Get rid of any collinear points the intersection operation may have left behind.

      // Read out the new cell geometry.
      cellPolygons[igen] = cellSet[0];
      newCellEdges[igen].clear();
      nverts = cellSet[0].size();
      for (unsigned j = 0; j != nverts; ++j) {
        const Point& v0 = *(cellSet[0].begin() + j);
        const Point& v1 = *(cellSet[0].begin() + ((j + 1) % nverts));
   
        // Insert vertex 0.
        p.x = v0.x();
        p.y = v0.y();
        const int j0 = node2id.index(p);
        if (j0 == newNodes.size()) {
          newNodes.push_back(p);
          newNodeCells.push_back(std::set<int>());
        }
        newNodeCells[j0].insert(igen);

        // Insert vertex 1.
        p.x = v1.x();
        p.y = v1.y();
        const int j1 = node2id.index(p);
        if (j1 == newNodes.size()) {
          newNodes.push_back(p);
          newNodeCells.push_back(std::set<int>());
        }
        newNodeCells[j1].insert(igen);
    
        // Now insert the edge if non-degenerate.
        if (j0 != j1) {
          const std::pair<int, int> edge = internal::hashEdge(j0, j1);
          POLY_ASSERT((edge.first == j0 and edge.second == j1) or
                      (edge.first == j1 and edge.second == j0));
          const int old_size = edge2id.size();
          const int e1 = internal::addKeyToMap(edge, edge2id);
          if (e1 == old_size) {
            POLY_ASSERT(e1 == newEdges.size());
            newEdges.push_back(edge);
          }
          if (edge.first == j0) {
            newCellEdges[igen].push_back(e1);
          } else {
            newCellEdges[igen].push_back(~e1);
          }
        }
      }
    }
  }
  POLY_ASSERT(newCellEdges.size() == ncells);

  // If we dealt with any orphans, there may be unused edges and nodes we should clear out before
  // copying the final topology.
  if (orphans.size() > 0) {

    // Count how many times each edge is used.
    vector<unsigned> edgeUseCount(newEdges.size(), 0);
    for (unsigned i = 0; i != ncells; ++i) {
      const unsigned nedges = newCellEdges[i].size();
      for (unsigned j = 0; j != nedges; ++j) {
        const int iedge = internal::positiveID(newCellEdges[i][j]);
        POLY_ASSERT(iedge < newEdges.size());
        ++edgeUseCount[iedge];
      }
    }

    // Look for unused edges, and build up the mapping of new edge IDs.
    vector<unsigned> edges2kill;
    for (unsigned i = 0; i != newEdges.size(); ++i) {
      POLY_ASSERT(edgeUseCount[i] == 0 or edgeUseCount[i] == 1 or edgeUseCount[i] == 2);
      if (edgeUseCount[i] == 0) edges2kill.push_back(i);
    }
    const vector<int> old2new_edgeIDs = internal::removeElements(newEdges, edges2kill);

    // Patch up the cell->edge indexing.
    for (unsigned i = 0; i != ncells; ++i) {
      const unsigned nedges = newCellEdges[i].size();
      for (unsigned j = 0; j != nedges; ++j) {
        if (newCellEdges[i][j] < 0) {
          POLY_ASSERT(old2new_edgeIDs[~newCellEdges[i][j]] != -1);
          newCellEdges[i][j] = ~old2new_edgeIDs[~newCellEdges[i][j]];
        } else {
          POLY_ASSERT(old2new_edgeIDs[newCellEdges[i][j]] != -1);
          newCellEdges[i][j] = old2new_edgeIDs[newCellEdges[i][j]];
        }
      }
    }

    // Count how many times each node is used.
    vector<unsigned> nodeUseCount(newNodes.size(), 0);
    for (unsigned i = 0; i != newEdges.size(); ++i) {
      ++nodeUseCount[newEdges[i].first];
      ++nodeUseCount[newEdges[i].second];
    }

    // Now remove any unused nodes.
    vector<unsigned> nodes2kill;
    for (unsigned i = 0; i != newNodes.size(); ++i) {
      if (nodeUseCount[i] == 0) nodes2kill.push_back(i);
    }
    const vector<int> old2new_nodeIDs = internal::removeElements(newNodes, nodes2kill);

    // Patch the edge->node indexing.
    for (unsigned i = 0; i != newEdges.size(); ++i) {
      POLY_ASSERT(old2new_nodeIDs[newEdges[i].first] != -1);
      POLY_ASSERT(old2new_nodeIDs[newEdges[i].second] != -1);
      newEdges[i].first = old2new_nodeIDs[newEdges[i].first];
      newEdges[i].second = old2new_nodeIDs[newEdges[i].second];
    }
  }

  // { // BLAGO!
  //   // Dump the current state of the mesh.
  //   typename QuantizedTessellation2d<IntType, RealType>::RealPoint rp;
  //   cerr << "********************************************************************************" << endl
  //        << "AFTER REMOVING UNUSED EDGES" << endl
  //        << "Nodes:" << endl;
  //   for (unsigned i = 0; i != newNodes.size(); ++i) {
  //     qmesh.dequantize(&newNodes[i].x, &rp.x);
  //     cerr << "    " << i << " " << newNodes[i] << " " << rp << endl;
  //   }
  //   cerr << "Edges:" << endl;
  //   for (unsigned i = 0; i != newEdges.size(); ++i) {
  //     cerr << "    " << i << " (" << newEdges[i].first << " " << newEdges[i].second << ") == [";
  //     qmesh.dequantize(&newNodes[newEdges[i].first].x, &rp.x);
  //     cerr << rp << " ";
  //     qmesh.dequantize(&newNodes[newEdges[i].second].x, &rp.x);
  //     cerr << rp << "]" << endl;
  //   }
  //   cerr << "Cells: " << endl;
  //   for (unsigned i = 0; i != newCellEdges.size(); ++i) {
  //     cerr << "    " << i << " Edges={";
  //     for (unsigned j = 0; j != newCellEdges[i].size(); ++j) cerr << internal::positiveID(newCellEdges[i][j]) << " ";
  //     cerr << "]" << endl
  //          << "    ";
  //     for (unsigned j = 0; j != newCellEdges[i].size(); ++j) {
  //       if (newCellEdges[i][j] < 0) {
  //         qmesh.dequantize(&newNodes[newEdges[~newCellEdges[i][j]].second].x, &rp.x);
  //       } else {
  //         qmesh.dequantize(&newNodes[newEdges[newCellEdges[i][j]].first].x, &rp.x);
  //       }
  //       cerr << " " << rp;
  //     }
  //     cerr << endl;
  //   }
  //   cerr << "********************************************************************************" << endl;
  // } // BLAGO!

  // Copy the topology and geometry back to the input QuantizedTessellation.
  qmesh.nodes = newNodes;
  qmesh.edges = newEdges;
  qmesh.cellEdges = newCellEdges;

  // { // BLAGO!
  //   // Check for any "boundary" edges on the interior.
  //   cerr << "Checking for internal boundary edges AFTER clipping." << endl;
  //   vector<unsigned> edgeCount(qmesh.edges.size());
  //   for (unsigned i = 0; i != qmesh.cellEdges.size(); ++i) {
  //     for (unsigned j = 0; j != qmesh.cellEdges[i].size(); ++j) {
  //       ++edgeCount[internal::positiveID(qmesh.cellEdges[i][j])];
  //     }
  //   }
  //   typename QuantizedTessellation2d<IntType, RealType>::RealPoint rp1, rp2;
  //   for (unsigned i = 0; i != edgeCount.size(); ++i) {
  //     if (edgeCount[i] < 2) {
  //       qmesh.dequantize(&qmesh.nodes[qmesh.edges[i].first].x, &rp1.x);
  //       qmesh.dequantize(&qmesh.nodes[qmesh.edges[i].second].x, &rp2.x);
  //       if ((rp1.x > -0.49 and rp1.x < 0.49 and rp1.y > -0.49 and rp1.y < 0.49) or
  //           (rp2.x > -0.49 and rp2.x < 0.49 and rp2.y > -0.49 and rp2.y < 0.49)) {
  //         cerr << "Internal edge : " << i << " " << rp1 << " " << rp2 << " "
  //              << qmesh.nodes[qmesh.edges[i].first] << " " << qmesh.nodes[qmesh.edges[i].second] << endl;
  //       }
  //     }
  //   }
  // } // BLAGO!

}

//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
template
void
clipQuantizedTessellation<int, double>(QuantizedTessellation2d<int, double>& qmesh,
                                       const std::vector<double>& PLCpoints,
                                       const PLC<2>& geometry,
                                       const Tessellator<2, double>& tessellator);

}
