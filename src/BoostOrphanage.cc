//------------------------------------------------------------------------
// BoostOrphanage
//------------------------------------------------------------------------
#include <set>
#include <map>
#include <vector>

#include "polytope.hh"
#include "BoostOrphanage.hh"
#include "PLC_Boost_2d.hh"

namespace polytope {

using namespace std;


//------------------------------------------------------------------------
// A collection of helper functions
//------------------------------------------------------------------------

namespace {



//------------------------------------------------------------------------
// Union a Boost.Geometry ring with a Boost.Geometry multi_polygon.
// The resulting multi_polygon is corrected to ensure it conforms to
// the proper geometric concept.
//------------------------------------------------------------------------
template<typename IntType>
void
createBGUnion(boost::geometry::model::ring<Point2<IntType>,false> ring, 
              boost::geometry::model::multi_polygon<
                boost::geometry::model::polygon<Point2<IntType>,false> >& multiPolygon ) {
  boost::geometry::model::multi_polygon<
    boost::geometry::model::polygon<Point2<IntType>,false> > temp;
  boost::geometry::union_(multiPolygon, ring, temp);
  boost::geometry::correct(temp);
  multiPolygon=temp;
}

//------------------------------------------------------------------------
// Union may produce cell rings containing faces broken into multiple
// pieces by redundant hanging nodes. This manifests as three (or more)
// collinear nodes. We clip them as a final connection step
//------------------------------------------------------------------------
template<typename IntType>
void
cleanRingEdges(boost::geometry::model::ring<Point2<IntType>,false>& ring) {
  // Pre-conditions
  POLY_ASSERT(ring.size() > 2);
  POLY_ASSERT(ring.front() == ring.back());
  //POLY_ASSERT(!boost::geometry::intersects(ring));

  // Initialize the temporary ring
  typedef boost::geometry::model::ring<Point2<IntType>,false> BGring;
  typedef IntType CoordHash;
  BGring tmpRing;
  
  // Check collinearity on the middle points
  bool collinear;
  for (typename BGring::const_iterator itr = ring.begin()+1;
       itr != ring.end()-1; 
       ++itr) {
    collinear = geometry::collinear<2,CoordHash>(&(*(itr-1)).x,
                                                 &(*(itr  )).x,
                                                 &(*(itr+1)).x,
                                                 1);
    if(!collinear) boost::geometry::append(tmpRing,*itr);
  }
  
  // Check the beginning/end point
  collinear = geometry::collinear<2,CoordHash>(&(*(ring.end()-2  )).x,
                                               &(*(ring.begin()  )).x,
                                               &(*(ring.begin()+1)).x,
                                               1);
  if(!collinear) boost::geometry::append(tmpRing,*(ring.begin()));

  // Close the ring
  boost::geometry::correct(tmpRing);

  // Post-conditions
  POLY_ASSERT(tmpRing.size() > 2);
  POLY_ASSERT(tmpRing.front() == tmpRing.back());
  //POLY_ASSERT(!boost::geometry::intersects(tmpRing));
  ring = tmpRing;
}

//------------------------------------------------------------------------
// Remove any collinear adjacent facets from a reduced PLC
//------------------------------------------------------------------------
template<typename RealType>
void
removeCollinearPointsFromPLC(ReducedPLC<2, RealType>& plc, const RealType tol) {
  ReducedPLC<2, RealType> result;
  const int nFacets = plc.facets.size();
  for (int i = 0; i < nFacets; ++i) {
    const int i1 = i;
    const int i2 = (i+1) % nFacets;
    POLY_ASSERT(plc.facets[i1].size() == 2 and plc.facets[i2].size() == 2);
    POLY_ASSERT(plc.facets[i1][1] == plc.facets[i2][0]);
    const int k1 = plc.facets[i1][0];
    const int k2 = plc.facets[i1][1];
    const int k3 = plc.facets[i2][1];
    POLY_ASSERT(k1 != k2 and k1 != k3 and k2 != k3);
    if (not geometry::collinear<2, RealType>(&plc.points[2*k1],
                                             &plc.points[2*k2],
                                             &plc.points[2*k3],
                                             tol)) {
      result.points.push_back(plc.points[2*k2  ]);
      result.points.push_back(plc.points[2*k2+1]);
    }
  }
  POLY_ASSERT(result.points.size()/2 > 0 and result.points.size()/2 <= nFacets);
  result.facets.resize(result.points.size()/2, std::vector<int>(2));
  for (int i = 0; i < result.points.size()/2; ++i) {
    result.facets[i][0] = i;
    result.facets[i][1] = (i+1) % (result.points.size()/2);
  }
  plc = result;
}

//------------------------------------------------------------------------
// Collapse a facet having points within a given tolerance
//------------------------------------------------------------------------
template<typename RealType>
void
simplifyPLC(ReducedPLC<2, RealType>& plc, const RealType tol) {
  typedef Point2<RealType> PointType;
  ReducedPLC<2, RealType> result;
  const int nFacets = plc.facets.size();
  for (int i = 0; i < nFacets; ++i) {
    POLY_ASSERT(plc.facets[i].size() == 2);
    const int i1 = plc.facets[i][0];
    const int i2 = plc.facets[i][1];
    POLY_ASSERT(i1 < plc.points.size()/2 and i2 < plc.points.size()/2);
    const PointType p1 = PointType(plc.points[2*i1], plc.points[2*i1+1]);
    const PointType p2 = PointType(plc.points[2*i2], plc.points[2*i2+1]);
    const double dx = double(p2.x - p1.x);
    const double dy = double(p2.y - p1.y);
    const double dist = sqrt(dx*dx + dy*dy);
    if (dist > tol) {
      result.points.push_back(plc.points[2*i1]);
      result.points.push_back(plc.points[2*i1+1]);
    }
  }
  POLY_ASSERT(result.points.size()/2 > 0 and result.points.size()/2 <= nFacets);
  result.facets.resize(result.points.size()/2, std::vector<int>(2));
  for (int i = 0; i < result.points.size()/2; ++i) {
    result.facets[i][0] = i;
    result.facets[i][1] = (i+1) % (result.points.size()/2);
  }
  plc = result;
}

//------------------------------------------------------------------------
// Extract the point from a 2D ReducedPLC as a polytope::Point2
//------------------------------------------------------------------------
template<typename RealType>
Point2<RealType>
getPoint(const ReducedPLC<2, RealType> plc, const int i) {
  POLY_ASSERT(i < plc.points.size()/2);
  typedef Point2<RealType> PointType;
  return PointType(plc.points[2*i], plc.points[2*i+1]);
}

//------------------------------------------------------------------------
// Extract the point from a 2D ReducedPLC as a polytope::Point2
//------------------------------------------------------------------------
template<typename RealType, typename IntType>
void
constructBoostVoronoiCells(boost::polygon::voronoi_diagram<RealType>& voronoi,
                           const vector<Point2<IntType> >& intPoints,
                           const Point2<RealType>& low,
                           const Point2<RealType>& high,
                           const RealType degeneracy) {
  typedef Point2<RealType> RealPoint;
  typedef Point2<IntType>  IntPoint;
  typedef boost::polygon::voronoi_diagram<RealType> VD;
  const int numGenerators = intPoints.size();
  const RealPoint center = (low + high)*0.5;
  int sortedIndex=0, cellIndex;
  IntPoint node, direction, pinf;
  RealPoint endpoint;
  map<IntPoint, int> node2id;
  map<int, IntPoint> id2node;
  // vector<unsigned> cellNodes(numGenerators);
  for (typename VD::const_cell_iterator cellItr = voronoi.cells().begin(); 
       cellItr != voronoi.cells().end(); 
       ++cellItr, ++sortedIndex) {
    const typename VD::edge_type* edge = cellItr->incident_edge();
    vector<unsigned> nodeChain;
    do {
      cellIndex = intPoints[sortedIndex].index;
      const typename VD::vertex_type* v0 = edge->vertex0();
      const typename VD::vertex_type* v1 = edge->vertex1();

      // Finite edge
      if (v0 and v1) {
        node = IntPoint(IntType(v0->x()), IntType(v0->y()));
        const unsigned old_size = node2id.size();
        const unsigned j = internal::addKeyToMap(node, node2id);
        nodeChain.push_back(j);
        if (j == old_size) id2node[j] = node;
      }
      
      // Infinite edge
      else {
        POLY_ASSERT(v0 or v1);
        const typename VD::vertex_type* vfin = v0 ? v0 : v1;
        node = IntPoint(IntType(vfin->x()), IntType(vfin->y()));
        
        // Determine the edge direction pointing to infinity
        const typename VD::cell_type* cell1 = edge->cell();
        const typename VD::cell_type* cell2 = edge->twin()->cell();
        POLY_ASSERT(cell1->contains_point() and cell2->contains_point());
        const size_t index1 = cell1->source_index();
        const size_t index2 = cell2->source_index();
        POLY_ASSERT(index1 < numGenerators and index2 < numGenerators);
        const unsigned cellIndex1 = intPoints[index1].index;
        const unsigned cellIndex2 = intPoints[index2].index;
        
        // Floating point endpoint
        const IntPoint r = IntPoint(intPoints[2*cellIndex2  ] - intPoints[2*cellIndex1  ],
                                    intPoints[2*cellIndex2+1] - intPoints[2*cellIndex1+1]);
        const IntPoint d = v0 ? IntPoint(-r.y, r.x) : IntPoint(r.y, -r.x);
        endpoint.x = center.x + (RealType(node.x) + 0.5)*degeneracy;
        endpoint.y = center.y + (RealType(node.y) + 0.5)*degeneracy;
        
        
        //TODO
        //TODO Set endpoint equal to node
        //TODO Compute a direction vector to infinity
        //TODO Project endpoint along direction vector to bounding box (i.e. +/- coordMax/2)
        //TODO
        
        
        // Vertex 0 is finite, vertex 1 is the projected node. Add them in order
        if (v0) {
          { // Vertex 0
            const unsigned old_size = node2id.size();
            const unsigned j = internal::addKeyToMap(node, node2id);
            nodeChain.push_back(j);
            if (j == old_size)  id2node[j] = node;
          }

          { // Vertex 1
            node = pinf;
            const unsigned old_size = node2id.size();
            const unsigned j = internal::addKeyToMap(node, node2id);
            nodeChain.push_back(j);
            if (j == old_size)  id2node[j] = node;
          }
        }
  
        // Vertex 0 is the projected infNode. Only add vertex 0.
        else {
          node = pinf;
          const unsigned old_size = node2id.size();
          const unsigned j = internal::addKeyToMap(node, node2id);
          nodeChain.push_back(j);
          if (j == old_size)  id2node[j] = node;
        }
      }

      edge = edge->next();
    } while (edge != cellItr->incident_edge());
    POLY_ASSERT(not nodeChain.empty());

    // Remove repeated node indices in the chain
    vector<unsigned>::iterator it = std::unique(nodeChain.begin(), nodeChain.end());
    nodeChain.resize(std::distance(nodeChain.begin(), it));
    if (nodeChain.front() == nodeChain.back()) nodeChain.resize(nodeChain.size()-1);
    POLY_ASSERT(not nodeChain.empty());

    // cellNodes[cellIndex] = nodeChain;
  }
}



} // end anonymous namespace





//------------------------------------------------------------------------------
template<typename RealType>
BoostOrphanage<RealType>::
BoostOrphanage(const Tessellator<2, RealType>* tessellator):
  OrphanageBase<2, RealType>(tessellator) {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
BoostOrphanage<RealType>::
~BoostOrphanage() {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<typename RealType>
void
BoostOrphanage<RealType>::
adoptOrphans(const vector<RealType>& points,
             const QuantizedCoordinates<2, RealType>& coords,
             vector<IntPLC>& intCells,
             vector<IntPLC>& orphans) const {
  // Pre-conditions
  POLY_ASSERT(not points.empty());
  POLY_ASSERT(not intCells.empty());
  POLY_ASSERT(not orphans.empty());
  POLY_ASSERT(points.size()/2 == intCells.size());

  // Create reduced set of orphans by unioning all that neighbor each other.
  vector<IntPLC> reducedOrphans = BG::boost_unionReduce(orphans);
  POLY_ASSERT2(reducedOrphans.size() > 0 and reducedOrphans.size() <= orphans.size(), 
               reducedOrphans.size() << " " << orphans.size());

  // Construct map from node points to cells
  std::map<IntPoint, std::set<int> > point2cells;
  for (int j = 0; j < intCells.size(); ++j) {
    for (int k = 0; k < intCells[j].points.size()/2; ++k) {
      point2cells[getPoint<CoordHash>(intCells[j], k)].insert(j);
    }
  }

  // Compute the adoption on each orphan in the reduced set
  for (int iorphan = 0; iorphan < reducedOrphans.size(); ++iorphan) {
    IntPLC orphan = reducedOrphans[iorphan];
    POLY_ASSERT(not orphan.facets.empty());
    POLY_ASSERT(not orphan.points.empty());
    POLY_ASSERT2(not BG::boost_intersects(orphan),
                 "Unioning orphans has produced self-intersections");
  
    // Determine neighbors of this orphan
    std::set<int> orphanNeighbors;
    std::map<IntPoint, std::set<int> > point2neighbor;
    for (int k = 0; k < orphan.points.size()/2; ++k) {
      IntPoint ip = getPoint<CoordHash>(orphan, k);
      typename std::map<IntPoint, std::set<int> >::iterator it = point2cells.find(ip);
      if (it != point2cells.end()) {
        orphanNeighbors.insert(it->second.begin(), it->second.end());
        point2neighbor[it->first].insert(it->second.begin(), it->second.end());
      }
    }
    POLY_ASSERT(not orphanNeighbors.empty());

    // If the orphan only has a single neighbor, we can skip a lot of work.
    // No need to tessellate - simply union the orphan with its neighbor cell.
    std::vector<IntPLC> nbh_cells;
    if (orphanNeighbors.size() > 1) {

      // Organize a list of the orphan and its neighboring cells
      vector<IntPLC> neighborhood;
      std::vector<RealType> nbh_points;
      neighborhood.push_back(orphan);
      for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
           nbItr != orphanNeighbors.end();
           ++nbItr) {
        nbh_points.push_back(points[2*(*nbItr)  ]);
        nbh_points.push_back(points[2*(*nbItr)+1]);
        neighborhood.push_back(intCells[*nbItr]);
      }
 
      // Union reduce the orphan and its neighboring cells
      neighborhood = BG::boost_unionReduce(neighborhood);
      POLY_ASSERT(neighborhood.size() == 1);
      POLY_ASSERT2(neighborhood.size() > 0, "Union produced empty set!");

      // Union the orphan with the neighborhood one more time as a fail-safe
      //
      // NOTE: The order in which we union the orphan with each neighbor matters.
      //       If the orphan and a neighboring cell only share a single node, and
      //       that union is computed first, boost.geometry will return two 
      //       polygons. The orphan ring may become lost in subsequent unions.
      neighborhood.push_back(orphan);
      neighborhood = BG::boost_unionReduce(neighborhood);
      POLY_ASSERT(neighborhood.size() == 1);

      // Some informative output if the neighborhood has more than one member
      if (neighborhood.size() > 1) {
        std::cerr << "Blago!" << std::endl;
        for (int i = 0; i != neighborhood.size(); ++i) {
          std::cerr << "Polygon " << i << " in the union contains" << std::endl
                    << neighborhood[i] << std::endl;
        }
        POLY_ASSERT(false);
      }


      POLY_BEGIN_CONTRACT_SCOPE;
      {
        for (int ipoint = 0; ipoint < neighborhood[0].points.size()/2; ++ipoint) {
          bool result = false;
          const IntPoint ip = getPoint<CoordHash>(neighborhood[0], ipoint);
          for (int jpoint = 0; jpoint < orphan.points.size()/2; ++jpoint) {
            const IntPoint jp = getPoint<CoordHash>(orphan, jpoint);
            result += (ip == jp);
          }
          for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
               nbItr != orphanNeighbors.end(); 
               ++nbItr) {
            for (int jpoint = 0; jpoint < intCells[*nbItr].points.size()/2; ++jpoint) {
              const IntPoint jp = getPoint<CoordHash>(intCells[*nbItr], jpoint);
              result += (ip == jp);
            }
          }
          POLY_ASSERT2(result, "Union error: there are points in the bounding "
                       << "neighborhood that are not in the orphan or its neighbor cells.");
        }
      }
      POLY_END_CONTRACT_SCOPE;

      IntPLC nbh_boundary = neighborhood[0];
      removeCollinearPointsFromPLC<CoordHash>(nbh_boundary, 1);
      POLY_ASSERT(nbh_boundary.facets.size() > 2);
      POLY_ASSERT(not BG::boost_intersects(nbh_boundary));

      this->callPrivateTessellate(nbh_points, nbh_boundary, coords, nbh_cells);
    }


    // We're only concerned with the cells in the sub-tessellation whose generators
    // are immediate neighbors of the orphaned chunk. These are the only cells which can
    // "adopt" the orphan based on a local Voronoi principle
    for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
         nbItr != orphanNeighbors.end(); 
         ++nbItr){
      std::set<int>::iterator it = orphanNeighbors.find(*nbItr);
      POLY_ASSERT(it != orphanNeighbors.end());
      
      // Index for the sub-tessellation
      int nbh_index = distance(orphanNeighbors.begin(), it);
      POLY_ASSERT(nbh_index < orphanNeighbors.size());

      // Index for the full tessellation
      int index = *nbItr;
      POLY_ASSERT(index < points.size()/2);
      
      // The full tessellation's cell
      IntPLC cell;

      if (orphanNeighbors.size() > 1) {
        cell = nbh_cells[nbh_index];
        simplifyPLC<CoordHash>(cell, 1);

        // // Blago!
        // cerr << "Old Cell " << index << ":" << endl << intCells[index] << endl;
        // cerr << "New Cell " << index << ":" << endl << cell << endl;
        // // Blago!

        // When calling the private tessellate routine, we may end up with new cell
        // rings that vertices that are close to the original ones, but off by a few
        // quantized grid spacings. Unioning the old and new rings together can produce
        // self-intersections and inconsistency with cells unaffected by the orphanage.
        // Solution: Walk the new ring and make sure its vertices are consistent with
        // the old ring.
        IntPLC newCell;
        newCell.points.resize(cell.points.size());
        for (unsigned i = 0; i < cell.facets.size(); ++i) {
          POLY_ASSERT(i < cell.points.size()/2);
          const IntPoint ip = IntPoint(cell.points[2*i], cell.points[2*i+1]);
          
          // First check for points that are exactly equal
          bool result = false;
          unsigned j = 0;
          IntPoint jp;
          while (j < intCells[index].facets.size() and not result) {
            jp = IntPoint(intCells[index].points[2*j], intCells[index].points[2*j+1]);
            result += (ip == jp);
            ++j;
          }

          // One point was exatly equal
          if (result) {
            newCell.points[2*i  ] = jp.x;
            newCell.points[2*i+1] = jp.y;
          }

          // No equal points
          else {
            RealPoint pOld, pNew = coords.dequantize(&ip.x);
            RealType dist = std::numeric_limits<RealType>::max();
            unsigned j = 0;
            IntPoint jp;
            const RealType tol = 10.0*coords.delta;
            while (j < intCells[index].facets.size() and dist > tol) {
              jp = IntPoint(intCells[index].points[2*j], intCells[index].points[2*j+1]);
              pOld = coords.dequantize(&jp.x);
              dist = geometry::distance<2, RealType>(&pOld.x, &pNew.x);
              ++j;
            }

            // A point is very close by. This is likely the point we were aiming for.
            // Replace the new cell vertex with the old vertex location.
            if (dist < tol) {
              newCell.points[2*i  ] = jp.x;
              newCell.points[2*i+1] = jp.y;
            }

            // No points are close by. This is likely a new cell vertex resulting from
            // the orphan, so we'll just add it.
            else {
              newCell.points[2*i  ] = ip.x;
              newCell.points[2*i+1] = ip.y;
            }
          }
        }
        cell.points = newCell.points;

        // // Blago!
        // cerr << "New New Cell " << index << ":" << endl << cell << endl;
        // // Blago!
      }

      // If the orphan has only a single neighbor, just compute its union with
      // that neighbor's cell ring from the full tessellation
      else {
        cell = orphan;
  
        vector<IntPLC> unionCell = BG::boost_union(cell, intCells[index]);
        if (unionCell.size() > 1) {
          std::cerr << "Blago!" << std::endl << "Unioned cell " << index
                    << " has more than one cell:" << std::endl;
          for (int j = 0; j < unionCell.size(); ++j) {
            std::cerr << std::endl << "Cell " << j << ":" << std::endl << unionCell[j] << std::endl;
          }
        }
        POLY_ASSERT(unionCell.size() == 1);
        cell = unionCell[0];

        POLY_BEGIN_CONTRACT_SCOPE;
        {
          for (int i = 0; i < cell.facets.size(); ++i) {
            bool result = false;
            const IntPoint ip = IntPoint(cell.points[2*i], cell.points[2*i+1]);
            for (int j = 0; j < orphan.facets.size(); ++j) {
              const IntPoint jp = IntPoint(orphan.points[2*j], orphan.points[2*j+1]);
              result += (ip == jp);
            }
            for (int j = 0; j < intCells[index].facets.size(); ++j) {
              const IntPoint jp = IntPoint(intCells[index].points[2*j], intCells[index].points[2*j+1]);
              result += (ip == jp);
            }
            if (not result) {
              cerr << "\nCELL:\n" << cell
                   << "\nORPHAN:\n" << orphan
                   << "\nCURRENT NEIGHBOR " << index << ":\n" << intCells[index] << endl;
            }
            POLY_ASSERT2(result, "Union error: there are points in the single-neighbor "
                         << "union cell that are not in the orphan or its neighbor.");
          }
        }
        POLY_END_CONTRACT_SCOPE;
      }
      POLY_ASSERT(cell.facets.size() > 2);
      
      // Union may produce cell rings containing edges broken into multiple
      // pieces (i.e. three or more sequential collinear nodes).
      removeCollinearPointsFromPLC<CoordHash>(cell, 1);
      POLY_ASSERT(cell.facets.size() > 2);
      POLY_ASSERT(not BG::boost_intersects(cell));

      // Replace the corresponding cell ring
      intCells[index] = cell;
    }
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class BoostOrphanage<double>;

} //end polytope namespace






// //------------------------------------------------------------------------------
// template<typename RealType>
// void
// BoostOrphanage<RealType>::
// adoptOrphans_OLD(const vector<RealType>& points,
//                  const QuantizedCoordinates<2, RealType>& coords,
//                  vector<BGring>& cellRings,
//                  vector<BGring>& orphans) const {
//   // Pre-conditions
//   POLY_ASSERT(!points.empty());
//   POLY_ASSERT(!cellRings.empty());
//   POLY_ASSERT(points.size()/2 == cellRings.size());
//   POLY_ASSERT(!orphans.empty());
  
//   int i,j;
  
//   // Create reduced set of orphans by unioning all that neighbor each other.
//   // A Boost.Geometry multi-polygon does this operation seemlessly
//   BGmulti_polygon reducedOrphanSet;
//   for (i = 0; i != orphans.size(); ++i) {
//     createBGUnion(orphans[i], reducedOrphanSet);
//   }  
  
//   // std::cerr << "Number of orphans = " << reducedOrphanSet.size() << std::endl;

//   // Compute the adoption on each orphan in the reduced set
//   for (i = 0; i < reducedOrphanSet.size(); ++i) {
//     BGring orphan = reducedOrphanSet[i].outer();
//     POLY_ASSERT(!orphan.empty());
//     POLY_ASSERT(orphan.front() == orphan.back());
//     POLY_ASSERT2(!boost::geometry::intersects(orphan),
//                  "Unioning orphans has produced self-intersections");
    
//     // // Blago!
//     // std::cerr << "Orphan " << i << ":" << std::endl;
//     // for (typename BGring::const_iterator itr = orphan.begin();
//     //      itr != orphan.end(); 
//     //      ++itr)  std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//     // std::cerr << std::endl;
//     // // Blago!

//     // Construct map from node points to neighboring cells
//     std::map<IntPoint, std::set<int> > point2neighbors;
//     for (j = 0; j != cellRings.size(); ++j) {
//       for (typename BGring::const_iterator itr = cellRings[j].begin();
// 	   itr != cellRings[j].end() - 1; ++itr) {
// 	point2neighbors[*itr].insert(j);
//       }
//     }
    
//     // Determine neighbors
//     std::set<int> orphanNeighbors;
//     for (typename BGring::const_iterator pointItr = orphan.begin();
//          pointItr != orphan.end()-1; 
//          ++pointItr) {
//       std::map<IntPoint, std::set<int> >::iterator it = 
//          point2neighbors.find(*pointItr);
//       if (it != point2neighbors.end()) {
//         std::set<int> neighborSet = it->second;
//         orphanNeighbors.insert(neighborSet.begin(), neighborSet.end());
//       }
//     }
//     POLY_ASSERT(orphanNeighbors.size() > 0);
    
//     // If the orphan only has a single neighbor, we can skip a lot of work.
//     // No need to tessellate - simply union the orphan with its neighbor cell.
//     std::vector<BGring> subCellRings;
//     if (orphanNeighbors.size() > 1) {
      
//       // Union the orphan and its neighbors to get a bounded polygon
//       BGmulti_polygon neighborCells;
//       createBGUnion(orphan, neighborCells);
//       POLY_ASSERT(neighborCells.size() == 1);
      
//       // // Blago!
//       // std::cerr << "Initial union polygon (just the orphan):" << std::endl;
//       // for (typename BGring::const_iterator itr = neighborCells[0].outer().begin();
//       //      itr != neighborCells[0].outer().end(); 
//       //      ++itr)  std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//       // std::cerr << std::endl;
//       // // Blago!

//       // Collect the neighboring generators
//       std::vector<RealType> subpoints;
//       for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
//            nbItr != orphanNeighbors.end(); 
//            ++nbItr) {
//         subpoints.push_back( points[2*(*nbItr)  ] );
//         subpoints.push_back( points[2*(*nbItr)+1] );
//         createBGUnion(cellRings[*nbItr], neighborCells);
//       }
//       POLY_ASSERT(neighborCells.size() == 1);
//       POLY_ASSERT2( neighborCells.size() > 0, "Union produced empty set!" );

//       // As a fail-safe, union with the orphan one last time.
//       //
//       // NOTE: The order in which we union the orphan with each neighbor matters.
//       //       If the orphan and a neighboring cell only share a single node, and
//       //       that union is computed first, boost.geometry will return two 
//       //       polygons. The orphan ring may become lost in subsequent unions.
//       createBGUnion(orphan, neighborCells);
//       POLY_ASSERT(neighborCells.size() == 1);

//       // std::cerr << "Number of neighbors = " << orphanNeighbors.size() << std::endl;
      
//       // Throw an error if the bounding polygon has more than one piece
//       if (neighborCells.size() > 1) {
//         std::cerr << "Blago!" << std::endl;
//         for (i = 0; i != neighborCells.size(); ++i){
//           std::cerr << "Polygon " << i << " in the union contains" << std::endl;
//           for (typename BGring::const_iterator itr = neighborCells[i].outer().begin();
//                itr != neighborCells[i].outer().end(); 
//                ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//           POLY_ASSERT(false);
//         }
//       }

//       // // Blago!
//       // std::cerr << "Bounding Union Polygon" << std::endl;
//       // for (typename BGring::const_iterator itr = neighborCells[0].outer().begin();
//       //      itr != neighborCells[0].outer().end(); 
//       //      ++itr)  std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//       // std::cerr << std::endl;
//       // // Blago!
      
// #ifndef NDEBUG
//       for (typename BGring::const_iterator itr = neighborCells[0].outer().begin();
// 	   itr != neighborCells[0].outer().end() - 1; 
// 	   ++itr ) {
// 	bool result = false;
// 	for (typename BGring::const_iterator pItr = orphan.begin();
// 	     pItr != orphan.end() - 1;
// 	     ++pItr) result += (*itr == *pItr);
// 	for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
// 	     nbItr != orphanNeighbors.end(); 
// 	     ++nbItr) {
// 	  for (typename BGring::const_iterator pItr = cellRings[*nbItr].begin();
// 	       pItr != cellRings[*nbItr].end() - 1;
// 	       ++pItr) result += (*itr == *pItr);
// 	}
// 	POLY_ASSERT2(result, "Union error: there are points in the bounding ring "
// 		     << "that are not in the orphan or any of its neighbors.");
//       }
// #endif
      
//       BGring boundaryRing = neighborCells[0].outer();
//       //BGpolygon boundaryPolygon = neighborCells[0];
//       cleanRingEdges(boundaryRing);
//       POLY_ASSERT(boundaryRing.size()  >  2);
//       POLY_ASSERT(boundaryRing.front() == boundaryRing.back());
//       POLY_ASSERT(!boost::geometry::intersects(boundaryRing));

//       // Extract the vertices of the bounding polygon
//       std::vector<CoordHash> subIntPLCpoints;
//       int nSides = 0;
//       for (typename BGring::const_iterator itr = boundaryRing.begin();
//            itr != boundaryRing.end() - 1; ++itr, ++nSides) {
//          subIntPLCpoints.push_back((*itr).x);
//          subIntPLCpoints.push_back((*itr).y);
//       }
      
//       // Form the bounding PLC
//       PLC<2> subPLC;
//       subPLC.facets.resize(nSides, std::vector<int>(2) );
//       for (unsigned ii = 0; ii < nSides; ++ii) {
//         subPLC.facets[ii][0] = ii;
//         subPLC.facets[ii][1] = (ii+1) % nSides;
//       }
      
//       // Tessellate this sub-region
//       vector<vector<vector<CoordHash> > > IntCells;
//       // this->callPrivateTessellate(subpoints, subIntPLCpoints, subPLC, coords, IntCells);
      
//       subCellRings.resize(IntCells.size());
//       for (unsigned ii = 0; ii != IntCells.size(); ++ii) {
//         vector<IntPoint> cellNodes;
//         for (unsigned jj = 0; jj != IntCells[ii].size(); ++jj) {
//           POLY_ASSERT(IntCells[ii][jj].size() == 2);
//           cellNodes.push_back(IntPoint(IntCells[ii][jj][0], IntCells[ii][jj][1]));
//         }
//         boost::geometry::assign(subCellRings[ii], BGring(cellNodes.begin(),
//                                                          cellNodes.end()));
//         POLY_ASSERT(IntCells[ii].size()      == subCellRings[ii].size());
//         POLY_ASSERT(subCellRings[ii].front() == subCellRings[ii].back());
//       }
      
      
//       // // Blago!
//       // std::cerr << std::endl << "SubCellRings:" << std::endl;
//       // for (int ii = 0; ii != subCellRings.size(); ++ii) {
//       //   std::cerr << "  Cell " << ii << std::endl;
//       //   for (typename BGring::const_iterator itr = subCellRings[ii].begin();
//       //        itr != subCellRings[ii].end(); 
//       //        ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//       //   std::cerr << std::endl;
//       // }
//       // // Blago!
//     }
    
    
//     // We're only concerned with the cells in the sub-tessellation whose generators
//     // are immediate neighbors of the orphaned chunk. These are the only cells which can
//     // "adopt" the orphan based on a local Voronoi principle
//     for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
//          nbItr != orphanNeighbors.end(); ++nbItr){
//       std::set<int>::iterator it = orphanNeighbors.find(*nbItr);
//       POLY_ASSERT(it != orphanNeighbors.end());
      
//       // Index for the sub-tessellation
//       int subIndex = distance(orphanNeighbors.begin(), it);  
//       POLY_ASSERT(subIndex < orphanNeighbors.size());
      
//       // Index for the full tessellation
//       int thisIndex = *nbItr;                                
//       POLY_ASSERT(thisIndex < points.size()/2);
      
//       // The full tessellation's cell ring
//       BGring thisRing;
      
//       if (orphanNeighbors.size() > 1) {
//         thisRing = subCellRings[subIndex];
        
//         // // Blago!
//         // std::cerr << "\nCell " << thisIndex << " before fixup" << std::endl;
//         // for (typename BGring::const_iterator itr = thisRing.begin();
//         //      itr != thisRing.end();
//         //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//         // std::cerr << "\nCurrent Cell " << thisIndex << std::endl;
//         // for (typename BGring::const_iterator itr = cellRings[thisIndex].begin();
//         //      itr != cellRings[thisIndex].end();
//         //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//         // // Blago!


//         // When calling the private tessellate routine, we may end up with new cell
//         // rings that vertices that are close to the original ones, but off by a few
//         // quantized grid spacings. Unioning the old and new rings together can produce
//         // self-intersections and inconsistency with cells unaffected by the orphanage.
//         // Solution: Walk the new ring and make sure its vertices are consistent with
//         // the old ring.
//         BGring newRing;
//         std::vector<IntPoint> consistentPoints;
//         for (typename BGring::const_iterator itr = thisRing.begin();
//              itr != thisRing.end() - 1;
//              ++itr) {

//           // First check for points that are exactly equal
//           bool result = false;
//           typename BGring::iterator oItr = cellRings[thisIndex].begin();
//           while (oItr != cellRings[thisIndex].end() and !result) {
//             result += (*itr == *oItr);
//             ++oItr;
//           }
//           if (result) consistentPoints.push_back(*(oItr-1));
          
//           // No equal points. Check for points that are VERY close by. 
//           else {
//             RealPoint pOld, pNew = coords.dequantize(&(*itr).x);
//             RealType dist = std::numeric_limits<RealType>::max();
//             oItr = cellRings[thisIndex].begin();
//             while (oItr != cellRings[thisIndex].end() and dist > 15.0*coords.delta) {
//               pOld = coords.dequantize(&(*oItr).x);
//               dist = geometry::distance<2,RealType>(&pOld.x, &pNew.x);
//               ++oItr;
//             }

//             // No points are close by. This is likely a new ring vertex resulting
//             // from the orphan, so we'll just add it.
//             if (oItr == cellRings[thisIndex].end())  consistentPoints.push_back(*itr     );

//             // A point is very close by. This is likely the point we were aiming for.
//             // Replace the new ring vertex with the old vertex location.
//             else                                     consistentPoints.push_back(*(oItr-1));
//           }
//         }
//         consistentPoints.push_back(consistentPoints[0]);
//         boost::geometry::assign(newRing, BGring(consistentPoints.begin(),
//                                                 consistentPoints.end()));
//         boost::geometry::unique(newRing);
//         thisRing = newRing;

//         // // Blago!
//         // std::cerr << "\nCell " << thisIndex << " after fixup" << std::endl;
//         // for (typename BGring::const_iterator itr = thisRing.begin();
//         //      itr != thisRing.end();
//         //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//         // // Blago!
//       }
      
//       // If the orphan has only a single neighbor, just compute its union with
//       // that neighbor's cell ring from the full tessellation
//       else{
//         thisRing = orphan;
	
// 	// Union orphan with its one neighbor
// 	std::vector<BGring> unionRing;
// 	boost::geometry::union_(thisRing, cellRings[thisIndex], unionRing);
// 	if (unionRing.size() > 1) {
// 	  std::cerr << "Blago!" << std::endl << "Cell " << thisIndex
// 		    << " has more than one cell ring:" << std::endl;
// 	  for (i = 0; i < unionRing.size(); ++i) {
// 	    std::cerr << std::endl << "Ring " << i << ":" << std::endl;
// 	    for (typename BGring::const_iterator itr = unionRing[i].begin();
// 		 itr != unionRing[i].end(); 
// 		 ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
// 	  }
// 	}
// 	POLY_ASSERT(unionRing.size() == 1);
// 	thisRing = unionRing[0];

// #ifndef NDEBUG
//         for (typename BGring::const_iterator itr = thisRing.begin();
//              itr != thisRing.end() - 1; 
//              ++itr ) {
//           bool result = false;
//           for (typename BGring::const_iterator pItr = orphan.begin();
//                pItr != orphan.end() - 1;
//                ++pItr) result += (*itr == *pItr);
//           for (typename BGring::const_iterator pItr = cellRings[thisIndex].begin();
// 	       pItr != cellRings[thisIndex].end() - 1;
// 	       ++pItr) result += (*itr == *pItr);
//           POLY_ASSERT2(result, "Union error: there are points in the single-neighbor "
//                        << "union ring that are not in the orphan or its neighbor.");
//         }
// #endif
//       }
//       POLY_ASSERT(thisRing.size()  >  2);
//       POLY_ASSERT(thisRing.front() == thisRing.back());
            
//       // Union may produce cell rings containing edges broken into multiple
//       // pieces (i.e. three or more sequential collinear nodes).
//       cleanRingEdges(thisRing);
//       POLY_ASSERT(thisRing.size()  >  2);
//       POLY_ASSERT(thisRing.front() == thisRing.back());
//       POLY_ASSERT2(!boost::geometry::intersects(thisRing), "Cell " << thisIndex);

//       // // Blago!
//       // std::cerr << "\nCell " << thisIndex << " after final cleaning" << std::endl;
//       // for (typename BGring::const_iterator itr = thisRing.begin();
//       //      itr != thisRing.end();
//       //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
//       // // Blago!
      
//       // Replace the corresponding cell ring
//       cellRings[thisIndex] = thisRing;
//     }
//   }
// }
// //------------------------------------------------------------------------------
