//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <set>
#include <list>
#include <map>
#include <limits>

#include "polytope.hh" // Pulls in POLY_ASSERT and TriangleTessellator.hh.
#include "convexHull_2d.hh"
#include "nearestPoint.hh"
#include "within.hh"
#include "intersect.hh"
#include "PLC_Boost_2d.hh"
#include "polytope_plc_canned_geometries.hh"

// Pull in triangle. Since triangle isn't built to work out-of-the-box with C++, 
// we slurp in its source here, bracketing it with the necessary dressing.
#define TRILIBRARY
#define REAL double
#define VOID void
#define ANSI_DECLARATORS 
#define CDT_ONLY // Conforming Delaunay triangulations only! 

// Because Hang Si has messed with some of Shewchuk's predicates and
// included them with his own Tetgen library, we need to rename some of 
// the symbols therein to prevent duplicate symbols from confusing the 
// linker. Gross.
#define exactinit triangle_exactinit
#define fast_expansion_sum_zeroelim triangle_fast_expansion_sum_zeroelim
#define scale_expansion_zeroelim triangle_scale_expansion_zeroelim
#define estimate triangle_estimate
#define orient3dadapt triangle_orient3dadapt
#define orient3d triangle_orient3d
#define incircleadapt triangle_incircleadapt
#define incircle triangle_incircle

extern "C"
{
#include "triangle.h"
}


// Fast predicate for determining colinearity of points.o
extern double orient2d(double* pa, double* pb, double* pc);

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------
// A collection of helper functions
//------------------------------------------------------------------------

namespace {

//------------------------------------------------------------------------------
// Wrap the call to triangle for computing the Delaunay.
//------------------------------------------------------------------------------
template<typename RealType>
void
computeDelaunay(const vector<RealType>& points,
                triangulateio& delaunay) {

  triangulateio in;
   
  // Find the range of the generator points.
  POLY_ASSERT(points.size() % 2 == 0);
  const unsigned numGenerators = points.size()/2;

  in.numberofpoints = numGenerators;
  in.pointlist = new RealType[2*in.numberofpoints];
  copy(points.begin(), points.end(), in.pointlist);
  in.numberofsegments = 0;

  // // Add the generators
  // in.numberofpoints = numGenerators;
  // in.pointlist = new RealType[2*in.numberofpoints];
  // copy(points.begin(), points.end(), in.pointlist);
  // in.numberofsegments = 0;
    
  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0;   // No point markers.
  in.segmentmarkerlist = 0;
  in.numberofholes = 0;
  in.holelist = 0;

  // No regions.
  in.numberofregions = 0;
  in.regionlist = 0;
  
  // Set up the structure for the triangulation.
  delaunay.pointlist = 0;
  delaunay.pointattributelist = 0;
  delaunay.pointmarkerlist = 0;
  delaunay.trianglelist = 0;
  delaunay.triangleattributelist = 0;
  delaunay.neighborlist = 0;
  delaunay.segmentlist = 0;
  delaunay.segmentmarkerlist = 0;
  delaunay.edgelist = 0;
  delaunay.edgemarkerlist = 0;
  delaunay.holelist = 0;

  // Do the triangulation. Switches pass to triangle are:
  // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // -z : Indices are all numbered from zero.
  // -e : Generates edges and places them in out.edgelist.
  // -c : Generates convex hull and places it in out.segmentlist.
  // -p : Uses the given PLC information.
  // triangulate((char*)"Qzec", &in, &delaunay, 0);
  triangulate((char*)"Qz", &in, &delaunay, 0);
  
  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != in.numberofpoints) {
    char err[1024];
    snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, (int)numGenerators);
    error(err);
  }
  
  // Clean up
  delete [] in.pointlist;
}

//------------------------------------------------------------------------------
// Given an array of 3 integers and 1 unique value, find the other two.
//------------------------------------------------------------------------------
void
findOtherTriIndices(const vector<int>& indices,
                    const int a,
                    int& b,
                    int& c) {
  POLY_ASSERT(indices.size() == 3);
  POLY_ASSERT(a == indices[0] or a == indices[1] or a == indices[2]);
  POLY_ASSERT(indices[0] != indices[1] and 
              indices[0] != indices[2] and 
              indices[1] != indices[2]);
  if (a != indices[0]) {
    b = indices[0];
    c = (a != indices[1] ? indices[1] : indices[2]);
  } else {
    b = indices[1];
    c = indices[2];
  } 
}

//------------------------------------------------------------------------------
// Sort a set of edges around a face so that sequential edges share nodes.
// We account for one break in the chain, representing an unbounded surface.
//------------------------------------------------------------------------------
vector<unsigned>
computeSortedFaceNodes(const vector<pair<int, int> >& edges) {
  typedef pair<int, int> EdgeHash;
  vector<unsigned> result;

  const unsigned nedges = edges.size();
  if (nedges > 1) {

    // Invert the mapping, from nodes to edges.
    std::map<int, std::set<EdgeHash> > nodes2edges;
    internal::CounterMap<int> nodeUseCount;
    unsigned i, j;
    for (i = 0; i != nedges; ++i) {
      nodes2edges[edges[i].first].insert(edges[i]);
      nodes2edges[edges[i].second].insert(edges[i]);
      ++nodeUseCount[edges[i].first];
      ++nodeUseCount[edges[i].second];
    }
    
    // Look for any edges with one node in the set.  There can be at most
    // two such edges, representing the two ends of the chain.  We will put 
    // the edges with those nodes first in the ordering, so that all remaining
    // edges should naturally hook together.
    std::vector<EdgeHash> orderedEdges;
    orderedEdges.reserve(nedges);
    int lastNode;
    bool hangingNodes = false;
    for (i = 0; i != nedges; ++i) {
      if ((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
          (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1)) {
        POLY_ASSERT((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
                    (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1));
        orderedEdges.push_back(edges[i]);
        nodes2edges[edges[i].first].erase(edges[i]);
        nodes2edges[edges[i].second].erase(edges[i]);
        lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
        hangingNodes = true;
      }
    }
    POLY_ASSERT(orderedEdges.size() == 0 or orderedEdges.size() == 2);

    // Pick a node to start the chain.
    if (hangingNodes) {
      POLY_ASSERT(nodeUseCount[orderedEdges.back().first] == 2 or
                  nodeUseCount[orderedEdges.back().second] == 2);
      lastNode = (nodeUseCount[orderedEdges.back().first] == 2 ? 
                  orderedEdges.back().first :
                  orderedEdges.back().second);
    } else {
      lastNode = edges[0].first;
    }

    // Walk the remaining edges
    while (orderedEdges.size() != nedges) {
      POLY_ASSERT(nodes2edges[lastNode].size() > 0);
      orderedEdges.push_back(*nodes2edges[lastNode].begin());
      nodes2edges[orderedEdges.back().first].erase(orderedEdges.back());
      nodes2edges[orderedEdges.back().second].erase(orderedEdges.back());
      lastNode = (orderedEdges.back().first == lastNode ? orderedEdges.back().second : orderedEdges.back().first);
    }

    // Read the nodes in order.
    if (hangingNodes) {
      result.push_back(nodeUseCount[orderedEdges[0].first] == 1 ? orderedEdges[0].first : orderedEdges[0].second);
      result.push_back(nodeUseCount[orderedEdges[1].first] == 1 ? orderedEdges[1].first : orderedEdges[1].second);
      i = 1;
    } else {
      i = 0;
    }
    for (; i != nedges; ++i) {
      j = (i + 1) % nedges;
      POLY_ASSERT(orderedEdges[i].first == orderedEdges[j].first or
                  orderedEdges[i].first == orderedEdges[j].second or
                  orderedEdges[i].second == orderedEdges[j].first or
                  orderedEdges[i].second == orderedEdges[j].second);
      result.push_back((orderedEdges[i].first == orderedEdges[j].first or orderedEdges[i].first == orderedEdges[j].second) ? 
                       orderedEdges[i].first : 
                       orderedEdges[i].second);
    }
    POLY_ASSERT2((hangingNodes and result.size() == nedges + 1) or
                 ((not hangingNodes) and result.size() == nedges), result.size());      

  } else {

    // There are either one or no edges, so the solution is pretty simple.
    if (nedges == 1) {
      result.push_back(edges[0].first);
      result.push_back(edges[0].second);
    }
  }
  // That's it.
  return result;
}

// //------------------------------------------------------------------------------
// // Exprss the Triangle-generated delaunay struct as a Tessellation
// //------------------------------------------------------------------------------
// template<typename RealType>
// void
// constructDelaunayMesh(const triangulateio& delaunay,
//                       Tessellation<2, RealType>& dmesh) {  
//   typedef pair<int, int> EdgeHash;

//   int i;
//   dmesh.nodes.resize(2*delaunay.numberofpoints);
//   for (i = 0; i < 2*delaunay.numberofpoints; ++i) {
//     dmesh.nodes[i] = delaunay.pointlist[i];
//   }

//   int iedge, p, q, r;
//   EdgeHash pq, qr, rp;
//   map<EdgeHash, int> edge2id;
//   map<int, vector<int> > edgeCells;
//   dmesh.cells.resize(delaunay.numberoftriangles, vector<int>(3));
//   for (i = 0; i < delaunay.numberoftriangles; ++i) {
//     p = delaunay.trianglelist[3*i  ];
//     q = delaunay.trianglelist[3*i+1];
//     r = delaunay.trianglelist[3*i+2];
//     pq = internal::hashEdge(p,q);
//     qr = internal::hashEdge(q,r);
//     rp = internal::hashEdge(r,p);

//     iedge = internal::addKeyToMap(pq, edge2id);
//     edgeCells[iedge].push_back(p < q ? i : ~i);
//     dmesh.cells[i][0] = (p < q ? iedge : ~iedge);

//     iedge = internal::addKeyToMap(qr, edge2id);
//     edgeCells[iedge].push_back(q < r ? i : ~i);
//     dmesh.cells[i][2] = (q < r ? iedge : ~iedge);

//     iedge = internal::addKeyToMap(rp, edge2id);
//     edgeCells[iedge].push_back(r < p ? i : ~i);
//     dmesh.cells[i][1] = (r < p ? iedge : ~iedge);
//   }

//   dmesh.faces.resize(edge2id.size(), vector<unsigned>(2));
//   for (typename map<EdgeHash, int>::const_iterator itr = edge2id.begin();
//        itr != edge2id.end(); 
//        ++itr) {
//     const EdgeHash& edge = itr->first;
//     i = itr->second;
//     dmesh.faces[i][0] = edge.first;
//     dmesh.faces[i][1] = edge.second;
//   }

//   dmesh.faceCells.resize(edge2id.size());
//   for (i = 0; i < dmesh.faces.size(); ++i) dmesh.faceCells[i] = edgeCells[i];
// }

// //------------------------------------------------------------------------------
// // Output a Triangle-generated delaunay mesh into a silo file for viewing.
// //------------------------------------------------------------------------------
// template<typename RealType>
// void
// dumpDelaunay(const triangulateio& delaunay,
//              string name,
//              const int cycle) {
//   Tessellation<2, RealType> dmesh;
//   constructDelaunayMesh(delaunay, dmesh);
//   map<string, double*> fields;
// #ifdef POLYTOPE_ENABLE_SILO
//   SiloWriter<2,RealType>::write(dmesh, fields, fields, fields, fields, name, cycle, 0.0);
// #endif
// }

} // end anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
TriangleTessellator():
  Tessellator<2, RealType>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
~TriangleTessellator() {
}

//------------------------------------------------------------------------------
// Compute the QuantizedTessellation
//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellateQuantized(QuantizedTessellation& result) const {

  // (Re)construct the floating generator values for use in the triangularization.
  const int numGenerators = result.generators.size();
  const int numTotalGenerators = numGenerators + result.guardGenerators.size();
  vector<RealType> gencoords(2*numTotalGenerators);
  for (unsigned i = 0; i != numGenerators; ++i) {
    result.dequantize(&result.generators[i].x, &gencoords[2*i]);
  }
  for (unsigned i = 0; i != result.guardGenerators.size(); ++i) {
    result.dequantize(&result.guardGenerators[i].x, &gencoords[2*(numGenerators + i)]);
  }

  // Compute the Delaunay triangularization of the generators.
  triangulateio delaunay;
  computeDelaunay(gencoords, delaunay);
  
  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  int i, k, p, q, r;
  EdgeHash pq, pr, qr;
  RealPoint rp;
  IntPoint ppoint, qpoint, rpoint;
  result.nodes.reserve(delaunay.numberoftriangles);
  vector<set<unsigned> > gen2tri(numTotalGenerators);                         // generators -> associated triangles
  vector<vector<int> > tri2gens(delaunay.numberoftriangles, vector<int>(3));  // triangles -> associated generators
  map<EdgeHash, vector<unsigned> > dedge2tris;                                // Delauany edges -> associated triangles
  vector<int> tri2id(delaunay.numberoftriangles, -1);                         // triangles -> quantized circumcenter ID
  k = 0;
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    p = delaunay.trianglelist[3*i  ];
    q = delaunay.trianglelist[3*i+1];
    r = delaunay.trianglelist[3*i+2];
    POLY_ASSERT(p != q and p != r and q != r);
    POLY_ASSERT2(orient2d(&delaunay.pointlist[2*p],
			  &delaunay.pointlist[2*q],
			  &delaunay.pointlist[2*r]) > 0,
		 "TriangleTessellator Error: Delaunay vertices are not in CCW order for triangle " << i);
    result.quantize(&delaunay.pointlist[2*p], &ppoint.x);
    result.quantize(&delaunay.pointlist[2*q], &qpoint.x);
    result.quantize(&delaunay.pointlist[2*r], &rpoint.x);
    if ((p < numGenerators or q < numGenerators or r < numGenerators) and
        (ppoint != qpoint and qpoint != rpoint and ppoint != rpoint)) {      // Second check to ensure not a degenerate triangle
      POLY_ASSERT(result.nodes.size() == k);
      geometry::computeCircumcenter(&delaunay.pointlist[2*p],
                                    &delaunay.pointlist[2*q],
                                    &delaunay.pointlist[2*r],
                                    &rp.x);
      result.nodes.push_back(IntPoint());
      result.quantize(&rp.x, &result.nodes[k].x);
      pq = internal::hashEdge(p,q);
      pr = internal::hashEdge(p,r);
      qr = internal::hashEdge(q,r);
      dedge2tris[pq].push_back(i);
      dedge2tris[pr].push_back(i);
      dedge2tris[qr].push_back(i);
      gen2tri[p].insert(i);
      gen2tri[q].insert(i);
      gen2tri[r].insert(i);
      tri2gens[i][0] = p;
      tri2gens[i][1] = q;
      tri2gens[i][2] = r;
      tri2id[i] = k;
      ++k;
    }
  }

  // Check some stuff.
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    for (map<EdgeHash, vector<unsigned> >::const_iterator itr = dedge2tris.begin();
         itr != dedge2tris.end();
         ++itr) {
      POLY_ASSERT(itr->second.size() == 1 or itr->second.size() == 2);
    }
    for (i = 0; i != numGenerators; ++i) {
      POLY_ASSERT2(gen2tri[i].size() >= 3, i << " " << gen2tri[i].size());
    }
  }
  POLY_END_CONTRACT_SCOPE;

  // Now we can walk the Delaunay triangles associated with each (non-guard) generator
  // to build the Voronoi cells and edges.
  unsigned ii, jj, n, old_size;
  int e1;
  vector<unsigned> cellNodes;
  EdgeHash edge;
  map<std::pair<int, int>, int> edge2id;
  result.cellEdges.resize(numGenerators);
  for (p = 0; p != numGenerators; ++p) {
    const set<unsigned>& tris = gen2tri[p];
    set<EdgeHash> meshEdges;
    for (set<unsigned>::const_iterator triItr = tris.begin();
         triItr != tris.end(); 
         ++triItr) {
      i = *triItr;
      ii = tri2id[i];
      POLY_ASSERT(i < delaunay.numberoftriangles);
      POLY_ASSERT(tri2id[i] >= 0 and tri2id[i] < result.nodes.size());

      // Get the other two indices for this triangle and hash their edges wrt. p
      POLY_ASSERT(tri2gens[i].size() == 3);
      findOtherTriIndices(tri2gens[i], p, q, r);
      POLY_ASSERT(p != q and p != r);
      pq = internal::hashEdge(p,q);
      pr = internal::hashEdge(p,r);
      POLY_ASSERT2(dedge2tris[pq].size() == 2 and (dedge2tris[pq][0] == i or dedge2tris[pq][1] == i),
                   "(" << pq.first << " " << pq.second << ") " << dedge2tris[pq].size());
      POLY_ASSERT2(dedge2tris[pr].size() == 2 and (dedge2tris[pr][0] == i or dedge2tris[pr][1] == i),
                   "(" << pr.first << " " << pr.second << ") " << dedge2tris[pr].size());

      // Voronoi edge orthogonal to pq.
      k = (dedge2tris[pq][0] == i ? dedge2tris[pq][1] : dedge2tris[pq][0]);
      POLY_ASSERT(tri2id[k] >= 0 and tri2id[k] < result.nodes.size());
      jj = tri2id[k];
      if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));

      // Voronoi edge orthogonal to pr.
      k = (dedge2tris[pr][0] == i ? dedge2tris[pr][1] : dedge2tris[pr][0]);
      POLY_ASSERT(tri2id[k] >= 0 and tri2id[k] < result.nodes.size());
      jj = tri2id[k];
      if (jj != ii) meshEdges.insert(internal::hashEdge(ii,jj));
    }

    // Hook together the hashed edges around this point.  This method actually returns the
    // sorted node indicies, so we have to put the edges back together.
    cellNodes = computeSortedFaceNodes(vector<EdgeHash>(meshEdges.begin(), meshEdges.end()));
    n = cellNodes.size();
    POLY_ASSERT(n >= 3);
    for (i = 0; i != n; ++i) {
      ii = cellNodes[i];
      jj = cellNodes[(i + 1) % n];
      POLY_ASSERT(ii != jj);
      edge = internal::hashEdge(ii, jj);
      POLY_ASSERT((edge.first == ii and edge.second == jj) or
                  (edge.first == jj and edge.second == ii));
      old_size = edge2id.size();
      e1 = internal::addKeyToMap(edge, edge2id);
      if (e1 == old_size) {
        POLY_ASSERT(e1 == result.edges.size());
        result.edges.push_back(edge);
      }
      if (edge.first == ii) {
        result.cellEdges[p].push_back(e1);
      } else {
        result.cellEdges[p].push_back(~e1);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

} // end polytope namespace
