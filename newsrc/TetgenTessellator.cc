
#include "TetgenTessellator.hh"
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

namespace polytope {

void TetgenTessellator::tessellateQuantized(const QuantPLC<3>& qplc,
                                            QT& result) const {
  tetgenio in = createTetgenPoints(qplc, result);
  tetgenio out;
  tetrahedralize((char*)"pzqQv", &in, &out);

  // Convert vorout to QuantTessellation format
  convertVoronoiToQuantTessellation(out, result);
}

void TetgenTessellator::tessellateQuantized(QT& result) const {
  tetgenio in = createTetgenPoints(result);
  tetgenio out;
  tetrahedralize((char*)"pzqQv", &in, &out);

  // Convert vorout to QuantTessellation format
  convertVoronoiToQuantTessellation(out, result);
}

void TetgenTessellator::convertVoronoiToQuantTessellation(
    const tetgenio& vorout, QT& result) const {

  result.clear();

  // 1. Add all Voronoi vertices (nodes) - these are already computed circumcenters
  POLY_ASSERT2(vorout.numberofvpoints > 0, "Voronoi has no points");
  result.m_nodes.reserve(vorout.numberofvpoints);
  for (int i = 0; i < vorout.numberofvpoints; ++i) {
    RealPoint p(vorout.vpointlist[3*i],
                vorout.vpointlist[3*i+1],
                vorout.vpointlist[3*i+2]);
    IntPoint ip = result.m_Q.quantize(p);
    result.m_nodes.push_back(ip);
  }

  // Track infinite (boundary) nodes
  std::set<int> infNodeSet;

  // 2. Process Voronoi edges and handle infinite rays
  //    Build a map from Voronoi edge index to pair of node indices
  std::map<int, std::pair<int, int>> vedge2nodes;

  for (int i = 0; i < vorout.numberofvedges; ++i) {
    int v1 = vorout.vedgelist[i].v1;
    int v2 = vorout.vedgelist[i].v2;

    if (v2 == -1) {
      // Infinite edge - create a far point along the ray direction
      RealPoint origin(vorout.vpointlist[3*v1],
                       vorout.vpointlist[3*v1+1],
                       vorout.vpointlist[3*v1+2]);
      RealPoint dir(vorout.vedgelist[i].vnormal[0],
                    vorout.vedgelist[i].vnormal[1],
                    vorout.vedgelist[i].vnormal[2]);

      // Project to a far distance (use quantizer bounds)
      RealPoint far_dist = result.m_Q.maxCoord.template type_cast<double>();
      RealPoint inf_point = origin + dir * far_dist;

      // Quantize and add as new node
      IntPoint iinf = result.m_Q.quantize(inf_point);
      v2 = result.m_nodes.size();
      result.m_nodes.push_back(iinf);
      infNodeSet.insert(v2);
    }

    vedge2nodes[i] = {v1, v2};
  }

  // 3. Build faces from Voronoi facets
  //    Each Voronoi facet corresponds to a face between two Voronoi cells
  result.m_faces.reserve(vorout.numberofvfacets);
  result.m_cells.resize(vorout.numberofvcells);

  for (int ifacet = 0; ifacet < vorout.numberofvfacets; ++ifacet) {
    int c1 = vorout.vfacetlist[ifacet].c1;
    int c2 = vorout.vfacetlist[ifacet].c2;
    int nedges = vorout.vfacetlist[ifacet].elist[0];

    // Build ordered node list by following edge connectivity
    std::vector<int> face_nodes;

    if (nedges >= 3) {
      // Build adjacency map for this face
      std::map<int, std::vector<int>> node_neighbors;
      for (int j = 1; j <= nedges; ++j) {
        int edge_idx = vorout.vfacetlist[ifacet].elist[j];
        auto [n1, n2] = vedge2nodes[edge_idx];
        node_neighbors[n1].push_back(n2);
        node_neighbors[n2].push_back(n1);
      }

      // Start with first edge and follow the cycle
      if (!node_neighbors.empty()) {
        int start_node = node_neighbors.begin()->first;
        int current = start_node;
        int prev = -1;
        std::set<int> visited;

        face_nodes.push_back(current);
        visited.insert(current);

        // Follow edges around the cycle
        while (face_nodes.size() < nedges) {
          bool found = false;
          for (int next : node_neighbors[current]) {
            if (next != prev && visited.find(next) == visited.end()) {
              face_nodes.push_back(next);
              visited.insert(next);
              prev = current;
              current = next;
              found = true;
              break;
            }
          }
          if (!found) break;
        }
      }
    }

    // Need at least 3 nodes to form a face
    if (face_nodes.size() < 3) continue;

    int face_id = result.m_faces.size();
    result.m_faces.push_back(face_nodes);

    // Assign face to cells with proper orientation
    // c1 always gets the face, c2 gets the negated face (opposite orientation)
    if (c1 >= 0 && c1 < result.m_cells.size()) {
      result.m_cells[c1].push_back(face_id);
    }

    if (c2 >= 0 && c2 < result.m_cells.size()) {
      // Negative index indicates opposite orientation
      result.m_cells[c2].push_back(~face_id);
    }
  }

  // Alternative: use vcelllist if available
  // The vcelllist directly gives us which facets belong to each cell
  // if (vorout.vcelllist != nullptr) {
  //   result.m_cells.clear();
  //   result.m_cells.resize(vorout.numberofvcells);

  //   for (int icell = 0; icell < vorout.numberofvcells; ++icell) {
  //     if (vorout.vcelllist[icell] == nullptr) continue;

  //     int nfacets = vorout.vcelllist[icell][0];
  //     for (int j = 1; j <= nfacets; ++j) {
  //       int facet_idx = vorout.vcelllist[icell][j];

  //       // Determine orientation based on whether this cell is c1 or c2
  //       if (vorout.vfacetlist[facet_idx].c1 == icell) {
  //         result.m_cells[icell].push_back(facet_idx);
  //       } else if (vorout.vfacetlist[facet_idx].c2 == icell) {
  //         result.m_cells[icell].push_back(~facet_idx);
  //       }
  //     }
  //   }
  // }
}

void TetgenTessellator::setTetgenFacet(tetgenio::facet& f,
                                       const std::vector<int>& verts) const {
  tetgenio::init(&(f));
  f.numberofpolygons = 1;
  f.polygonlist = new tetgenio::polygon[1];
  f.numberofholes = 0;
  f.holelist = nullptr;

  tetgenio::polygon& p = f.polygonlist[0];
  tetgenio::init(&(p));
  p.numberofvertices = static_cast<int>(verts.size());
  p.vertexlist = new int[p.numberofvertices];
  for (int i = 0; i < p.numberofvertices; ++i) {
    p.vertexlist[i] = verts[i];
  }
}

// Create Tetgen class
tetgenio TetgenTessellator::createTetgenPoints(const QT& quant) const {
  // Since no PLC boundary is provided, we will create an external
  // boundary based on the maximum coordinates
  std::vector<IntPoint> guards = quant.generateGuards();
  PLC<3> plc;
  plc.facets = shapes::createCubeFaces();
  QuantPLC<3> qplc(plc, quant.m_Q, guards);
  tetgenio in = createTetgenPoints(qplc, quant);
  return in;
}

// Create Tetgen class for a given PLC
tetgenio TetgenTessellator::createTetgenPoints(const QPLC& qplc,
                                               const QT& quant) const {
  tetgenio in;
  auto points = quant.getRealPoints();
  auto bpoints = qplc.getRealPoints();
  const auto N_gen = points.size();
  const auto N_bound = bpoints.size();
  const auto N_total = N_gen + N_bound;
  in.pointlist = new REAL[N_total*3];
  for (int i = 0; i < N_gen; ++i) {
    auto rpoint = points[i].template type_cast<double>();
    for (int d = 0; d < 3; ++d) {
      in.pointlist[3*i+d] = rpoint[d];
    }
  }
  unsigned kk = N_gen;
  for (int i = 0; i < N_bound; ++i) {
    auto rpoint = bpoints[i].template type_cast<double>();
    for (int d = 0; d < 3; ++d) {
      in.pointlist[3*kk+d] = rpoint[d];
    }
    kk++;
  }
  in.numberofpoints = N_total;
  size_t N_facets = qplc.facets.size();
  for (const auto& hole : qplc.holes) {
    N_facets += hole.size();
  }
  in.numberoffacets = N_facets;
  in.facetlist = new tetgenio::facet[N_facets];
  in.facetmarkerlist = new int[N_facets];
  int f = 0;
  auto shiftedFacet = [&](const std::vector<int>& localFacet)
  {
    std::vector<int> globalFacet;
    globalFacet.reserve(localFacet.size());
    for (int idx : localFacet) {
      globalFacet.push_back(N_gen + idx);
    }
    return globalFacet;
  };
  // Outer boundary facets
  for (size_t i = 0; i < qplc.facets.size(); ++i) {
    std::vector<int> verts = shiftedFacet(qplc.facets[i]);
    setTetgenFacet(in.facetlist[f], verts);
    in.facetmarkerlist[f] = 1;
    ++f;
  }
  for (size_t h = 0; h < qplc.holes.size(); ++h) {
    for (size_t j = 0; j < qplc.holes[h].size(); ++j) {
      std::vector<int> verts = shiftedFacet(qplc.holes[h][j]);
      setTetgenFacet(in.facetlist[f], verts);
      in.facetmarkerlist[f] = 100 + static_cast<int>(h);
      ++f;
    }
  }
  // TODO: Set seeds
  return in;
}

}
