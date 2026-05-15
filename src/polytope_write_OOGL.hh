//------------------------------------------------------------------------------
// Output OOGL (Object Oriented Graphics Library) files from polytope types.
// Mostly useful visualizing with Geomview.
//------------------------------------------------------------------------------
#ifndef __Polytope_writeOOGL__
#define __Polytope_writeOOGL__

#include <vector>
#include <set>
#include <string>
#include <fstream>

#include "PLC.hh"
#include "polytope_internal.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Write a 2D OFF polylist file.
//------------------------------------------------------------------------------
template<typename RealType>
void writePLCtoOFF(const PLC<2>& plc,
                   const std::vector<RealType>& coords,
                   const std::string filename) {
  
  typedef std::pair<int, int> EdgeHash;
  POLY_ASSERT(coords.size() % 2 == 0);

  // Ope the file for output.
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out);

  // Count the numbers of elements we have.
  const unsigned npoints = coords.size()/2, nfacets = plc.facets.size();
  std::set<EdgeHash> edges;
  for (unsigned i = 0; i != plc.facets.size(); ++i) {
    POLY_ASSERT(plc.facets[i].size() == 2);
    edges.insert(internal::hashEdge(plc.facets[i][0], plc.facets[i][1]));
  }
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    for (unsigned i = 0; i != plc.holes[ihole].size(); ++i) {
      POLY_ASSERT(plc.holes[ihole][i].size() == 2);
      edges.insert(internal::hashEdge(plc.holes[ihole][i][0], plc.holes[ihole][i][1]));
    }
  }
  const unsigned nedges = edges.size();

  // Write the header.
  // First we write the OFF string to indicate the format, and then
  // the number of point, facets, and edges.
  outfile << "OFF" << std::endl;
  outfile << npoints << " " << nfacets << " " << nedges << std::endl;

  // Write the vertex coordinates.
  for (unsigned i = 0; i != npoints; ++i) {
    outfile << coords[2*i] << " " << coords[2*i+1] << std::endl;
  }

  // Write the outer facets.
  for (unsigned i = 0; i != plc.facets.size(); ++i) {
    const unsigned n = plc.facets[i].size();
    POLY_ASSERT(n == 2);
    outfile << n;
    for (unsigned j = 0; j != n; ++j) {
      outfile << " " << plc.facets[i][j];
    }
    outfile << std::endl;
  }

  // Write the holes.
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    for (unsigned i = 0; i != plc.holes[ihole].size(); ++i) {
      const unsigned n = plc.holes[ihole][i].size();
      POLY_ASSERT(n == 2);
      outfile << n;
      for (unsigned j = 0; j != n; ++j) {
        outfile << " " << plc.holes[ihole][i][j];
      }
      outfile << std::endl;
    }
  }

  // Close the files.
  outfile.close();
}

//------------------------------------------------------------------------------
// Write a 3D OFF polylist file.
//------------------------------------------------------------------------------
template<typename RealType>
void writePLCtoOFF(const PLC<3>& plc,
                   const std::vector<RealType>& coords,
                   const std::string filename) {

  typedef std::pair<int, int> EdgeHash;
  POLY_ASSERT(coords.size() % 3 == 0);

  // Ope the file for output.
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out);

  // Count the numbers of elements we have.
  const unsigned npoints = coords.size()/3, nfacets = plc.facets.size();
  std::set<EdgeHash> edges;
  // for (unsigned i = 0; i != plc.facets.size(); ++i) {
  //   const unsigned n = plc.facets[i].size();
  //   for (unsigned j = 0; j != n; ++j) {
  //     const unsigned k = (j + 1) % n;
  //     edges.insert(internal::hashEdge(plc.facets[i][j], plc.facets[i][k]));
  //   }
  // }
  // for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
  //   for (unsigned i = 0; i != plc.holes[ihole].size(); ++i) {
  //     const unsigned n = plc.holes[ihole][i].size();
  //     for (unsigned j = 0; j != n; ++j) {
  //       const unsigned k = (j + 1) % n;
  //       edges.insert(internal::hashEdge(plc.holes[ihole][i][j], plc.holes[ihole][i][k]));
  //     }
  //   }
  // }
  const unsigned nedges = edges.size();

  // Write the header.
  // First we write the OFF string to indicate the format, and then
  // the number of point, facets, and edges.
  outfile << "OFF" << std::endl;
  outfile << npoints << " " << nfacets << " " << nedges << std::endl;

  // Write the vertex coordinates.
  for (unsigned i = 0; i != npoints; ++i) {
    outfile << coords[3*i] << " " << coords[3*i+1] << " " << coords[3*i+2] << std::endl;
  }

  // Write the outer facets.
  for (unsigned i = 0; i != plc.facets.size(); ++i) {
    const unsigned n = plc.facets[i].size();
    outfile << n;
    for (unsigned j = 0; j != n; ++j) {
      outfile << " " << plc.facets[i][j];
    }
    outfile << std::endl;
  }

  // Write the holes.
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    for (unsigned i = 0; i != plc.holes[ihole].size(); ++i) {
      const unsigned n = plc.holes[ihole][i].size();
      outfile << n;
      for (unsigned j = 0; j != n; ++j) {
        outfile << " " << plc.holes[ihole][i][j];
      }
      outfile << std::endl;
    }
  }

  // Close the files.
  outfile.close();
}

}

#endif
