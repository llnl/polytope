#ifndef POLYTOPE_PLC_HH
#define POLYTOPE_PLC_HH

#include <cstddef>
#include <vector>
#include <iostream>

namespace polytope {

//! \class PLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
//! Graph (PSLG) in 2D.
template<int Dimension, typename IndexType>
class PLCBase {
public:

  //! This array defines the topology of the facets of the 
  //! piecewise linear complex in terms of connections to generating points. 
  //! A facet has an arbitrary number of points in 3D and 2 points in 2D. 
  //! facets[i] gives the index of the generating point of the ith 
  //! facet.
  std::vector<IndexType > facets;

  //! This array defines the topology of the inner facets
  //! or holes in the geometry.  The outer-most dimension is the number of 
  //! holes, and the remaining are facets using the same convention as the
  //! the "facets" member.  In other words, holes[k][i] is the
  //! generating point of the ith facet of the kth hole.
  std::vector<std::vector<IndexType > > holes;
  
  //! Clears facets and holes to empty the PLC
  void clear() {
    facets.clear();
    holes.clear();
  }

  //! Returns true if this PLC is empty, false otherwise.
  bool empty() const {
    return (facets.empty() and holes.empty());
  }

  virtual bool valid() const = 0;
};

template<int Dimension>
class PLC : public PLCBase<Dimension, std::vector<int>> {
public:
  using PLCBase<Dimension, std::vector<int>>::facets;
  using PLCBase<Dimension, std::vector<int>>::holes;
  //! Returns true if this PLC is valid (at first glance), false if it 
  //! is obviously invalid. This is not a rigorous check!
  bool valid() const {
    if (Dimension == 2) {
      // In 2D all facets must have at least 2 points.
      for (std::size_t f = 0; f < facets.size(); ++f) {
        if (facets[f].size() != 2) {
          return false;
        }
      }
      for (std::size_t h = 0; h < holes.size(); ++h) {
        for (std::size_t f = 0; f < holes[h].size(); ++f) {
          if (holes[h][f].size() != 2) {
            return false;
          }
        }
      }
    } else if (Dimension == 3) {
      // In 3D all facets must have at least 3 points.
      for (std::size_t f = 0; f < facets.size(); ++f) {
        if (facets[f].size() < 3) {
          return false;
        }
      }
      for (std::size_t h = 0; h < holes.size(); ++h) {
        for (std::size_t f = 0; f < holes[h].size(); ++f) {
          if (holes[h][f].size() < 3) {
            return false;
          }
        }
      }
    }
    return true;
  }

  //! output operator.
  friend std::ostream& operator<<(std::ostream& s, const PLC& plc) {
    s << "PLC (" << Dimension << "D):" << std::endl;
    s << plc.facets.size() << " facets:" << std::endl;
    for (std::size_t f = 0; f < plc.facets.size(); ++f) {
      s << " " << f << ": (";
      for (std::size_t p = 0; p < plc.facets[f].size(); ++p) {
        if (p < plc.facets[f].size()-1) {
          s << plc.facets[f][p] << ", ";
        } else {
          s << plc.facets[f][p];
        }
      }
      s << ")" << std::endl;
    }
    s << std::endl;
    s << plc.holes.size() << " holes:" << std::endl;
    for (std::size_t h = 0; h < plc.holes.size(); ++h) {
      s << "Hole #" << h << std::endl;
      for (std::size_t f = 0; f < plc.holes[h].size(); ++f) {
        s << "    " << f << ": (";
        for (std::size_t p = 0; p < plc.holes[h][f].size(); ++p) {
          if (p < plc.holes[h][f].size()-1) {
            s << plc.holes[h][f][p] << ", ";
          } else {
            s << plc.holes[h][f][p];
          }
        }
        s << ")" << std::endl;
      }
    }
    return s;
  }

};

}

#endif
