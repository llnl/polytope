#ifndef POLYTOPE_TETGENTESSELLATOR_HH
#define POLYTOPE_TETGENTESSELLATOR_HH

#include "QuantPLC.hh"
#include "QuantTessellation.hh"
#include "Tessellator.hh"
#include "Shapes.hh"

#define TETLIBRARY
#include "tetgen.h"

namespace polytope {

class TetgenTessellator : public Tessellator<3, double> {
public:
  using RealType = double;
  using QT = QuantTessellation<3>;
  using QPLC = QuantPLC<3>;
  using IntType = QT::IntType;
  using IntPoint = QT::IntPoint;
  using RealPoint = Point<3, double>;

  TetgenTessellator() = default;
  TetgenTessellator(const Quantizer<3>& Q) :
    Tessellator(Q) {}

  virtual ~TetgenTessellator() = default;

  virtual void tessellateQuantized(QT& result) const;
  virtual void tessellateQuantized(const QPLC& qplc, QT& result) const;

  std::string name() const { return "TetgenTessellator"; }
protected:
  void setTetgenFacet(tetgenio::facet& f, const std::vector<int>& verts) const;

  tetgenio createTetgenPoints(const QT& quant) const;
  tetgenio createTetgenPoints(const QPLC& qplc, const QT& quant) const;

  void convertVoronoiToQuantTessellation(const tetgenio& vorout, QT& result) const;
};

}
#endif
