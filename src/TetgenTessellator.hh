#ifndef __Polytope_TetgenTessellator__
#define __Polytope_TetgenTessellator__

#include "QuantPLC.hh"
#include "QuantTessellation.hh"
#include "Tessellator.hh"

// Forward declaration
class tetgenio;

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
  virtual ~TetgenTessellator() = default;

  virtual void tessellateQuantizedImpl(QT& result) override;

  std::string name() const { return "TetgenTessellator"; }
protected:
  tetgenio createTetgenPoints(const QT& result) const;
};

}
#endif
