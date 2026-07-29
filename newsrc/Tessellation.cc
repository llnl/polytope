
#include "Tessellation.hh"
#include "GeomUtils.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Compute the centroid of the mesh cell.
//------------------------------------------------------------------------------
// template<int Dimension, typename RealType>
// void
// Tessellation<Dimension, RealType>::
// computeCellCentroid(const unsigned ci,
//                     RealType* ccent) const {
//   POLY_ASSERT(ci < cells.size());
//   const unsigned nf = cells[ci].size();
//   unsigned i, j, k;
//   for (j = 0; j != Dimension; ++j) ccent[j] = 0.0;
//   std::vector<unsigned> uniqueNodes;
//   for (k = 0; k != nf; ++k) {
//     const std::vector<unsigned>& faceNodes = faces[internal::positiveID(cells[ci][k])];
//     std::copy(faceNodes.begin(), faceNodes.end(), std::back_inserter(uniqueNodes));
//   }

//   // Now reduce to the unique vertices.
//   std::sort(uniqueNodes.begin(), uniqueNodes.end());
//   std::vector<unsigned>::iterator endItr = std::unique(uniqueNodes.begin(), uniqueNodes.end());
//   for (std::vector<unsigned>::iterator itr = uniqueNodes.begin(); itr != endItr; ++itr) {
//     i = *itr;
//     POLY_ASSERT2(i < nodes.size()/Dimension, i << " " << nodes.size()/Dimension);
//     for (j = 0; j != Dimension; ++j) ccent[j] += nodes[Dimension*i + j];
//   }
//   const unsigned n = std::distance(uniqueNodes.begin(), endItr);
//   POLY_ASSERT(n > 0);
//   for (j = 0; j != Dimension; ++j) ccent[j] /= n;
// }

//------------------------------------------------------------------------------
// Compute the cell centroid and signed area.
// Taken from http://www.wikipedia.org/wiki/Centroid
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
template<int D>
std::enable_if_t<D == 2, void>
Tessellation<Dimension, RealType>::
computeCellCentroidAndSignedArea(const unsigned ci,
				 const RealType& tol,
				 RealType* ccent,
				 RealType& area) const {
  POLY_ASSERT(ci < cells.size());
  unsigned iface, n0, n1;
  RealType d, x0, x1, y0, y1;
  ccent[0] = 0.0; ccent[1] = 0.0;  area = 0.0;
  for (std::vector<int>::const_iterator itr = cells[ci].begin();
       itr != cells[ci].end(); ++itr) {
    iface = (*itr < 0) ? ~(*itr) : *itr;
    POLY_ASSERT(iface < faces.size());
    POLY_ASSERT(faces[iface].size() == 2);
    n0 = (*itr < 0) ? faces[iface][1] : faces[iface][0];
    n1 = (*itr < 0) ? faces[iface][0] : faces[iface][1];
    POLY_ASSERT(n0 < nodes.size() and n1 < nodes.size());
    POLY_ASSERT(n0 != n1);
    x0 = nodes[n0].x;  y0 = nodes[n0].y;
    x1 = nodes[n1].x;  y1 = nodes[n1].y;
    d = x0*y1 - y0*x1;
    area     += d;
    ccent[0] += d*(x0+x1);
    ccent[1] += d*(y0+y1);
  }
  POLY_ASSERT(std::abs(area) > tol);
  area     /= 2.0;
  ccent[0] /= (6*area);
  ccent[1] /= (6*area);
}

//------------------------------------------------------------------------------
// Compute the centroid and unit normal of a Tessellation face.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
template<int D>
std::enable_if_t<D == 3, void>
Tessellation<Dimension, RealType>::
computeFaceCentroidAndNormal(const unsigned fi,
                             RealType* fcent,
                             RealType* fhat) const {
  POLY_ASSERT(fi < faces.size());
  const unsigned n = faces[fi].size();
  POLY_ASSERT(n >= 3);
  unsigned i, ni;
  std::vector<unsigned> verts;
  const double degeneracy = 1.0e-10;

  // Compute the centroid, and look for three vertices in the face that are
  // not collinear.
  fcent[0] = 0.0; fcent[1] = 0.0; fcent[2] = 0.0;
  for (i = 0; i != n; ++i) {
    ni = faces[fi][i];
    POLY_ASSERT(ni < nodes.size());
    fcent[0] += nodes[ni].x;
    fcent[1] += nodes[ni].y;
    fcent[2] += nodes[ni].z;
    if (verts.size() < 2 or
        (verts.size() == 2 and not collinear<3, RealType>(&nodes[verts[0]][0],
                                                          &nodes[verts[1]][0],
                                                          &nodes[ni],
                                                          degeneracy))) verts.push_back(ni);
  }
  POLY_ASSERT(n > 0);
  fcent[0] /= n; fcent[1] /= n; fcent[2] /= n;

  // Now we can compute the unit normal.
  POLY_ASSERT2(verts.size() == 3, verts.size());
  RealType ab[3], ac[3];
  ab[0] = nodes[verts[1]].x - nodes[verts[0]].x;
  ab[1] = nodes[verts[1]].y - nodes[verts[0]].y;
  ab[2] = nodes[verts[1]].z - nodes[verts[0]].z;
  ac[0] = nodes[verts[2]].x - nodes[verts[0]].x;
  ac[1] = nodes[verts[2]].y - nodes[verts[0]].y;
  ac[2] = nodes[verts[2]].z - nodes[verts[0]].z;
  cross<3, RealType>(ab, ac, fhat);
  UnitVector<3, RealType>(fhat);
}

}
