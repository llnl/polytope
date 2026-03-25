//------------------------------------------------------------------------------
// A bunch of useful canned methods to generate PLC's of various geometries.
//------------------------------------------------------------------------------
#ifndef __Polytope_canned_geometries__
#define __Polytope_canned_geometries__

#include "ReducedPLC.hh"
#include "Point.hh"
#include "convexHull_3d.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Create a PLC representation of a circle.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, RealType>
plc_circle(const Point2<RealType>& center,
           const RealType radius,
           const unsigned ntheta) {
  POLY_ASSERT(ntheta > 2);

  // Prepare the result.
  ReducedPLC<2, RealType> result;

  // Put a bunch of points on the sphere.
  const double dtheta = 2.0*M_PI/ntheta;
  result.facets.resize(ntheta);
  for (unsigned itheta = 0; itheta != ntheta; ++itheta) {
    result.facets[itheta].push_back(itheta);
    result.facets[itheta].push_back((itheta + 1) % ntheta);
    const double theta = (itheta + 0.5)*dtheta;
    result.points.push_back(center.x + radius*cos(theta));
    result.points.push_back(center.y + radius*sin(theta));
  }
  return result;
}

//------------------------------------------------------------------------------
// Create a PLC representation of a sphere.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
plc_sphere(const Point3<RealType>& center,
           const RealType radius,
           const unsigned nphi) {
  POLY_ASSERT(nphi > 0);

  // Prepare the result.
  ReducedPLC<3, RealType> result;

  // Put a bunch of points on the sphere.
  const double dphi = M_PI/nphi;
  for (unsigned iphi = 0; iphi != nphi; ++iphi) {
    const double phi = (iphi + 0.5)*dphi;
    const double dl = radius*dphi;
    const double rxy = radius*sin(phi);
    const double circ = 2.0*M_PI*rxy;
    const unsigned ntheta = std::max(1U, unsigned(circ/dl + 0.5));
    const double dtheta = 2.0*M_PI/ntheta;
    for (unsigned itheta = 0; itheta != ntheta; ++itheta) {
      const double theta = (itheta + 0.5)*dtheta;
      result.points.push_back(center.x + radius*cos(theta)*sin(phi));
      result.points.push_back(center.y + radius*sin(theta)*sin(phi));
      result.points.push_back(center.z + radius*cos(phi));
    }
  }

  // Build the convex hull of our points.
  const RealType low[3] = {center.x - 1.1*radius,
                           center.y - 1.1*radius,
                           center.z - 1.1*radius};
  const RealType dx = radius/(1 << 21);
  const PLC<3, RealType> hull = convexHull_3d(result.points, low, dx);

  // Put that topology in our result and we're done.
  result.facets = hull.facets;
  return result;
}

//------------------------------------------------------------------------------
// Create a PLC representation of a cylinder.
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType>
plc_cylinder(const Point3<RealType>& center,
             const RealType radius,
             const RealType length,
             const unsigned ncirc) {
  POLY_ASSERT(ncirc >= 3);

  // Prepare the result.
  ReducedPLC<3, RealType> result;

  // Put points on the two end circles.
  const double dtheta = 2.0*M_PI/ncirc;
  for (unsigned itheta = 0; itheta != ncirc; ++itheta) {
    const double theta = (itheta + 0.5)*dtheta;
    result.points.push_back(center.x + radius*cos(theta));
    result.points.push_back(center.y + radius*sin(theta));
    result.points.push_back(center.z + 0.5*length);
    result.points.push_back(center.x + radius*cos(theta));
    result.points.push_back(center.y + radius*sin(theta));
    result.points.push_back(center.z - 0.5*length);
  }

  // Build the convex hull of our points.
  const RealType low[3] = {center.x - 1.1*radius,
                           center.y - 1.1*radius,
                           center.z - 0.6*length};
  const RealType dx = length/(1 << 21);
  const PLC<3, RealType> hull = convexHull_3d(result.points, low, dx);

  // Put that topology in our result and we're done.
  result.facets = hull.facets;
  return result;
}

//------------------------------------------------------------------------------
// Return a PLC box.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
ReducedPLC<Dimension, RealType>
plc_box(const RealType* low, const RealType* high) {
  ReducedPLC<Dimension, RealType> box;

  if (Dimension == 2) {
    // 2D
    // Add the new generators to points.
    box.points.push_back(low[0]);
    box.points.push_back(low[1]);

    box.points.push_back(high[0]);
    box.points.push_back(low[1]);

    box.points.push_back(high[0]);
    box.points.push_back(high[1]);

    box.points.push_back(low[0]);
    box.points.push_back(high[1]);

    // Construct the box.
    box.facets.resize(4);

    // -y face.
    box.facets[0].resize(2);
    box.facets[0][0] = 0;
    box.facets[0][1] = 1;

    // +x face.
    box.facets[1].resize(2);
    box.facets[1][0] = 1;
    box.facets[1][1] = 2;

    // +y face.
    box.facets[2].resize(2);
    box.facets[2][0] = 2;
    box.facets[2][1] = 3;

    // -x face.
    box.facets[3].resize(2);
    box.facets[3][0] = 3;
    box.facets[3][1] = 0;

    return box;

  } else {

    // 3D
    const RealType x1 = low[0], x2 = high[0],
      y1 = low[1], y2 = high[1],
      z1 = low[2], z2 = high[2];

    // Create the piecewise linear complex representing the box. Note that 
    // the box consists of facets that are defined by their connections to 
    // generating points.
    // Should look like the following:
    //
    //        6--------7            y
    //       /        /|            |
    //      /        / |            |
    //     2--------3  |             ------x
    //     |  .     |  |           /
    //     |  4.....|..5          z
    //     | .      | / 
    //     |.       |/
    //     0--------1             
    //
    // Create the vertices for our bounding surface.
    box.points.resize(3*8);
    box.points[3*0+0] = x1; box.points[3*0+1] = y1; box.points[3*0+2] = z2;
    box.points[3*1+0] = x2; box.points[3*1+1] = y1; box.points[3*1+2] = z2;
    box.points[3*2+0] = x1; box.points[3*2+1] = y2; box.points[3*2+2] = z2;
    box.points[3*3+0] = x2; box.points[3*3+1] = y2; box.points[3*3+2] = z2;
    box.points[3*4+0] = x1; box.points[3*4+1] = y1; box.points[3*4+2] = z1;
    box.points[3*5+0] = x2; box.points[3*5+1] = y1; box.points[3*5+2] = z1;
    box.points[3*6+0] = x1; box.points[3*6+1] = y2; box.points[3*6+2] = z1;
    box.points[3*7+0] = x2; box.points[3*7+1] = y2; box.points[3*7+2] = z1;

    // 6 facets
    box.facets.resize(6);

    // facet 0 -- bottom face.
    box.facets[0].resize(4);
    box.facets[0][0] = 0;
    box.facets[0][1] = 4;
    box.facets[0][2] = 5;
    box.facets[0][3] = 1;

    // facet 1 -- top face.
    box.facets[1].resize(4);
    box.facets[1][0] = 2;
    box.facets[1][1] = 3;
    box.facets[1][2] = 7;
    box.facets[1][3] = 6;

    // facet 2 -- left face.
    box.facets[2].resize(4);
    box.facets[2][0] = 0;
    box.facets[2][1] = 2;
    box.facets[2][2] = 6;
    box.facets[2][3] = 4;

    // facet 3 -- right face.
    box.facets[3].resize(4);
    box.facets[3][0] = 1;
    box.facets[3][1] = 5;
    box.facets[3][2] = 7;
    box.facets[3][3] = 3;

    // facet 4 -- front face.
    box.facets[4].resize(4);
    box.facets[4][0] = 0;
    box.facets[4][1] = 1;
    box.facets[4][2] = 3;
    box.facets[4][3] = 2;

    // facet 5 -- back face.
    box.facets[5].resize(4);
    box.facets[5][0] = 5;
    box.facets[5][1] = 4;
    box.facets[5][2] = 6;
    box.facets[5][3] = 7;

    // That's it.
    return box;
  }
}

}

#endif
