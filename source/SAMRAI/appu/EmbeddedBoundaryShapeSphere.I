/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Sphere embedded boundary shape 
 *
 ************************************************************************/

#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace appu {

SAMRAI_INLINE_KEYWORD
bool EmbeddedBoundaryShapeSphere::isInside(
   const double* xyz) const
{
   /*
    * Compare the radius of the point from center to the radius
    * of the ball from center.
    */
   double rad_squared = 0.;
   for (int i = 0; i < d_dim.getValue(); i++) {
      rad_squared += (xyz[i] - d_center[i]) * (xyz[i] - d_center[i]);
   }
   return (d_radius * d_radius - rad_squared) > 0.;
}

/*
 * Array based version of isInside - set the node-centered "inout"
 * array on a mesh defined by the provided nx, dx, and origin.
 */
SAMRAI_INLINE_KEYWORD
void EmbeddedBoundaryShapeSphere::isInside(
   const int* nx,
   const double* dx,
   const double* origin,
   int* inout) const
{
   /*
    * Compare the radius of the point from center to the radius
    * of the ball from center.
    */
   int kstr = 0;
   int jstr = 0;
   int istr = 0;
   int ijk;
   double rad_diff;
   double xyz[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Compute the node centered inout points
    */
   int khi = 1;
   if (d_dim == tbox::Dimension(3)) {
      khi = nx[2];
   }
   for (int kk = 0; kk < khi; kk++) {
      kstr = kk * nx[1] * nx[0];
      if (d_dim == tbox::Dimension(3)) {
         xyz[d_dim.getValue()
             - 1] = origin[d_dim.getValue() - 1] + dx[d_dim.getValue() - 1] * (double)kk;
      }
      for (int jj = 0; jj < nx[1]; jj++) {
         jstr = jj * nx[0];
         xyz[1] = origin[1] + dx[1] * (double)jj;
         for (int ii = 0; ii < nx[0]; ii++) {
            istr = ii;
            xyz[0] = origin[0] + dx[0] * (double)ii;

            double rad_squared = 0.;
            for (int d = 0; d < d_dim.getValue(); ++d) {
               rad_squared +=
                  (xyz[d] - d_center[d]) * (xyz[d] - d_center[d]);
            }

            ijk = kstr + jstr + istr;
            rad_diff = d_radius * d_radius - rad_squared;

            if (tbox::MathUtilities<double>::equalEps(rad_diff, 0.)) {
               inout[ijk] = EmbeddedBoundaryDefines::ONBOUNDARY;
            } else if (rad_diff > 0.) {
               inout[ijk] = EmbeddedBoundaryDefines::INSIDE;
            } else {
               inout[ijk] = EmbeddedBoundaryDefines::OUTSIDE;
            }
         }
      }
   }

}

SAMRAI_INLINE_KEYWORD
const std::string& EmbeddedBoundaryShapeSphere::getObjectName() const
{
   return d_object_name;
}

}
}
