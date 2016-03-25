//
// File:        EmbeddedBoundaryShapeSphere.C
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 605 $
// Modified:    $Date: 2005-09-09 15:39:55 -0700 (Fri, 09 Sep 2005) $
// Description: Sphere embedded boundary shape
//              
// 

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace appu {

template<int DIM> inline
bool EmbeddedBoundaryShapeSphere<DIM>::isInside(const double* xyz) const
{   
   /*
    * Compare the radius of the point from center to the radius
    * of the ball from center.
    */
    double rad_squared = 0.;
    for (int i=0; i < DIM; i++) {
       rad_squared += (xyz[i] - d_center[i])*(xyz[i] - d_center[i]);
    }
    return((d_radius*d_radius - rad_squared) > 0.);
}

/*
 * Array based version of isInside - set the node-centered "inout"
 * array on a mesh defined by the provided nx, dx, and origin.
 */
template<int DIM> inline
void EmbeddedBoundaryShapeSphere<DIM>::isInside(const int* nx,
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
   double xyz[DIM];
   
   /*
    * Compute the node centered inout points
    */
   int khi = 1;
   if (DIM == 3) {
      khi = nx[2];
   }
   for (int kk = 0; kk < khi; kk++) {
      kstr = kk*nx[1]*nx[0];
      if (DIM == 3) {
         xyz[DIM-1] = origin[DIM-1] + dx[DIM-1]*(double)kk;
      }
      for (int jj = 0; jj < nx[1]; jj++) {
         jstr = jj*nx[0];
         xyz[1] = origin[1] + dx[1]*(double)jj;
         for (int ii = 0; ii < nx[0]; ii++) {
            istr = ii;
            xyz[0] = origin[0] + dx[0]*(double)ii;

            double rad_squared = 0.;
            for (int dim = 0; dim < DIM; dim++) {
               rad_squared += 
                  (xyz[dim] - d_center[dim])*(xyz[dim] - d_center[dim]);
            }

            ijk = kstr + jstr + istr;
            rad_diff = d_radius*d_radius - rad_squared;

            if (tbox::Utilities::deq(rad_diff,0.)) {
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

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
