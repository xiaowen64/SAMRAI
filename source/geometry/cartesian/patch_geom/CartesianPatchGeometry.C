//
// File:        CartesianPatchGeometry.C
// Package:     SAMRAI geometry package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Simple Cartesian grid geometry for an AMR hierarchy.
//

#ifndef included_geom_CartesianPatchGeometry_C
#define included_geom_CartesianPatchGeometry_C

#include "CartesianPatchGeometry.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#include "CartesianPatchGeometry.I"
#endif
namespace SAMRAI {
    namespace geom {

/*
*************************************************************************
*                                                                       *
* Constructor for CartesianPatchGeometry allocates and sets        *
* patch coordinate system information.                                  *
*                                                                       *
*************************************************************************
*/
template<int DIM>  CartesianPatchGeometry<DIM>::CartesianPatchGeometry(
   const hier::IntVector<DIM>& ratio_to_level_zero,
   const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
   const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry,
   const double* dx,
   const double* x_lo,
   const double* x_up)
: hier::PatchGeometry<DIM>(ratio_to_level_zero,
                           touches_regular_bdry, touches_periodic_bdry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(dx == (double*)NULL));
   assert(!(x_lo == (double*)NULL));
   assert(!(x_up == (double*)NULL));
#endif
   for (int id = 0; id < DIM; id++) {
      d_dx[id]   = dx[id];
      d_x_lo[id] = x_lo[id];
      d_x_up[id] = x_up[id];
   }
}


/*
*************************************************************************
*                                                                       *
* Destructor for CartesianPatchGeometry deallocates dx array.      *
*                                                                       *
*************************************************************************
*/
template<int DIM>  CartesianPatchGeometry<DIM>::~CartesianPatchGeometry()
{
}


/*
*************************************************************************
*                                                                       *
* Print CartesianPatchGeometry class data.                         *
*                                                                       *
*************************************************************************
*/
template<int DIM> void CartesianPatchGeometry<DIM>::printClassData(ostream& os) const
{
   os << "Printing CartesianPatchGeometry data: this = "
      << (CartesianPatchGeometry*)this << endl;
   os << "x_lo = ";
   for (int id1 = 0; id1 < DIM; id1++) {
      os << d_x_lo[id1] << "   ";
   }
   os << endl;
   os << "x_up = ";
   for (int id2 = 0; id2 < DIM; id2++) {
      os << d_x_up[id2] << "   ";
   }
   os << endl;
   os << "dx = ";
   for (int id3 = 0; id3 < DIM; id3++) {
      os << d_dx[id3] << "   ";
   }
   os << endl;
 
   hier::PatchGeometry<DIM>::printClassData(os);
}

}
}
#endif
