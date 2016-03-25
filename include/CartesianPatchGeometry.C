//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/cartesian/patch_geom/CartesianPatchGeometry.C $
// Package:     SAMRAI geometry package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Simple Cartesian grid geometry for an AMR hierarchy.
//

#ifndef included_geom_CartesianPatchGeometry_C
#define included_geom_CartesianPatchGeometry_C

#include "CartesianPatchGeometry.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
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
   TBOX_ASSERT(!(dx == (double*)NULL));
   TBOX_ASSERT(!(x_lo == (double*)NULL));
   TBOX_ASSERT(!(x_up == (double*)NULL));
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
template<int DIM> void CartesianPatchGeometry<DIM>::printClassData(std::ostream& os) const
{
   os << "Printing CartesianPatchGeometry data: this = "
      << (CartesianPatchGeometry*)this << std::endl;
   os << "x_lo = ";
   for (int id1 = 0; id1 < DIM; id1++) {
      os << d_x_lo[id1] << "   ";
   }
   os << std::endl;
   os << "x_up = ";
   for (int id2 = 0; id2 < DIM; id2++) {
      os << d_x_up[id2] << "   ";
   }
   os << std::endl;
   os << "dx = ";
   for (int id3 = 0; id3 < DIM; id3++) {
      os << d_dx[id3] << "   ";
   }
   os << std::endl;
 
   hier::PatchGeometry<DIM>::printClassData(os);
}

}
}
#endif
