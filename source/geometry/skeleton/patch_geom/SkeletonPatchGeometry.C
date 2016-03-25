//
// File:        SkeletonPatchGeometry.C
// Package:     SAMRAI geometry package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 332 $
// Modified:    $Date: 2005-05-03 11:39:46 -0700 (Tue, 03 May 2005) $
// Description: Skeleton grid geometry for an AMR hierarchy.
//

#ifndef included_geom_SkeletonPatchGeometry_C
#define included_geom_SkeletonPatchGeometry_C

#include "SkeletonPatchGeometry.h"

namespace SAMRAI {
    namespace geom {

/*
*************************************************************************
*                                                                       *
* Constructor for SkeletonPatchGeometry.                           *
* variable.                                                             *
*                                                                       *
*************************************************************************
*/
template<int DIM>  SkeletonPatchGeometry<DIM>::SkeletonPatchGeometry(
   const hier::IntVector<DIM>& ratio_to_level_zero,
   const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
   const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry)
: hier::PatchGeometry<DIM>(ratio_to_level_zero,
                           touches_regular_bdry,
                           touches_periodic_bdry)
{
}


/*
*************************************************************************
*                                                                       *
* Destructor for SkeletonPatchGeometry.                            *
*                                                                       *
*************************************************************************
*/
template<int DIM>  SkeletonPatchGeometry<DIM>::~SkeletonPatchGeometry()
{
}


/*
*************************************************************************
*                                                                       *
* Print SkeletonPatchGeometry class data.                          *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SkeletonPatchGeometry<DIM>::printClassData(ostream& os) const
{
   os << "Printing SkeletonPatchGeometry data: this = "
      << (SkeletonPatchGeometry<DIM>*)this << endl;
 
   hier::PatchGeometry<DIM>::printClassData(os);
}

}
}
#endif
