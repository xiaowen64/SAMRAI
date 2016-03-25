//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/skeleton/patch_geom/SkeletonPatchGeometry.C $
// Package:     SAMRAI geometry package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
template<int DIM> void SkeletonPatchGeometry<DIM>::printClassData(std::ostream& os) const
{
   os << "Printing SkeletonPatchGeometry data: this = "
      << (SkeletonPatchGeometry<DIM>*)this << std::endl;
 
   hier::PatchGeometry<DIM>::printClassData(os);
}

}
}
#endif
