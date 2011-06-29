/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Skeleton grid geometry for an AMR hierarchy. 
 *
 ************************************************************************/

#ifndef included_geom_SkeletonPatchGeometry_C
#define included_geom_SkeletonPatchGeometry_C

#include "SAMRAI/geom/SkeletonPatchGeometry.h"

namespace SAMRAI {
namespace geom {

// using namespace std;

/*
 *************************************************************************
 *                                                                       *
 * Constructor for SkeletonPatchGeometry.                           *
 * variable.                                                             *
 *                                                                       *
 *************************************************************************
 */
SkeletonPatchGeometry::SkeletonPatchGeometry(
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const TwoDimBool& touches_periodic_bdry):
   hier::PatchGeometry(ratio_to_level_zero,
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
SkeletonPatchGeometry::~SkeletonPatchGeometry()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Print SkeletonPatchGeometry class data.                          *
 *                                                                       *
 *************************************************************************
 */
void SkeletonPatchGeometry::printClassData(
   std::ostream& os) const
{
   os << "Printing SkeletonPatchGeometry data: this = "
      << (SkeletonPatchGeometry *)this << std::endl;

   hier::PatchGeometry::printClassData(os);
}

}
}
#endif
