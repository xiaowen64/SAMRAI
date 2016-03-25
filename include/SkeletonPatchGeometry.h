//
// File:	SkeletonPatchGeometry.h
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Simple Skeleton grid geometry for an AMR hierarchy.
//

#ifndef included_geom_SkeletonPatchGeometry
#define included_geom_SkeletonPatchGeometry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchGeometry
#include "PatchGeometry.h"
#endif

namespace SAMRAI {
    namespace geom {

/**
 * Class SkeletonPatchGeometry implements geometry management 
 * for a single patch in an AMR hierarchy with no information about
 * the physical characteristics of the problem.  This is intended for
 * use in an application that will manage the physical geometry in
 * user-defined code.
 * 
 * The grid data is set by SkeletonGridGeometry class.  This patch
 * geometry class is derived from hier::PatchGeometry<DIM> base class.
 *
 * @see hier::BoundaryBox
 * @see hier::PatchGeometry
 * @see geom::SkeletonGridGeometry
 */

template<int DIM> class SkeletonPatchGeometry 
: public hier::PatchGeometry<DIM>
{
public:

   /**
    * Constructor for SkeletonPatchGeometry class.  It simply passes 
    * patch boundary information and the ratio to the coarsest level to
    * hier::PatchGeometry constructor.
    */
   SkeletonPatchGeometry<DIM>(
      const hier::IntVector<DIM>& ratio_to_level_zero,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry);

   /**
    * Destructor for SkeletonPatchGeometry.
    */
   ~SkeletonPatchGeometry<DIM>();

   /**
    * Print SkeletonPatchGeometry class data.
    */
   void printClassData(ostream& os) const;

private:
   // These are not implemented. 
   SkeletonPatchGeometry(const SkeletonPatchGeometry<DIM>&); 
   void operator=(const SkeletonPatchGeometry<DIM>&);

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SkeletonPatchGeometry.C"
#endif
