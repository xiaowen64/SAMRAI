/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple Skeleton grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_SkeletonPatchGeometry
#define included_geom_SkeletonPatchGeometry

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"

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
 * geometry class is derived from hier::PatchGeometry base class.
 *
 * @see hier::BoundaryBox
 * @see hier::PatchGeometry
 * @see geom::SkeletonGridGeometry
 */

class SkeletonPatchGeometry:
   public hier::PatchGeometry
{
public:
   typedef hier::PatchGeometry::TwoDimBool TwoDimBool;

   /**
    * Constructor for SkeletonPatchGeometry class.  It simply passes
    * patch boundary information and the ratio to the coarsest level to
    * hier::PatchGeometry constructor.
    */
   SkeletonPatchGeometry(
      const hier::IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const TwoDimBool& touches_periodic_bdry);

   /**
    * Destructor for SkeletonPatchGeometry.
    */
   ~SkeletonPatchGeometry();

   /**
    * Print SkeletonPatchGeometry class data.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

private:
   // These are not implemented.
   SkeletonPatchGeometry(
      const SkeletonPatchGeometry&);
   void
   operator = (
      const SkeletonPatchGeometry&);

};

}
}

#endif
