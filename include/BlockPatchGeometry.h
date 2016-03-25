//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/geometry/skeleton/patch_geom/BlockPatchGeometry.h $
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description: Patch geometry for multiblock.
//

#ifndef included_geom_BlockPatchGeometry
#define included_geom_BlockPatchGeometry

#include "SAMRAI_config.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"

namespace SAMRAI {
    namespace geom {

/**
 * Class BlockPatchGeometry implements geometry management 
 * for a single patch in an AMR hierarchy with no information about
 * the physical characteristics of the problem.  This is intended for
 * use in an application that will manage the physical geometry in
 * user-defined code.
 * 
 * The grid data is set by BlockGridGeometry class.  This patch
 * geometry class is derived from hier::PatchGeometry<DIM> base class.
 *
 * @see hier::BoundaryBox
 * @see hier::PatchGeometry
 * @see geom::BlockGridGeometry
 */

template<int DIM> class BlockPatchGeometry 
: public hier::PatchGeometry<DIM>
{
public:

   /**
    * Constructor for BlockPatchGeometry class.  It simply passes 
    * patch boundary information and the ratio to the coarsest level to
    * hier::PatchGeometry constructor.
    */
   BlockPatchGeometry<DIM>(
      const hier::IntVector<DIM>& ratio_to_level_zero,
      const int block_number,
      const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
      const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry);

   /**
    * Destructor for BlockPatchGeometry.
    */
   ~BlockPatchGeometry<DIM>();

   /**
    * Get the block number for the block on which the patch lies. 
    */
   int getBlockNumber() const
   {
      return (d_block_number);
   }

   /**
    * Print BlockPatchGeometry class data.
    */
   virtual void printClassData(std::ostream& os) const;

private:
   // These are not implemented. 
   BlockPatchGeometry(const BlockPatchGeometry<DIM>&); 
   void operator=(const BlockPatchGeometry<DIM>&);

   int d_block_number;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BlockPatchGeometry.C"
#endif
