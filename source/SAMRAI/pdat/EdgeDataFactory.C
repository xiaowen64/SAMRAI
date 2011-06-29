/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating edge data objects 
 *
 ************************************************************************/

#ifndef included_pdat_EdgeDataFactory_C
#define included_pdat_EdgeDataFactory_C

#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/MemoryUtilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/EdgeDataFactory.I"
#endif

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *									*
 * The constructor simply caches the default ghost cell width and depth.	*
 *									*
 *************************************************************************
 */

template<class TYPE>
EdgeDataFactory<TYPE>::EdgeDataFactory(
   int depth,
   const hier::IntVector& ghosts,
   bool fine_boundary_represents_var):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_mb_trans(NULL)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
EdgeDataFactory<TYPE>::~EdgeDataFactory()
{
   if (d_mb_trans) {
      delete d_mb_trans;
   }
}

/*
 *************************************************************************
 *									*
 * Clone the factory and copy the default parameters to the new factory.	*
 *									*
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchDataFactory>
EdgeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return tbox::Pointer<hier::PatchDataFactory>(new EdgeDataFactory<TYPE>(
                                                   d_depth,
                                                   ghosts,
                                                   d_fine_boundary_represents_var));
}

/*
 *************************************************************************
 *
 * Allocate the concrete edge data classes.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchData>
EdgeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* patchdata =
      new EdgeData<TYPE>(patch.getBox(), d_depth, this->d_ghosts);
   return tbox::Pointer<hier::PatchData>(patchdata);
}

/*
 *************************************************************************
 *									*
 * Return the box geometry type for edge data objects.			*
 *									*
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::BoxGeometry>
EdgeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   hier::BoxGeometry* boxgeometry = new EdgeGeometry(box, this->d_ghosts);
   return tbox::Pointer<hier::BoxGeometry>(boxgeometry);
}

/*
 *************************************************************************
 *									*
 * Calculate the amount of memory needed to allocate the data object.	*
 *									*
 *************************************************************************
 */

template<class TYPE>
size_t EdgeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(EdgeData<TYPE>));
   const size_t data =
      EdgeData<TYPE>::getSizeOfData(box, d_depth, this->d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *									*
 * Determine whether this is a valid copy operation to/from EdgeData     *
 * between the supplied datatype.                                        *
 *									*
 *************************************************************************
 */

template<class TYPE>
bool EdgeDataFactory<TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are EdgeData and OuteredgeData.
    */
   if (!valid_copy) {
      tbox::Pointer<EdgeDataFactory<TYPE> > edf = dst_pdf;
      if (!edf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer<OuteredgeDataFactory<TYPE> > oedf = dst_pdf;
      if (!oedf.isNull()) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
