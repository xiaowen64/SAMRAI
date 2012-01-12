/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating cell data objects
 *
 ************************************************************************/

#ifndef included_pdat_CellDataFactory_C
#define included_pdat_CellDataFactory_C

#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/hier/Patch.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/CellDataFactory.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the default ghost cell width and depth.
 *
 *************************************************************************
 */

template<class TYPE>
CellDataFactory<TYPE>::CellDataFactory(
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_mb_trans(NULL)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   d_mb_trans = NULL;
}

template<class TYPE>
CellDataFactory<TYPE>::~CellDataFactory()
{
   if (d_mb_trans) {
      delete d_mb_trans;
   }
}

/*
 *************************************************************************
 *
 * Clone the factory and copy the default parameters to the new factory.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchDataFactory>
CellDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return tbox::Pointer<hier::PatchDataFactory>(new CellDataFactory<TYPE>(
                                                   d_depth, ghosts));
}

/*
 *************************************************************************
 *
 * Allocate the concrete cell data classes.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchData>
CellDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* patchdata =
      new CellData<TYPE>(patch.getBox(), this->d_depth, this->d_ghosts);
   return tbox::Pointer<hier::PatchData>(patchdata);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for cell data objects.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::BoxGeometry>
CellDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   hier::BoxGeometry* boxgeometry = new CellGeometry(box, this->d_ghosts);
   return tbox::Pointer<hier::BoxGeometry>(boxgeometry);
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory needed to allocate the data object.
 *
 *************************************************************************
 */

template<class TYPE>
size_t CellDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(CellData<TYPE>));
   const size_t data =
      CellData<TYPE>::getSizeOfData(box, d_depth, this->d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from CellData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool CellDataFactory<TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Only valid option is CellData.
    */
   tbox::Pointer<CellDataFactory<TYPE> > cdf(
      dst_pdf,
      tbox::__dynamic_cast_tag());
   if (cdf) {
      valid_copy = true;
   }
   return valid_copy;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
