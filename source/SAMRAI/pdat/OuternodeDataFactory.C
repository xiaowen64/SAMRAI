/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outernode data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeDataFactory_C
#define included_pdat_OuternodeDataFactory_C

#include "SAMRAI/pdat/OuternodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/pdat/OuternodeData.h"
#include "SAMRAI/pdat/OuternodeGeometry.h"
#include "SAMRAI/hier/Patch.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/OuternodeDataFactory.I"
#endif

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the depth of the patch data.
 *
 *************************************************************************
 */

template<class TYPE>
OuternodeDataFactory<TYPE>::OuternodeDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim))
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OuternodeDataFactory<TYPE>::~OuternodeDataFactory()
{
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
OuternodeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return tbox::Pointer<hier::PatchDataFactory>(new OuternodeDataFactory<TYPE>(
                                                   ghosts.getDim(), d_depth));
}

/*
 *************************************************************************
 *
 * Allocate the concrete outernode data classes.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchData>
OuternodeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* patchdata =
      new OuternodeData<TYPE>(patch.getBox(), d_depth);
   return tbox::Pointer<hier::PatchData>(patchdata);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outernode data objects.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::BoxGeometry>
OuternodeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   hier::BoxGeometry* boxgeometry = new OuternodeGeometry(box, zero_vector);
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
size_t OuternodeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OuternodeData<TYPE>));
   const size_t data = OuternodeData<TYPE>::getSizeOfData(box,
         d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from OuternodeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool OuternodeDataFactory<TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are NodeData and OuternodeData.
    */
   if (!valid_copy) {
      tbox::Pointer<NodeDataFactory<TYPE> > ndf = dst_pdf;
      if (!ndf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer<OuternodeDataFactory<TYPE> > ondf = dst_pdf;
      if (!ondf.isNull()) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
