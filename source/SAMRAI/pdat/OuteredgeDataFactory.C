/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outeredge data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeDataFactory_C
#define included_pdat_OuteredgeDataFactory_C

#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
#include "SAMRAI/pdat/OuteredgeGeometry.h"
#include "SAMRAI/hier/Patch.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/OuteredgeDataFactory.I"
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
OuteredgeDataFactory<TYPE>::OuteredgeDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim))
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OuteredgeDataFactory<TYPE>::~OuteredgeDataFactory()
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
OuteredgeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return tbox::Pointer<hier::PatchDataFactory>(new OuteredgeDataFactory<TYPE>(
                                                   ghosts.getDim(), d_depth));
}

/*
 *************************************************************************
 *
 * Allocate the concrete outeredge data classes.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchData>
OuteredgeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* patchdata =
      new OuteredgeData<TYPE>(patch.getBox(), d_depth);
   return tbox::Pointer<hier::PatchData>(patchdata);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outeredge data objects.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::BoxGeometry>
OuteredgeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   hier::BoxGeometry* boxgeometry = new OuteredgeGeometry(box, zero_vector);
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
size_t OuteredgeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OuteredgeData<TYPE>));
   const size_t data = OuteredgeData<TYPE>::getSizeOfData(box,
         d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from EdgeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool OuteredgeDataFactory<TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are EdgeData and OuteredgeData.
    */
   if (!valid_copy) {
      tbox::Pointer<EdgeDataFactory<TYPE> > edf(
         dst_pdf,
         tbox::__dynamic_cast_tag());
      if (edf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer<OuteredgeDataFactory<TYPE> > oedf(
         dst_pdf,
         tbox::__dynamic_cast_tag());
      if (oedf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
