/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outerside data objects
 *
 ************************************************************************/

#ifndef included_pdat_OutersideDataFactory_C
#define included_pdat_OutersideDataFactory_C

#include "SAMRAI/pdat/OutersideDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/SideDataFactory.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/OutersideDataFactory.I"
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
OutersideDataFactory<TYPE>::OutersideDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim))
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OutersideDataFactory<TYPE>::~OutersideDataFactory()
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
boost::shared_ptr<hier::PatchDataFactory>
OutersideDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return boost::shared_ptr<hier::PatchDataFactory>(
      new OutersideDataFactory<TYPE>(ghosts.getDim(), d_depth));
}

/*
 *************************************************************************
 *
 * Allocate the concrete outerside data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
OutersideDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* patchdata =
      new OutersideData<TYPE>(patch.getBox(), d_depth);
   return boost::shared_ptr<hier::PatchData>(patchdata);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outerside data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
OutersideDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   hier::BoxGeometry* boxgeometry = new OutersideGeometry(box, zero_vector);
   return boost::shared_ptr<hier::BoxGeometry>(boxgeometry);
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory needed to allocate the data object.
 *
 *************************************************************************
 */

template<class TYPE>
size_t OutersideDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OutersideData<TYPE>));
   const size_t data = OutersideData<TYPE>::getSizeOfData(box, d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from NodeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool OutersideDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are SideData and OutersideData.
    */
   if (!valid_copy) {
      boost::shared_ptr<SideDataFactory<TYPE> > sdf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (sdf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      boost::shared_ptr<OutersideDataFactory<TYPE> > osdf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (osdf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
