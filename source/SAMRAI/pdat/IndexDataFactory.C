/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   IndexDataFactory implementation
 *
 ************************************************************************/

#ifndef included_pdat_IndexDataFactory_C
#define included_pdat_IndexDataFactory_C

#include "SAMRAI/pdat/IndexDataFactory.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MemoryUtilities.h"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the default ghost cell width.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
IndexDataFactory<TYPE, BOX_GEOMETRY>::IndexDataFactory(
   const hier::IntVector& ghosts):
   hier::PatchDataFactory(ghosts)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataFactory<TYPE, BOX_GEOMETRY>::~IndexDataFactory()
{
}

/*
 *************************************************************************
 *
 * Clone the factory and copy the default parameters to the new factory.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
tbox::Pointer<hier::PatchDataFactory>
IndexDataFactory<TYPE, BOX_GEOMETRY>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return tbox::Pointer<hier::PatchDataFactory>(new IndexDataFactory<TYPE,
                                                                     BOX_GEOMETRY>(
                                                   ghosts));
}

/*
 *************************************************************************
 *
 * Allocate the concrete irregular data class.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
tbox::Pointer<hier::PatchData>
IndexDataFactory<TYPE, BOX_GEOMETRY>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* pd = new IndexData<TYPE, BOX_GEOMETRY>(patch.getBox(), this->d_ghosts);
   return tbox::Pointer<hier::PatchData>(pd);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for index data objects.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
tbox::Pointer<hier::BoxGeometry>
IndexDataFactory<TYPE, BOX_GEOMETRY>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   hier::BoxGeometry* boxgeometry = new BOX_GEOMETRY(box, this->d_ghosts);
   return tbox::Pointer<hier::BoxGeometry>(boxgeometry);
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory needed to allocate the object.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
size_t IndexDataFactory<TYPE, BOX_GEOMETRY>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   NULL_USE(box);
   return tbox::MemoryUtilities::align(sizeof(IndexData<TYPE, BOX_GEOMETRY>));
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from IndexData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
bool IndexDataFactory<TYPE, BOX_GEOMETRY>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid option is another IndexData object of the same dimension
    * and type.
    */
   if (!valid_copy) {
      tbox::Pointer<IndexDataFactory<TYPE, BOX_GEOMETRY> > idf(
         dst_pdf, tbox::__dynamic_cast_tag());
      if (idf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
