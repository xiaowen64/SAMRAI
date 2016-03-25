/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating face data objects
 *
 ************************************************************************/

#ifndef included_pdat_FaceDataFactory_C
#define included_pdat_FaceDataFactory_C

#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/OuterfaceDataFactory.h"
#include "SAMRAI/hier/Patch.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/FaceDataFactory.I"
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
FaceDataFactory<TYPE>::FaceDataFactory(
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
FaceDataFactory<TYPE>::~FaceDataFactory()
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
FaceDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return tbox::Pointer<hier::PatchDataFactory>(new FaceDataFactory<TYPE>(
                                                   d_depth,
                                                   ghosts,
                                                   d_fine_boundary_represents_var));
}

/*
 *************************************************************************
 *
 * Allocate the concrete face data classes.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::PatchData>
FaceDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   hier::PatchData* patchdata =
      new FaceData<TYPE>(patch.getBox(), d_depth, this->d_ghosts);
   return tbox::Pointer<hier::PatchData>(patchdata);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for face data objects.
 *
 *************************************************************************
 */

template<class TYPE>
tbox::Pointer<hier::BoxGeometry>
FaceDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   hier::BoxGeometry* boxgeometry = new FaceGeometry(box, this->d_ghosts);
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
size_t FaceDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(FaceData<TYPE>));
   const size_t data =
      FaceData<TYPE>::getSizeOfData(box, d_depth, this->d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from FaceData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool FaceDataFactory<TYPE>::validCopyTo(
   const tbox::Pointer<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are FaceData and OuterfaceData.
    */
   if (!valid_copy) {
      tbox::Pointer<FaceDataFactory<TYPE> > fdf = dst_pdf;
      if (!fdf.isNull()) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      tbox::Pointer<OuterfaceDataFactory<TYPE> > ofdf = dst_pdf;
      if (!ofdf.isNull()) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
