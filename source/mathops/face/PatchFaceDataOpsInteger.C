//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/face/PatchFaceDataOpsInteger.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Operations for integer face-centered patch data.
//

#ifndef included_math_PatchFaceDataOpsInteger_C
#define included_math_PatchFaceDataOpsInteger_C

#include "PatchFaceDataOpsInteger.h"
#include "FaceGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchFaceDataOpsInteger<DIM>::PatchFaceDataOpsInteger()
{
}

template<int DIM>  PatchFaceDataOpsInteger<DIM>::~PatchFaceDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchFaceDataOpsInteger<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::FaceData<DIM,int> >& data,
   const hier::Box<DIM>& box) const
{
   int retval = 0;
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> dbox = pdat::FaceGeometry<DIM>::toFaceBox(ibox, d);
      retval += (dbox.size() * data->getDepth());
   }
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* General operations for integer face-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchFaceDataOpsInteger<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::FaceData<DIM,int> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::FaceData<DIM,int> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchFaceDataOpsInteger<DIM>::printData(
   const tbox::Pointer< pdat::FaceData<DIM,int> >& data,
   const hier::Box<DIM>& box,
   std::ostream& s) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

template<int DIM> void PatchFaceDataOpsInteger<DIM>::copyData(
   tbox::Pointer< pdat::FaceData<DIM,int> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), face_box);
   }
}

template<int DIM> void PatchFaceDataOpsInteger<DIM>::setToScalar(
   tbox::Pointer< pdat::FaceData<DIM,int> >& dst,
   const int& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

template<int DIM> void PatchFaceDataOpsInteger<DIM>::abs(
   tbox::Pointer< pdat::FaceData<DIM,int> >& dst,
   const tbox::Pointer< pdat::FaceData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> face_box = pdat::FaceGeometry<DIM>::toFaceBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
                      src->getArrayData(d),
                      face_box);
   }
}


}
}
#endif
