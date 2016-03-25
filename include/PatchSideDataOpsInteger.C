//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/side/PatchSideDataOpsInteger.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Operations for integer side-centered patch data.
//

#ifndef included_math_PatchSideDataOpsInteger_C
#define included_math_PatchSideDataOpsInteger_C

#include "PatchSideDataOpsInteger.h"
#include "SideGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchSideDataOpsInteger<DIM>::PatchSideDataOpsInteger()
{
}

template<int DIM>  PatchSideDataOpsInteger<DIM>::~PatchSideDataOpsInteger()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchSideDataOpsInteger<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::SideData<DIM,int> >& data,
   const hier::Box<DIM>& box) const
{
   int retval = 0;
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   const hier::IntVector<DIM>& directions = data->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> dbox = pdat::SideGeometry<DIM>::toSideBox(ibox, d);
         retval += (dbox.size() * data->getDepth());
      }
   }
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* General operations for integer side-centered patch data.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchSideDataOpsInteger<DIM>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::SideData<DIM,int> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::SideData<DIM,int> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getDirectionVector() == d2->getDirectionVector()); 
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::printData(
   const tbox::Pointer< pdat::SideData<DIM,int> >& data,
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

template<int DIM> void PatchSideDataOpsInteger<DIM>::copyData(
   tbox::Pointer< pdat::SideData<DIM,int> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         (dst->getArrayData(d)).copy(src->getArrayData(d), side_box);
      }
   }
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::setToScalar(
   tbox::Pointer< pdat::SideData<DIM,int> >& dst,
   const int& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

template<int DIM> void PatchSideDataOpsInteger<DIM>::abs(
   tbox::Pointer< pdat::SideData<DIM,int> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,int> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
#endif
   const hier::IntVector<DIM>& directions = dst->getDirectionVector();
   for (int d = 0; d < DIM; d++) {
      if (directions(d)) {
         const hier::Box<DIM> side_box = pdat::SideGeometry<DIM>::toSideBox(box, d);
         d_array_ops.abs(dst->getArrayData(d),
                         src->getArrayData(d),
                         side_box);
      }
   }
}


}
}
#endif
