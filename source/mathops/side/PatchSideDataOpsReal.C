//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/side/PatchSideDataOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated operations for real side-centered patch data.
//

#ifndef included_math_PatchSideDataOpsReal_C
#define included_math_PatchSideDataOpsReal_C

#include "PatchSideDataOpsReal.h"
#include "SideGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchSideDataOpsReal<DIM,TYPE>::PatchSideDataOpsReal()
{
}

template<int DIM, class TYPE>
PatchSideDataOpsReal<DIM,TYPE>::~PatchSideDataOpsReal()
{
}

/*
*************************************************************************
*                                                                       *
* The const constructor and assignment operator are not actually used   *
* but are defined here for compilers that require an implementation for *
* every declaration.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
PatchSideDataOpsReal<DIM,TYPE>::PatchSideDataOpsReal(
   const PatchSideDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void PatchSideDataOpsReal<DIM,TYPE>::operator=(
   const PatchSideDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* General templated operations for real side-centered patch data.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchSideDataOpsReal<DIM,TYPE>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::SideData<DIM,TYPE> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::SideData<DIM,TYPE> > d2 = patch->getPatchData(data2_id);
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

template<int DIM, class TYPE>
void PatchSideDataOpsReal<DIM,TYPE>::printData(
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& data,
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

template<int DIM, class TYPE>
void PatchSideDataOpsReal<DIM,TYPE>::copyData(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::SideData<DIM,TYPE> >& src,
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

template<int DIM, class TYPE>
void PatchSideDataOpsReal<DIM,TYPE>::setToScalar(
   tbox::Pointer< pdat::SideData<DIM,TYPE> >& dst,
   const TYPE& alpha,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

}
}
#endif
