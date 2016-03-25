//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/mathops/edge/PatchEdgeDataOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Templated operations for real edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsReal_C
#define included_math_PatchEdgeDataOpsReal_C

#include "PatchEdgeDataOpsReal.h"
#include "EdgeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchEdgeDataOpsReal<DIM,TYPE>::PatchEdgeDataOpsReal()
{
}

template<int DIM, class TYPE>
PatchEdgeDataOpsReal<DIM,TYPE>::~PatchEdgeDataOpsReal()
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
PatchEdgeDataOpsReal<DIM,TYPE>::PatchEdgeDataOpsReal(
   const PatchEdgeDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

template<int DIM, class TYPE>
void PatchEdgeDataOpsReal<DIM,TYPE>::operator=(
   const PatchEdgeDataOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo); 
}

/*
*************************************************************************
*                                                                       *
* General templated operations for real edge-centered patch data.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void PatchEdgeDataOpsReal<DIM,TYPE>::swapData(
   tbox::Pointer< hier::Patch<DIM> > patch,
   const int data1_id,
   const int data2_id) const
{
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d1 = patch->getPatchData(data1_id);
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> > d2 = patch->getPatchData(data2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d1.isNull() && !d2.isNull());
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox() == d2->getBox());
   TBOX_ASSERT(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM, class TYPE>
void PatchEdgeDataOpsReal<DIM,TYPE>::printData(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
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
void PatchEdgeDataOpsReal<DIM,TYPE>::copyData(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   for (int d = 0; d < DIM; d++) {
      const hier::Box<DIM> edge_box = pdat::EdgeGeometry<DIM>::toEdgeBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), edge_box);
   }
}

template<int DIM, class TYPE>
void PatchEdgeDataOpsReal<DIM,TYPE>::setToScalar(
   tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
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
