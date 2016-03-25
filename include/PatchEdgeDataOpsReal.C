//
// File:	PatchEdgeDataOpsReal.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated operations for real edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsReal_C
#define included_math_PatchEdgeDataOpsReal_C

#include "PatchEdgeDataOpsReal.h"
#include "EdgeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
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
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchEdgeDataOpsReal<DIM,TYPE>::operator=(
   const PatchEdgeDataOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
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
   assert(!d1.isNull() && !d2.isNull());
   assert(d1->getDepth() && d2->getDepth());
   assert(d1->getBox() == d2->getBox());
   assert(d1->getGhostBox() == d2->getGhostBox());
#endif
   patch->setPatchData( data1_id, d2 );
   patch->setPatchData( data2_id, d1 );
}

template<int DIM, class TYPE>
void PatchEdgeDataOpsReal<DIM,TYPE>::printData(
   const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   ostream& s) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   s << "Data box = " << box << endl;
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
   assert(!dst.isNull() && !src.isNull());
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
   assert(!dst.isNull());
#endif
   dst->fillAll(alpha, box);
}

}
}
#endif
