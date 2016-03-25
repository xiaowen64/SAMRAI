//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/node/PatchNodeDataMiscellaneousOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated miscellaneous operations for real node-centered data.
//

#ifndef included_math_PatchNodeDataMiscellaneousOpsReal_C
#define included_math_PatchNodeDataMiscellaneousOpsReal_C

#include "PatchNodeDataMiscellaneousOpsReal.h"
#include "NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::PatchNodeDataMiscellaneousOpsReal()
{
}

template<int DIM, class TYPE>
PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::~PatchNodeDataMiscellaneousOpsReal()
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
PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::PatchNodeDataMiscellaneousOpsReal(
   const PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

template<int DIM, class TYPE>
void PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::operator=(
   const PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   NULL_USE(foo);
}

/*
*************************************************************************
*                                                                       *
* Templated miscellaneous opertions for real node-centered data.        * 
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPos(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   int retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.computeConstrProdPos(data1->getArrayData(),
                                                data2->getArrayData(),
                                                node_box);
   } else {
      retval = d_array_ops.computeConstrProdPosWithControlVolume(
                           data1->getArrayData(),
                           data2->getArrayData(),
                           cvol->getArrayData(),
                           node_box);
   }
   return( retval );
}

template<int DIM, class TYPE>
void PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalar(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      d_array_ops.compareToScalar(dst->getArrayData(),
                                  src->getArrayData(),
                                  alpha,
                                  node_box);
   } else {
      d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(),
                                                   src->getArrayData(),
                                                   alpha,
                                                   cvol->getArrayData(),
                                                   node_box);
   }
}

template<int DIM, class TYPE>
int PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocal(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   int retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.testReciprocal(dst->getArrayData(),
                                          src->getArrayData(),
                                          node_box);
   } else {
      retval = d_array_ops.testReciprocalWithControlVolume(
                           dst->getArrayData(),
                           src->getArrayData(),
                           cvol->getArrayData(),
                           node_box);
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   retval = d_array_ops.maxPointwiseDivide(numer->getArrayData(),
					   denom->getArrayData(),
					   node_box);
   return( retval );
}

template <int DIM, class TYPE>
TYPE PatchNodeDataMiscellaneousOpsReal<DIM,TYPE>::minPointwiseDivide(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   retval = d_array_ops.minPointwiseDivide(numer->getArrayData(),
					   denom->getArrayData(),
					   node_box);
   return( retval );
}


}
}
#endif
