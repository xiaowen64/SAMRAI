//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/cell/PatchCellDataMiscellaneousOpsReal.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated miscellaneous operations for real cell-centered data.
//

#ifndef included_math_PatchCellDataMiscellaneousOpsReal_C
#define included_math_PatchCellDataMiscellaneousOpsReal_C

#include "PatchCellDataMiscellaneousOpsReal.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::PatchCellDataMiscellaneousOpsReal()
{
}

template<int DIM, class TYPE>
PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::~PatchCellDataMiscellaneousOpsReal()
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
PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::PatchCellDataMiscellaneousOpsReal(
   const PatchCellDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::operator=(
   const PatchCellDataMiscellaneousOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* Templated miscellaneous operations for real cell-centered data.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::computeConstrProdPos(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   int retval;
   if (cvol.isNull()) {
      retval = d_array_ops.computeConstrProdPos(data1->getArrayData(),
                                                data2->getArrayData(),
                                                box);
   } else {
      retval = d_array_ops.computeConstrProdPosWithControlVolume(
                           data1->getArrayData(),
                           data2->getArrayData(),
                           cvol->getArrayData(),
                           box);
   }
   return( retval );
}

template<int DIM, class TYPE>
void PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::compareToScalar(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const TYPE& alpha,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   if (cvol.isNull()) {
      d_array_ops.compareToScalar(dst->getArrayData(),
                                  src->getArrayData(),
                                  alpha,
                                  box);
   } else {
      d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(),
                                                   src->getArrayData(),
                                                   alpha,
                                                   cvol->getArrayData(),
                                                   box);
   }
}

template<int DIM, class TYPE>
int PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::testReciprocal(
   tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   int retval;
   if (cvol.isNull()) {
      retval = d_array_ops.testReciprocal(dst->getArrayData(),
                                          src->getArrayData(),
                                          box);
   } else {
      retval = d_array_ops.testReciprocalWithControlVolume(
                           dst->getArrayData(),
                           src->getArrayData(),
                           cvol->getArrayData(),
                           box);
   }
   return( retval );
}

template<int DIM, class TYPE>
TYPE PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::maxPointwiseDivide(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   retval = d_array_ops.maxPointwiseDivide(numer->getArrayData(),
					   denom->getArrayData(),
					   box);
   return( retval );
}


template <int DIM, class TYPE>
TYPE PatchCellDataMiscellaneousOpsReal<DIM,TYPE>::minPointwiseDivide(
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& numer,
   const tbox::Pointer< pdat::CellData<DIM,TYPE> >& denom,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!numer.isNull() && !denom.isNull());
#endif
   TYPE retval;
   retval = d_array_ops.minPointwiseDivide(numer->getArrayData(),
					   denom->getArrayData(),
					   box);
   return( retval );
}

}
}
#endif
