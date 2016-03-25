//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/cell/PatchCellDataNormOpsComplex.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Norm operations for complex cell-centered patch data.
//

#ifndef included_math_PatchCellDataNormOpsComplex_C
#define included_math_PatchCellDataNormOpsComplex_C

#include "PatchCellDataNormOpsComplex.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchCellDataNormOpsComplex<DIM>::PatchCellDataNormOpsComplex()
{
}

template<int DIM>  PatchCellDataNormOpsComplex<DIM>::~PatchCellDataNormOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchCellDataNormOpsComplex<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box) const
{
   const hier::Box<DIM> ibox = box * data->getGhostBox();
   int retval = ibox.size() * data->getDepth();
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* Norm operations for complex cell-centered data.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::sumControlVolumes(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::CellData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
#endif
   return( d_array_ops.sumControlVolumes(data->getArrayData(),
                                         cvol->getArrayData(),
                                         box) );
}

template<int DIM> void PatchCellDataNormOpsComplex<DIM>::abs(
   tbox::Pointer< pdat::CellData<DIM,double> >& dst,
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   d_array_ops.abs(dst->getArrayData(),
                   src->getArrayData(),
                   box);
}

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::L1Norm(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.L1Norm(data->getArrayData(), box);
   } else {
      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
                                                   cvol->getArrayData(),
                                                   box);
   }
   return( retval );
}

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::L2Norm(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.L2Norm(data->getArrayData(), box);
   } else {
      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
                                                   cvol->getArrayData(),
                                                   box);
   }
   return( retval );
}

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::weightedL2Norm(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
                                          weight->getArrayData(),
                                          box);
   } else {
      retval = d_array_ops.weightedL2NormWithControlVolume(
                           data->getArrayData(),
                           weight->getArrayData(),
                           cvol->getArrayData(),
                           box);
   }
   return( retval );
}

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::RMSNorm(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::weightedRMSNorm(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM> double PatchCellDataNormOpsComplex<DIM>::maxNorm(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval;
   if (cvol.isNull()) {
      retval = d_array_ops.maxNorm(data->getArrayData(), box);
   } else {
      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
                                                    cvol->getArrayData(),
                                                    box);
   }
   return( retval );
}

template<int DIM> dcomplex PatchCellDataNormOpsComplex<DIM>::dot(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data1,
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   dcomplex retval;
   if (cvol.isNull()) {
      retval = d_array_ops.dot(data1->getArrayData(),
                               data2->getArrayData(),
                               box);
   } else {
      retval = d_array_ops.dotWithControlVolume(
                           data1->getArrayData(),
                           data2->getArrayData(),
                           cvol->getArrayData(),
                           box);
   }
   return( retval );
}

template<int DIM> dcomplex PatchCellDataNormOpsComplex<DIM>::integral(
   const tbox::Pointer< pdat::CellData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::CellData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   dcomplex retval = dcomplex(0.0,0.0);

   retval = d_array_ops.integral(data->getArrayData(),
                                 vol->getArrayData(),
                                 box);

   return( retval );
}


}
}
#endif
