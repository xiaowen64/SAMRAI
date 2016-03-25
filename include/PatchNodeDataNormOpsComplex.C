//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/node/PatchNodeDataNormOpsComplex.C $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Norm operations for complex node-centered patch data.
//

#ifndef included_math_PatchNodeDataNormOpsComplex_C
#define included_math_PatchNodeDataNormOpsComplex_C

#include "PatchNodeDataNormOpsComplex.h"
#include "NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace math {

template<int DIM>  PatchNodeDataNormOpsComplex<DIM>::PatchNodeDataNormOpsComplex()
{
}

template<int DIM>  PatchNodeDataNormOpsComplex<DIM>::~PatchNodeDataNormOpsComplex()
{
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchNodeDataNormOpsComplex<DIM>::numberOfEntries(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box) const
{
   const hier::Box<DIM> ibox = 
                   pdat::NodeGeometry<DIM>::toNodeBox(box * data->getGhostBox());
   int retval = ibox.size() * data->getDepth();
   return( retval );
}

/*
*************************************************************************
*                                                                       *
* Norm operations for complex node-centered data.                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::sumControlVolumes(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::NodeData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !cvol.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   return( d_array_ops.sumControlVolumes(data->getArrayData(),
                                         cvol->getArrayData(),
                                         node_box) );
}

template<int DIM> void PatchNodeDataNormOpsComplex<DIM>::abs(
   tbox::Pointer< pdat::NodeData<DIM,double> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.abs(dst->getArrayData(),
                   src->getArrayData(),
                   node_box);
}

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::L1Norm(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.L1Norm(data->getArrayData(), node_box);
   } else {
      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
                                                   cvol->getArrayData(),
                                                   node_box);
   }
   return( retval );
}

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::L2Norm(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.L2Norm(data->getArrayData(), node_box);
   } else {
      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
                                                   cvol->getArrayData(),
                                                   node_box);
   }
   return( retval );
}

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::weightedL2Norm(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull() && !weight.isNull());
#endif
   double retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
                                          weight->getArrayData(),
                                          node_box);
   } else {
      retval = d_array_ops.weightedL2NormWithControlVolume(
                           data->getArrayData(),
                           weight->getArrayData(),
                           cvol->getArrayData(),
                           node_box);
   }
   return( retval );
}

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::RMSNorm(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
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

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::weightedRMSNorm(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
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

template<int DIM> double PatchNodeDataNormOpsComplex<DIM>::maxNorm(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   double retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.maxNorm(data->getArrayData(), node_box);
   } else {
      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
                                                    cvol->getArrayData(),
                                                    node_box);
   }
   return( retval );
}

template<int DIM> dcomplex PatchNodeDataNormOpsComplex<DIM>::dot(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data1,
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data1.isNull() && !data2.isNull());
#endif
   dcomplex retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   if (cvol.isNull()) {
      retval = d_array_ops.dot(data1->getArrayData(),
                               data2->getArrayData(),
                               node_box);
   } else {
      retval = d_array_ops.dotWithControlVolume(
                           data1->getArrayData(),
                           data2->getArrayData(),
                           cvol->getArrayData(),
                           node_box);
   }
   return( retval );
}

template<int DIM> dcomplex PatchNodeDataNormOpsComplex<DIM>::integral(
   const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   dcomplex retval;
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);

   retval = d_array_ops.integral(
                        data->getArrayData(),
                        vol->getArrayData(),
                        node_box);

   return( retval );
}

}
}
#endif
