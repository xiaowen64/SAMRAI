//
// File:	PatchNodeDataNormOpsReal.C
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated norm operations for real node-centered patch data.
//

#ifndef included_math_PatchNodeDataNormOpsReal_C
#define included_math_PatchNodeDataNormOpsReal_C

#include "PatchNodeDataNormOpsReal.h"
#include "NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace math {

template<int DIM, class TYPE>
PatchNodeDataNormOpsReal<DIM,TYPE>::PatchNodeDataNormOpsReal()
{
}

template<int DIM, class TYPE>
PatchNodeDataNormOpsReal<DIM,TYPE>::~PatchNodeDataNormOpsReal()
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
PatchNodeDataNormOpsReal<DIM,TYPE>::PatchNodeDataNormOpsReal(
   const PatchNodeDataNormOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

template<int DIM, class TYPE>
void PatchNodeDataNormOpsReal<DIM,TYPE>::operator=(
   const PatchNodeDataNormOpsReal<DIM,TYPE>& foo)
{
   (void) foo;  // not implemented (but needed by some compilers)
}

/*
*************************************************************************
*                                                                       *
* Compute the number of data entries on a patch in the given box.       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
int PatchNodeDataNormOpsReal<DIM,TYPE>::numberOfEntries(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
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
* Templated norm operations for real node-centered data.                *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::sumControlVolumes(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::NodeData<DIM,double> >& cvol,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !cvol.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   return( d_array_ops.sumControlVolumes(data->getArrayData(),
                                         cvol->getArrayData(),
                                         node_box) );
}

template<int DIM, class TYPE>
void PatchNodeDataNormOpsReal<DIM,TYPE>::abs(
   tbox::Pointer< pdat::NodeData<DIM,TYPE> >& dst,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& src,
   const hier::Box<DIM>& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst.isNull() && !src.isNull());
#endif
   const hier::Box<DIM> node_box = pdat::NodeGeometry<DIM>::toNodeBox(box);
   d_array_ops.abs(dst->getArrayData(),
                   src->getArrayData(),
                   node_box);
}

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::L1Norm(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
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

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::L2Norm(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
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

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::weightedL2Norm(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !weight.isNull());
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

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::RMSNorm(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   double retval = L2Norm(data, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::weightedRMSNorm(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& weight,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull() && !weight.isNull());
#endif
   double retval = weightedL2Norm(data, weight, box, cvol);
   if (cvol.isNull()) {
      retval /= sqrt( (double)numberOfEntries(data, box) );
   } else {
      retval /= sqrt( sumControlVolumes(data, cvol, box) );
   }
   return( retval );
}

template<int DIM, class TYPE>
double PatchNodeDataNormOpsReal<DIM,TYPE>::maxNorm(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
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

template<int DIM, class TYPE>
TYPE PatchNodeDataNormOpsReal<DIM,TYPE>::dot(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data1,
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data2,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > cvol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data1.isNull() && !data2.isNull());
#endif
   TYPE retval;
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

template<int DIM, class TYPE>
TYPE PatchNodeDataNormOpsReal<DIM,TYPE>::integral(
   const tbox::Pointer< pdat::NodeData<DIM,TYPE> >& data,
   const hier::Box<DIM>& box,
   const tbox::Pointer< pdat::NodeData<DIM,double> > vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!data.isNull());
#endif
   TYPE retval;
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
