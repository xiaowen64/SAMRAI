/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex edge-centered patch data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace math {

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
PatchEdgeDataNormOpsComplex::numberOfEntries(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const int data_depth = data->getDepth();
   for (int d = 0; d < dimVal; d++) {
      retval += (pdat::EdgeGeometry::toEdgeBox(ibox, d).size() * data_depth);
   }
   return retval;
}

SAMRAI_INLINE_KEYWORD
double
PatchEdgeDataNormOpsComplex::sumControlVolumes(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data && cvol);
#endif
   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   for (int d = 0; d < dimVal; d++) {
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
            cvol->getArrayData(d),
            pdat::EdgeGeometry::toEdgeBox(box, d));
   }
   return retval;
}

SAMRAI_INLINE_KEYWORD
void
PatchEdgeDataNormOpsComplex::abs(
   const boost::shared_ptr<pdat::EdgeData<double> >& dst,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; d++) {
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         pdat::EdgeGeometry::toEdgeBox(box, d));
   }
}

SAMRAI_INLINE_KEYWORD
dcomplex
PatchEdgeDataNormOpsComplex::integral(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   int dimVal = box.getDim().getValue();

   dcomplex retval = dcomplex(0.0, 0.0);

   for (int d = 0; d < dimVal; d++) {
      retval += d_array_ops.integral(
            data->getArrayData(d),
            vol->getArrayData(d),
            pdat::EdgeGeometry::toEdgeBox(box, d));
   }

   return retval;
}

}
}
