/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex side-centered patch data.
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
PatchSideDataNormOpsComplex::numberOfEntries(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const hier::IntVector& directions = data->getDirectionVector();
   const int data_depth = data->getDepth();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         retval += ((pdat::SideGeometry::toSideBox(ibox, d).size()) * data_depth);
      }
   }
   return retval;
}

SAMRAI_INLINE_KEYWORD
double
PatchSideDataNormOpsComplex::sumControlVolumes(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& data,
   const boost::shared_ptr<pdat::SideData<double> >& cvol,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data && cvol);
#endif
   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, cvol->getDirectionVector()));
#endif
   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
               cvol->getArrayData(d),
               pdat::SideGeometry::toSideBox(box, d));
      }
   }
   return retval;
}

SAMRAI_INLINE_KEYWORD
void
PatchSideDataNormOpsComplex::abs(
   const boost::shared_ptr<pdat::SideData<double> >& dst,
   const boost::shared_ptr<pdat::SideData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         d_array_ops.abs(dst->getArrayData(d),
            src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}

SAMRAI_INLINE_KEYWORD
dcomplex
PatchSideDataNormOpsComplex::integral(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   int dimVal = box.getDim().getValue();

   dcomplex retval = dcomplex(0.0, 0.0);
   const hier::IntVector& directions = data->getDirectionVector();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, vol->getDirectionVector()));
#endif
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         retval += d_array_ops.integral(
               data->getArrayData(d),
               vol->getArrayData(d),
               pdat::SideGeometry::toSideBox(box, d));
      }
   }

   return retval;
}

}
}
