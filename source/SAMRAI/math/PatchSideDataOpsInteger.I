/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer side-centered patch data.
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
PatchSideDataOpsInteger::numberOfEntries(
   const boost::shared_ptr<pdat::SideData<int> >& data,
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
void PatchSideDataOpsInteger::printData(
   const boost::shared_ptr<pdat::SideData<int> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

SAMRAI_INLINE_KEYWORD
void PatchSideDataOpsInteger::copyData(
   const boost::shared_ptr<pdat::SideData<int> >& dst,
   const boost::shared_ptr<pdat::SideData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         dst->getArrayData(d).copy(src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}

SAMRAI_INLINE_KEYWORD
void PatchSideDataOpsInteger::setToScalar(
   const boost::shared_ptr<pdat::SideData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

SAMRAI_INLINE_KEYWORD
void PatchSideDataOpsInteger::abs(
   const boost::shared_ptr<pdat::SideData<int> >& dst,
   const boost::shared_ptr<pdat::SideData<int> >& src,
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

}
}
