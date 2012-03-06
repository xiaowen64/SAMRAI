/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer edge-centered patch data.
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
PatchEdgeDataOpsInteger::numberOfEntries(
   const boost::shared_ptr<pdat::EdgeData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const int data_depth = data->getDepth();
   for (int d = 0; d < dimVal; d++) {
      retval += ((pdat::EdgeGeometry::toEdgeBox(ibox, d).size()) * data_depth);
   }
   return retval;
}

SAMRAI_INLINE_KEYWORD
void
PatchEdgeDataOpsInteger::printData(
   const boost::shared_ptr<pdat::EdgeData<int> >& data,
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
void
PatchEdgeDataOpsInteger::copyData(
   const boost::shared_ptr<pdat::EdgeData<int> >& dst,
   const boost::shared_ptr<pdat::EdgeData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; d++) {
      dst->getArrayData(d).copy(src->getArrayData(d),
         pdat::EdgeGeometry::toEdgeBox(box, d));
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchEdgeDataOpsInteger::setToScalar(
   const boost::shared_ptr<pdat::EdgeData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

SAMRAI_INLINE_KEYWORD
void
PatchEdgeDataOpsInteger::abs(
   const boost::shared_ptr<pdat::EdgeData<int> >& dst,
   const boost::shared_ptr<pdat::EdgeData<int> >& src,
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

}
}
