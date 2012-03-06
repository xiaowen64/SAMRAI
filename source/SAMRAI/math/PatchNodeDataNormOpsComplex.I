/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex node-centered patch data.
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
PatchNodeDataNormOpsComplex::numberOfEntries(
   const boost::shared_ptr<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box ibox =
      pdat::NodeGeometry::toNodeBox(box * data->getGhostBox());
   return ibox.size() * data->getDepth();
}

/*
 *************************************************************************
 *
 * Norm operations for complex node-centered data.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
double
PatchNodeDataNormOpsComplex::sumControlVolumes(
   const boost::shared_ptr<pdat::NodeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::NodeData<double> >& cvol,
   const hier::Box& box) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data && cvol);
#endif
   return d_array_ops.sumControlVolumes(data->getArrayData(),
      cvol->getArrayData(),
      pdat::NodeGeometry::toNodeBox(box));
}

SAMRAI_INLINE_KEYWORD
void
PatchNodeDataNormOpsComplex::abs(
   const boost::shared_ptr<pdat::NodeData<double> >& dst,
   const boost::shared_ptr<pdat::NodeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      pdat::NodeGeometry::toNodeBox(box));
}

SAMRAI_INLINE_KEYWORD
dcomplex
PatchNodeDataNormOpsComplex::integral(
   const boost::shared_ptr<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::NodeData<double> >& vol) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   return d_array_ops.integral(
      data->getArrayData(),
      vol->getArrayData(),
      pdat::NodeGeometry::toNodeBox(box));
}

}
}
