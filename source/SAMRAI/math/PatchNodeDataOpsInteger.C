/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer node-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataOpsInteger_C
#define included_math_PatchNodeDataOpsInteger_C

#include "SAMRAI/math/PatchNodeDataOpsInteger.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include "SAMRAI/tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace math {

PatchNodeDataOpsInteger::PatchNodeDataOpsInteger()
{
}

PatchNodeDataOpsInteger::~PatchNodeDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int PatchNodeDataOpsInteger::numberOfEntries(
   const boost::shared_ptr<pdat::NodeData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box ibox =
      pdat::NodeGeometry::toNodeBox(box * data->getGhostBox());
   int retval = ibox.size() * data->getDepth();
   return retval;
}

/*
 *************************************************************************
 *
 * General operations for integer node-centered patch data.
 *
 *************************************************************************
 */

void PatchNodeDataOpsInteger::swapData(
   boost::shared_ptr<hier::Patch> patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::NodeData<int> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::NodeData<int> > d2(
      patch->getPatchData(data2_id),
      boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));
#endif
   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void PatchNodeDataOpsInteger::printData(
   const boost::shared_ptr<pdat::NodeData<int> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void PatchNodeDataOpsInteger::copyData(
   boost::shared_ptr<pdat::NodeData<int> >& dst,
   const boost::shared_ptr<pdat::NodeData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   (dst->getArrayData()).copy(src->getArrayData(), node_box);
}

void PatchNodeDataOpsInteger::setToScalar(
   boost::shared_ptr<pdat::NodeData<int> >& dst,
   const int& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

void PatchNodeDataOpsInteger::abs(
   boost::shared_ptr<pdat::NodeData<int> >& dst,
   const boost::shared_ptr<pdat::NodeData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      node_box);
}

}
}
#endif
