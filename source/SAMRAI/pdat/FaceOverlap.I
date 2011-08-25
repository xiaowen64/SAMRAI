/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
const hier::BoxList&
FaceOverlap::getDestinationBoxList(
   const int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < d_dst_boxes.getSize()));

   return d_dst_boxes[axis];
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& FaceOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

SAMRAI_INLINE_KEYWORD
const hier::Transformation& FaceOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
