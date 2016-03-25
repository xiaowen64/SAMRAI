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
const hier::BoxContainer&
EdgeOverlap::getDestinationBoxContainer(
   const int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < d_dst_boxes.getSize()));

   return d_dst_boxes[axis];
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& EdgeOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

SAMRAI_INLINE_KEYWORD
const hier::Transformation& EdgeOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
