/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
const hier::BoxList&
CellOverlap::getDestinationBoxList() const
{
   return d_dst_boxes;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& CellOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

SAMRAI_INLINE_KEYWORD
const hier::Transformation& CellOverlap::getTransformation() const
{
   return d_transformation;
}


}
}
