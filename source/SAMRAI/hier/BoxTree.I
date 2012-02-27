/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Binary tree of Boxes for overlap searches.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
void
BoxTree::clear()
{
   d_bounding_box.setEmpty();
   d_left_child.reset();
   d_right_child.reset();
   d_boxes.clear();
   d_center_child.reset();
}

SAMRAI_INLINE_KEYWORD
bool
BoxTree::isInitialized() const
{
   return !d_bounding_box.empty();
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
BoxTree::getDim() const
{
   return d_dim;
}

}
}
