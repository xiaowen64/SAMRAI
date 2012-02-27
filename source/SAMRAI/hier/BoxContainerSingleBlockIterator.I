/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
bool
BoxContainerSingleBlockIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_iter != d_mapped_boxes->end() &&
          d_iter->getBlockId() == d_block_id;
}

SAMRAI_INLINE_KEYWORD
BoxContainerSingleBlockIterator&
BoxContainerSingleBlockIterator::operator = (
   const BoxContainerSingleBlockIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_block_id = r.d_block_id;
   return *this;
}

SAMRAI_INLINE_KEYWORD
const Box&
BoxContainerSingleBlockIterator::operator * () const
{
   return *d_iter;
}

SAMRAI_INLINE_KEYWORD
const Box*
BoxContainerSingleBlockIterator::operator -> () const
{
   return &(*d_iter);
}

SAMRAI_INLINE_KEYWORD
bool
BoxContainerSingleBlockIterator::operator == (
   const BoxContainerSingleBlockIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes &&
          d_block_id == r.d_block_id &&
          d_iter == r.d_iter;
}

SAMRAI_INLINE_KEYWORD
bool
BoxContainerSingleBlockIterator::operator != (
   const BoxContainerSingleBlockIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes ||
          d_block_id != r.d_block_id ||
          d_iter != r.d_iter;
}

}
}
