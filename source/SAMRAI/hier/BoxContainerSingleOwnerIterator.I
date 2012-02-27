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
BoxContainerSingleOwnerIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_iter != d_mapped_boxes->end() &&
          d_iter->getOwnerRank() == d_owner_rank;
}

SAMRAI_INLINE_KEYWORD
BoxContainerSingleOwnerIterator&
BoxContainerSingleOwnerIterator::operator = (
   const BoxContainerSingleOwnerIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_owner_rank = r.d_owner_rank;
   return *this;
}

SAMRAI_INLINE_KEYWORD
const Box&
BoxContainerSingleOwnerIterator::operator * () const
{
   return *d_iter;
}

SAMRAI_INLINE_KEYWORD
const Box*
BoxContainerSingleOwnerIterator::operator -> () const
{
   return &(*d_iter);
}

SAMRAI_INLINE_KEYWORD
bool
BoxContainerSingleOwnerIterator::operator == (
   const BoxContainerSingleOwnerIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes &&
          d_owner_rank == r.d_owner_rank &&
          d_iter == r.d_iter;
}

SAMRAI_INLINE_KEYWORD
bool
BoxContainerSingleOwnerIterator::operator != (
   const BoxContainerSingleOwnerIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes ||
          d_owner_rank != r.d_owner_rank ||
          d_iter != r.d_iter;
}

}
}
