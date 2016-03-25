/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerSingleOwnerIterator_C
#define included_hier_BoxContainerSingleOwnerIterator_C

#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"

namespace SAMRAI {
namespace hier {

BoxContainerSingleOwnerIterator::BoxContainerSingleOwnerIterator(
   const BoxContainer& mapped_boxes,
   const int& owner_rank):
   d_mapped_boxes(&mapped_boxes),
   d_owner_rank(owner_rank),
   d_iter(d_mapped_boxes->begin())
{
   while (d_iter != d_mapped_boxes->end() && d_iter->getOwnerRank() != d_owner_rank) {
      ++d_iter;
   }
}

BoxContainerSingleOwnerIterator::~BoxContainerSingleOwnerIterator()
{
   d_mapped_boxes = NULL;
}

bool BoxContainerSingleOwnerIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_iter != d_mapped_boxes->end() &&
          d_iter->getOwnerRank() == d_owner_rank;
}

BoxContainerSingleOwnerIterator& BoxContainerSingleOwnerIterator::operator = (
   const BoxContainerSingleOwnerIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_owner_rank = r.d_owner_rank;
   return *this;
}

const Box& BoxContainerSingleOwnerIterator::operator * () const
{
   return *d_iter;
}

const Box *BoxContainerSingleOwnerIterator::operator -> () const
{
   return &(*d_iter);
}

bool BoxContainerSingleOwnerIterator::operator == (
   const BoxContainerSingleOwnerIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes &&
          d_owner_rank == r.d_owner_rank &&
          d_iter == r.d_iter;
}

bool BoxContainerSingleOwnerIterator::operator != (
   const BoxContainerSingleOwnerIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes ||
          d_owner_rank != r.d_owner_rank ||
          d_iter != r.d_iter;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

BoxContainerSingleOwnerIterator& BoxContainerSingleOwnerIterator::operator ++ ()
{
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

BoxContainerSingleOwnerIterator BoxContainerSingleOwnerIterator::operator ++ (
   int)
{
   BoxContainerSingleOwnerIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return saved;
}

}
}
#endif
