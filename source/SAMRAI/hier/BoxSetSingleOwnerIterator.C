/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxSet.
 *
 ************************************************************************/
#ifndef included_hier_BoxSetSingleOwnerIterator_C
#define included_hier_BoxSetSingleOwnerIterator_C

#include "SAMRAI/hier/BoxSetSingleOwnerIterator.h"

namespace SAMRAI {
namespace hier {

BoxSetSingleOwnerIterator::BoxSetSingleOwnerIterator(
   const BoxSet& mapped_boxes,
   const int& owner_rank):
   d_mapped_boxes(&mapped_boxes),
   d_owner_rank(owner_rank),
   d_iter(d_mapped_boxes->orderedBegin())
{
   while (d_iter != d_mapped_boxes->orderedEnd() && d_iter->getOwnerRank() != d_owner_rank) {
      ++d_iter;
   }
}

BoxSetSingleOwnerIterator::~BoxSetSingleOwnerIterator()
{
   d_mapped_boxes = NULL;
}

bool BoxSetSingleOwnerIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_iter != d_mapped_boxes->orderedEnd() &&
          d_iter->getOwnerRank() == d_owner_rank;
}

BoxSetSingleOwnerIterator& BoxSetSingleOwnerIterator::operator = (
   const BoxSetSingleOwnerIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_owner_rank = r.d_owner_rank;
   return *this;
}

const Box& BoxSetSingleOwnerIterator::operator * () const
{
   return *d_iter;
}

const Box *BoxSetSingleOwnerIterator::operator -> () const
{
   return &(*d_iter);
}

bool BoxSetSingleOwnerIterator::operator == (
   const BoxSetSingleOwnerIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes &&
          d_owner_rank == r.d_owner_rank &&
          d_iter == r.d_iter;
}

bool BoxSetSingleOwnerIterator::operator != (
   const BoxSetSingleOwnerIterator& r) const
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

BoxSetSingleOwnerIterator& BoxSetSingleOwnerIterator::operator ++ ()
{
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->orderedEnd() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

BoxSetSingleOwnerIterator BoxSetSingleOwnerIterator::operator ++ (
   int)
{
   BoxSetSingleOwnerIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->orderedEnd() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return saved;
}

}
}
#endif
