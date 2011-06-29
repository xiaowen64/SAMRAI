/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for MappedBoxSet.
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxSetSingleOwnerIterator_C
#define included_hier_MappedBoxSetSingleOwnerIterator_C

#include "SAMRAI/hier/MappedBoxSetSingleOwnerIterator.h"

namespace SAMRAI {
namespace hier {

MappedBoxSetSingleOwnerIterator::MappedBoxSetSingleOwnerIterator(
   const MappedBoxSet& mapped_boxes,
   const int &owner_rank):
   d_mapped_boxes(&mapped_boxes),
   d_owner_rank(owner_rank)
{
   d_iter = d_mapped_boxes->begin();
   while (d_iter != d_mapped_boxes->end() && d_iter->getOwnerRank() != d_owner_rank) {
      ++d_iter;
   }
}

MappedBoxSetSingleOwnerIterator::~MappedBoxSetSingleOwnerIterator()
{
   d_mapped_boxes = NULL;
}

bool MappedBoxSetSingleOwnerIterator::isValid() const
{
   return
      d_mapped_boxes != NULL &&
      d_iter != d_mapped_boxes->end() &&
      d_iter->getOwnerRank() == d_owner_rank;
}

MappedBoxSetSingleOwnerIterator& MappedBoxSetSingleOwnerIterator::operator = (
   const MappedBoxSetSingleOwnerIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_owner_rank = r.d_owner_rank;
   return *this;
}

const MappedBox& MappedBoxSetSingleOwnerIterator::operator * () const
{
   return *d_iter;
}

const MappedBox *MappedBoxSetSingleOwnerIterator::operator -> () const
{
   return &(*d_iter);
}

bool MappedBoxSetSingleOwnerIterator::operator == (
   const MappedBoxSetSingleOwnerIterator& r) const
{
   return
      d_mapped_boxes == r.d_mapped_boxes &&
      d_owner_rank == r.d_owner_rank &&
      d_iter == r.d_iter;
}

bool MappedBoxSetSingleOwnerIterator::operator != (
   const MappedBoxSetSingleOwnerIterator& r) const
{
   return
      d_mapped_boxes != r.d_mapped_boxes ||
      d_owner_rank != r.d_owner_rank ||
      d_iter != r.d_iter;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

MappedBoxSetSingleOwnerIterator &MappedBoxSetSingleOwnerIterator::operator ++ ()
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

MappedBoxSetSingleOwnerIterator MappedBoxSetSingleOwnerIterator::operator ++ (
   int)
{
   MappedBoxSetSingleOwnerIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return saved;
}

}
}
#endif
