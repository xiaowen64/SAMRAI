/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerSingleBlockIterator_C
#define included_hier_BoxContainerSingleBlockIterator_C

#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"

namespace SAMRAI {
namespace hier {

BoxContainerSingleBlockIterator::BoxContainerSingleBlockIterator(
   const BoxContainer& mapped_boxes,
   const BlockId& block_id):
   d_mapped_boxes(&mapped_boxes),
   d_block_id(block_id),
   d_iter(d_mapped_boxes->begin())
{
   while (d_iter != d_mapped_boxes->end() && d_iter->getBlockId() != d_block_id) {
      ++d_iter;
   }
}

BoxContainerSingleBlockIterator::~BoxContainerSingleBlockIterator()
{
   d_mapped_boxes = NULL;
}

bool BoxContainerSingleBlockIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_iter != d_mapped_boxes->end() &&
          d_iter->getBlockId() == d_block_id;
}

BoxContainerSingleBlockIterator& BoxContainerSingleBlockIterator::operator = (
   const BoxContainerSingleBlockIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_block_id = r.d_block_id;
   return *this;
}

const Box& BoxContainerSingleBlockIterator::operator * () const
{
   return *d_iter;
}

const Box *BoxContainerSingleBlockIterator::operator -> () const
{
   return &(*d_iter);
}

bool BoxContainerSingleBlockIterator::operator == (
   const BoxContainerSingleBlockIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes &&
          d_block_id == r.d_block_id &&
          d_iter == r.d_iter;
}

bool BoxContainerSingleBlockIterator::operator != (
   const BoxContainerSingleBlockIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes ||
          d_block_id != r.d_block_id ||
          d_iter != r.d_iter;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

BoxContainerSingleBlockIterator& BoxContainerSingleBlockIterator::operator ++ ()
{
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getBlockId() != d_block_id);
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

BoxContainerSingleBlockIterator BoxContainerSingleBlockIterator::operator ++ (
   int)
{
   BoxContainerSingleBlockIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getBlockId() != d_block_id);
   return saved;
}

int BoxContainerSingleBlockIterator::count() const
{
   int ct = 0;
   BoxContainerSingleBlockIterator iter(*d_mapped_boxes, d_block_id);
   while (iter.isValid()) {
      ++ct;
      ++iter;
   }
   return ct;
}

}
}
#endif
