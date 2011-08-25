/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxSet.
 *
 ************************************************************************/
#ifndef included_hier_BoxSetSingleBlockIterator_C
#define included_hier_BoxSetSingleBlockIterator_C

#include "SAMRAI/hier/BoxSetSingleBlockIterator.h"

namespace SAMRAI {
namespace hier {

BoxSetSingleBlockIterator::BoxSetSingleBlockIterator(
   const BoxSet& mapped_boxes,
   const BlockId& block_id):
   d_mapped_boxes(&mapped_boxes),
   d_block_id(block_id)
{
   d_iter = d_mapped_boxes->begin();
   while (d_iter != d_mapped_boxes->end() && d_iter->getBlockId() != d_block_id) {
      ++d_iter;
   }
}

BoxSetSingleBlockIterator::~BoxSetSingleBlockIterator()
{
   d_mapped_boxes = NULL;
}

bool BoxSetSingleBlockIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_iter != d_mapped_boxes->end() &&
          d_iter->getBlockId() == d_block_id;
}

BoxSetSingleBlockIterator& BoxSetSingleBlockIterator::operator = (
   const BoxSetSingleBlockIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_block_id = r.d_block_id;
   return *this;
}

const Box& BoxSetSingleBlockIterator::operator * () const
{
   return *d_iter;
}

const Box *BoxSetSingleBlockIterator::operator -> () const
{
   return &(*d_iter);
}

bool BoxSetSingleBlockIterator::operator == (
   const BoxSetSingleBlockIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes &&
          d_block_id == r.d_block_id &&
          d_iter == r.d_iter;
}

bool BoxSetSingleBlockIterator::operator != (
   const BoxSetSingleBlockIterator& r) const
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

BoxSetSingleBlockIterator& BoxSetSingleBlockIterator::operator ++ ()
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

BoxSetSingleBlockIterator BoxSetSingleBlockIterator::operator ++ (
   int)
{
   BoxSetSingleBlockIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getBlockId() != d_block_id);
   return saved;
}

int BoxSetSingleBlockIterator::count() const
{
   int ct = 0;
   BoxSetSingleBlockIterator iter(*d_mapped_boxes, d_block_id);
   while (iter.isValid()) {
      ++ct;
      ++iter;
   }
   return ct;
}

}
}
#endif
