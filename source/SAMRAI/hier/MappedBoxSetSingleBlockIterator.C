/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for MappedBoxSet.
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxSetSingleBlockIterator_C
#define included_hier_MappedBoxSetSingleBlockIterator_C

#include "SAMRAI/hier/MappedBoxSetSingleBlockIterator.h"

namespace SAMRAI {
namespace hier {

MappedBoxSetSingleBlockIterator::MappedBoxSetSingleBlockIterator(
   const MappedBoxSet& mapped_boxes,
   const BlockId &block_id):
   d_mapped_boxes(&mapped_boxes),
   d_block_id(block_id)
{
   d_iter = d_mapped_boxes->begin();
   while (d_iter != d_mapped_boxes->end() && d_iter->getBlockId() != d_block_id) {
      ++d_iter;
   }
}

MappedBoxSetSingleBlockIterator::~MappedBoxSetSingleBlockIterator()
{
   d_mapped_boxes = NULL;
}

bool MappedBoxSetSingleBlockIterator::isValid() const
{
   return
      d_mapped_boxes != NULL &&
      d_iter != d_mapped_boxes->end() &&
      d_iter->getBlockId() == d_block_id;
}

MappedBoxSetSingleBlockIterator& MappedBoxSetSingleBlockIterator::operator = (
   const MappedBoxSetSingleBlockIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_iter = r.d_iter;
   d_block_id = r.d_block_id;
   return *this;
}

const MappedBox& MappedBoxSetSingleBlockIterator::operator * () const
{
   return *d_iter;
}

const MappedBox *MappedBoxSetSingleBlockIterator::operator -> () const
{
   return &(*d_iter);
}

bool MappedBoxSetSingleBlockIterator::operator == (
   const MappedBoxSetSingleBlockIterator& r) const
{
   return
      d_mapped_boxes == r.d_mapped_boxes &&
      d_block_id == r.d_block_id &&
      d_iter == r.d_iter;
}

bool MappedBoxSetSingleBlockIterator::operator != (
   const MappedBoxSetSingleBlockIterator& r) const
{
   return
      d_mapped_boxes != r.d_mapped_boxes ||
      d_block_id != r.d_block_id ||
      d_iter != r.d_iter;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

MappedBoxSetSingleBlockIterator &MappedBoxSetSingleBlockIterator::operator ++ ()
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

MappedBoxSetSingleBlockIterator MappedBoxSetSingleBlockIterator::operator ++ (
   int)
{
   MappedBoxSetSingleBlockIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getBlockId() != d_block_id);
   return saved;
}

int MappedBoxSetSingleBlockIterator::count() const
{
   int ct = 0;
   MappedBoxSetSingleBlockIterator iter(*d_mapped_boxes, d_block_id);
   while (iter.isValid()) {
      ++ct;
      ++iter;
   }
   return ct;
}

}
}
#endif
