/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for node centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_NodeIterator_C
#define included_pdat_NodeIterator_C

#include "SAMRAI/pdat/NodeIterator.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/NodeIterator.I"
#endif

namespace SAMRAI {
namespace pdat {

NodeIterator::NodeIterator(
   const hier::Box& box):
   d_index(box.lower(), hier::IntVector::getZero(box.getDim())),
   d_box(NodeGeometry::toNodeBox(box))
{
}

NodeIterator::NodeIterator(
   const NodeIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

NodeIterator::~NodeIterator()
{
}

NodeIterator::operator bool () const
{
   bool retval = true;
   for (int i = 0; i < d_box.getDim().getValue(); i++) {
      if (d_index(i) > d_box.upper(i)) {
         retval = false;
         break;
      }
   }

   return retval;
}

void
NodeIterator::operator ++ (
   int)
{
   d_index(0)++;
   for (int i = 0; i < d_box.getDim().getValue() - 1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i + 1)++;
      } else {
         break;
      }
   }
}

}
}
#endif
