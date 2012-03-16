/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for edge centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_EdgeIterator_C
#define included_pdat_EdgeIterator_C

#include "SAMRAI/pdat/EdgeIterator.h"

namespace SAMRAI {
namespace pdat {

EdgeIterator::EdgeIterator(
   const hier::Box& box,
   const int axis):
   d_index(box.lower(), axis, 0),
   d_box(EdgeGeometry::toEdgeBox(box, axis))
{
}

EdgeIterator::EdgeIterator(
   const EdgeIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

EdgeIterator::~EdgeIterator()
{
}

EdgeIterator::operator bool () const
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
EdgeIterator::operator ++ (
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
