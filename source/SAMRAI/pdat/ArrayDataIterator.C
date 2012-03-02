/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for array patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataIterator_C
#define included_pdat_ArrayDataIterator_C

#include "SAMRAI/pdat/ArrayDataIterator.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/pdat/ArrayDataIterator.I"
#endif

namespace SAMRAI {
namespace pdat {

ArrayDataIterator::ArrayDataIterator(
   const hier::Box& box):
   d_index(box.lower()),
   d_box(box)
{
}

ArrayDataIterator::ArrayDataIterator(
   const ArrayDataIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

ArrayDataIterator::~ArrayDataIterator()
{
}

ArrayDataIterator::operator bool () const
{
   const tbox::Dimension& dim(d_box.getDim());
   bool retval = true;
   for (int i = 0; i < dim.getValue(); i++) {
      if (d_index(i) > d_box.upper(i)) {
         retval = false;
         break;
      }
   }

   return retval;
}

void
ArrayDataIterator::operator ++ (
   int)
{
   const tbox::Dimension& dim(d_box.getDim());
   d_index(0)++;
   for (int i = 0; i < dim.getValue() - 1; i++) {
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
