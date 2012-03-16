/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for cell centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_CellIterator_C
#define included_pdat_CellIterator_C

#include "SAMRAI/pdat/CellIterator.h"

namespace SAMRAI {
namespace pdat {

CellIterator::CellIterator(
   const hier::Box& box):
   d_index(box.lower()),
   d_box(box)
{
}

CellIterator::CellIterator(
   const CellIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

CellIterator::~CellIterator()
{
}

CellIterator::operator bool () const
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
CellIterator::operator ++ (
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
