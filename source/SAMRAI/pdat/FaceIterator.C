/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for face centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_FaceIterator_C
#define included_pdat_FaceIterator_C

#include "SAMRAI/pdat/FaceIterator.h"

namespace SAMRAI {
namespace pdat {

FaceIterator::FaceIterator(
   const hier::Box& box,
   const int axis):
   d_index(box.lower(), axis, FaceIndex::Lower),
   d_box(FaceGeometry::toFaceBox(box, axis))
{
}

FaceIterator::FaceIterator(
   const FaceIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

FaceIterator::~FaceIterator()
{
}

FaceIterator::operator bool () const
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
FaceIterator::operator ++ (
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
