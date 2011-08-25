/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for face centered patch data types
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
FaceIterator::FaceIterator(
   const hier::Box& box,
   const int axis):
   d_index(box.lower(), axis, FaceIndex::Lower),
   d_box(FaceGeometry::toFaceBox(box, axis))
{
}

SAMRAI_INLINE_KEYWORD
FaceIterator::FaceIterator(
   const FaceIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SAMRAI_INLINE_KEYWORD
FaceIterator&
FaceIterator::operator = (
   const FaceIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
FaceIterator::~FaceIterator()
{
}

SAMRAI_INLINE_KEYWORD
const FaceIndex& FaceIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const FaceIndex& FaceIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
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

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
FaceIterator::operator const void * () const
{
   return FaceIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool FaceIterator::operator ! () const
{
   return !FaceIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void FaceIterator::operator ++ (
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

SAMRAI_INLINE_KEYWORD
bool FaceIterator::operator == (
   const FaceIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool FaceIterator::operator != (
   const FaceIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
