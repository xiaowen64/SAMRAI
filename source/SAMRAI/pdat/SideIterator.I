/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Iterator for side centered patch data types 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
SideIterator::SideIterator(
   const hier::Box& box,
   const int axis):
   d_index(box.lower(), axis, SideIndex::Lower),
   d_box(SideGeometry::toSideBox(box, axis))
{
}

SAMRAI_INLINE_KEYWORD
SideIterator::SideIterator(
   const SideIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SAMRAI_INLINE_KEYWORD
SideIterator&
SideIterator::operator = (
   const SideIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
SideIterator::~SideIterator()
{
}

SAMRAI_INLINE_KEYWORD
const SideIndex& SideIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const SideIndex& SideIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
SideIterator::operator bool () const
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
SideIterator::operator const void* () const
{
   return SideIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool SideIterator::operator ! () const
{
   return !SideIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void SideIterator::operator ++ (
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
bool SideIterator::operator == (
   const SideIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool SideIterator::operator != (
   const SideIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
