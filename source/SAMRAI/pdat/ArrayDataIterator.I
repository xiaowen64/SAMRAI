/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for array patch data types
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
ArrayDataIterator::ArrayDataIterator(
   const hier::Box& box):
   d_index(box.lower()),
   d_box(box)
{
}

SAMRAI_INLINE_KEYWORD
ArrayDataIterator::ArrayDataIterator(
   const ArrayDataIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SAMRAI_INLINE_KEYWORD
ArrayDataIterator&
ArrayDataIterator::operator = (
   const ArrayDataIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
ArrayDataIterator::~ArrayDataIterator()
{
}

SAMRAI_INLINE_KEYWORD
const hier::Index& ArrayDataIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const hier::Index& ArrayDataIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
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

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
ArrayDataIterator::operator const void * () const
{
   return ArrayDataIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool ArrayDataIterator::operator ! () const
{
   return !ArrayDataIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void ArrayDataIterator::operator ++ (
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

SAMRAI_INLINE_KEYWORD
bool ArrayDataIterator::operator == (
   const ArrayDataIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool ArrayDataIterator::operator != (
   const ArrayDataIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
