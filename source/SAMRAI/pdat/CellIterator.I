/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for cell centered patch data types 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
CellIterator::CellIterator(
   const hier::Box& box):
   d_index(box.lower()),
   d_box(box)
{
}

SAMRAI_INLINE_KEYWORD
CellIterator::CellIterator(
   const CellIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SAMRAI_INLINE_KEYWORD
CellIterator&
CellIterator::operator = (
   const CellIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
CellIterator::~CellIterator()
{
}

SAMRAI_INLINE_KEYWORD
const CellIndex& CellIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const CellIndex& CellIterator::operator () () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
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

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
CellIterator::operator const void* () const
{
   return CellIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool CellIterator::operator ! () const
{
   return !CellIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
void CellIterator::operator ++ (
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
bool CellIterator::operator == (
   const CellIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool CellIterator::operator != (
   const CellIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
