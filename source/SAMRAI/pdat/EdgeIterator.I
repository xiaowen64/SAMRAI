/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for edge centered patch data types
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

SAMRAI_INLINE_KEYWORD
EdgeIterator&
EdgeIterator::operator = (
   const EdgeIterator& iter)
{
   d_index = iter.d_index;
   d_box = iter.d_box;
   return *this;
}

SAMRAI_INLINE_KEYWORD
const EdgeIndex&
EdgeIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const EdgeIndex&
EdgeIterator::operator () () const
{
   return d_index;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
EdgeIterator::operator const void * () const
{
   return EdgeIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool
EdgeIterator::operator ! () const
{
   return !EdgeIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
bool
EdgeIterator::operator == (
   const EdgeIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool
EdgeIterator::operator != (
   const EdgeIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
