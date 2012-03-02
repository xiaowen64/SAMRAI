/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator for side centered patch data types
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

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
const SideIndex&
SideIterator::operator * () const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const SideIndex&
SideIterator::operator () () const
{
   return d_index;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
SAMRAI_INLINE_KEYWORD
SideIterator::operator const void * () const
{
   return SideIterator::operator bool () ? this : NULL;
}
#endif

SAMRAI_INLINE_KEYWORD
bool
SideIterator::operator ! () const
{
   return !SideIterator::operator bool ();
}

SAMRAI_INLINE_KEYWORD
bool
SideIterator::operator == (
   const SideIterator& iter) const
{
   return d_index == iter.d_index;
}

SAMRAI_INLINE_KEYWORD
bool
SideIterator::operator != (
   const SideIterator& iter) const
{
   return d_index != iter.d_index;
}

}
}
