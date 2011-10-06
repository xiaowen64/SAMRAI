/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#include "SAMRAI/hier/BoxContainer.h"


namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainerOrderedConstIterator::BoxContainerOrderedConstIterator(
   const BoxContainerOrderedConstIterator& other):
   d_box_container(other.d_box_container),
   d_set_iter(other.d_set_iter)
{
}


SAMRAI_INLINE_KEYWORD
BoxContainerOrderedConstIterator::~BoxContainerOrderedConstIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerOrderedConstIterator::operator * () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box* BoxContainerOrderedConstIterator::operator -> () const
{
   return *d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerOrderedConstIterator::operator () () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerOrderedConstIterator BoxContainerOrderedConstIterator::operator ++ (
   int)
{
   BoxContainerOrderedConstIterator return_iter(*this);
   ++d_set_iter;
   return return_iter;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerOrderedConstIterator& BoxContainerOrderedConstIterator::operator ++ ()
{
   ++d_set_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxContainerOrderedConstIterator BoxContainerOrderedConstIterator::operator -- (
   int)
{
   BoxContainerOrderedConstIterator return_iter(*this);
   --d_set_iter;
   return return_iter;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerOrderedConstIterator& BoxContainerOrderedConstIterator::operator -- ()
{
   --d_set_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerOrderedConstIterator::operator == (
   const BoxContainerOrderedConstIterator& other) const
{
   return d_set_iter == other.d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerOrderedConstIterator::operator != (
   const BoxContainerOrderedConstIterator& other) const
{
   return d_set_iter != other.d_set_iter;
}

}
}
