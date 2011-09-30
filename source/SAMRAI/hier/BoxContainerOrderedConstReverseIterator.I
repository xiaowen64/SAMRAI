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
BoxContainerOrderedConstReverseIterator::BoxContainerOrderedConstReverseIterator(
   const BoxContainerOrderedConstReverseIterator& other):
   d_box_container(other.d_box_container),
   d_set_iter(other.d_set_iter)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerOrderedConstReverseIterator::~BoxContainerOrderedConstReverseIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerOrderedConstReverseIterator::operator * () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box* BoxContainerOrderedConstReverseIterator::operator -> () const
{
   return *d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerOrderedConstReverseIterator::operator () () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerOrderedConstReverseIterator::operator bool () const
{
   return d_set_iter != d_box_container->setREnd().d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerOrderedConstReverseIterator::operator ! () const
{
   return d_set_iter == d_box_container->setREnd().d_set_iter;
}

SAMRAI_INLINE_KEYWORD
void BoxContainerOrderedConstReverseIterator::operator ++ (
   int)
{
   d_set_iter++;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerOrderedConstReverseIterator& BoxContainerOrderedConstReverseIterator::operator ++ ()
{
   ++d_set_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerOrderedConstReverseIterator::operator == (
   const BoxContainerOrderedConstReverseIterator& other) const
{
   return d_set_iter == other.d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerOrderedConstReverseIterator::operator != (
   const BoxContainerOrderedConstReverseIterator& other) const
{
   return d_set_iter != other.d_set_iter;
}

}
}
