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
BoxContainerIterator::BoxContainerIterator(
   const BoxContainerIterator& other):
   d_box_container(other.d_box_container),
   d_list_iter(other.d_list_iter)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerIterator::~BoxContainerIterator()
{
}

SAMRAI_INLINE_KEYWORD
Box& BoxContainerIterator::operator * () const
{
   return *d_list_iter;
}

SAMRAI_INLINE_KEYWORD
Box& BoxContainerIterator::operator () () const
{
   return *d_list_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerIterator::operator bool () const
{
   return d_list_iter != d_box_container->end().d_list_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerIterator::operator ! () const
{
   return d_list_iter == d_box_container->end().d_list_iter;
}

SAMRAI_INLINE_KEYWORD
void BoxContainerIterator::operator ++ (
   int)
{
   d_list_iter++;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerIterator& BoxContainerIterator::operator ++ ()
{
   ++d_list_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerIterator::operator == (
   const BoxContainerIterator& other) const
{
   return d_list_iter == other.d_list_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerIterator::operator != (
   const BoxContainerIterator& other) const
{
   return d_list_iter != other.d_list_iter;
}

}
}
