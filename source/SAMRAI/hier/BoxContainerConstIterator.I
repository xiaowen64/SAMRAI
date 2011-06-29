/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A const iterator over containers of boxes.
 *
 ************************************************************************/

#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainerConstIterator& other):
   d_box_container(other.d_box_container),
   d_list_iter(other.d_list_iter)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::~BoxContainerConstIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerConstIterator::operator * () const
{
   return *d_list_iter;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerConstIterator::operator () () const
{
   return *d_list_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::operator bool () const
{
   return d_list_iter != d_box_container->end().d_list_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerConstIterator::operator ! () const
{
   return d_list_iter == d_box_container->end().d_list_iter;
}

SAMRAI_INLINE_KEYWORD
void BoxContainerConstIterator::operator ++ (
   int)
{
   d_list_iter++;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerConstIterator& BoxContainerConstIterator::operator ++ ()
{
   ++d_list_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerConstIterator::operator == (
   const BoxContainerConstIterator& other) const
{
   return d_list_iter == other.d_list_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerConstIterator::operator != (
   const BoxContainerConstIterator& other) const
{
   return d_list_iter != other.d_list_iter;
}

}
}
