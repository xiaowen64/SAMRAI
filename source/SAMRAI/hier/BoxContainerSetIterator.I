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
#include "SAMRAI/hier/BoxContainerSetConstIterator.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainerSetIterator::BoxContainerSetIterator(
   const BoxContainerSetIterator& other):
   d_box_container(other.d_box_container),
   d_set_iter(other.d_set_iter)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerSetIterator::~BoxContainerSetIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerSetIterator::operator * () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box* BoxContainerSetIterator::operator -> () const
{
   return *d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerSetIterator::operator () () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerSetIterator::operator bool () const
{
   return d_set_iter != d_box_container->setEnd().d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerSetIterator::operator ! () const
{
   return d_set_iter == d_box_container->setEnd().d_set_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerSetIterator BoxContainerSetIterator::operator ++ (
   int)
{
   BoxContainerSetIterator return_iter(*this);
   ++d_set_iter;
   return return_iter;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerSetIterator& BoxContainerSetIterator::operator ++ ()
{
   ++d_set_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxContainerSetIterator BoxContainerSetIterator::operator -- (
   int)
{
   BoxContainerSetIterator return_iter(*this);
   --d_set_iter;
   return return_iter;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerSetIterator& BoxContainerSetIterator::operator -- ()
{
   --d_set_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerSetIterator::operator == (
   const BoxContainerSetIterator& other) const
{
   return d_set_iter == other.d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerSetIterator::operator != (
   const BoxContainerSetIterator& other) const
{
   return d_set_iter != other.d_set_iter;
}

}
}
