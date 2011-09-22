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

#include "SAMRAI/hier/BoxContainerSetIterator.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainerSetConstIterator::BoxContainerSetConstIterator(
   const BoxContainerSetConstIterator& other):
   d_box_container(other.d_box_container),
   d_set_iter(other.d_set_iter)
{
}


SAMRAI_INLINE_KEYWORD
BoxContainerSetConstIterator::BoxContainerSetConstIterator(
   const BoxContainerSetIterator& other):
   d_box_container(other.d_box_container),
   d_set_iter(other.d_set_iter)
{
}


SAMRAI_INLINE_KEYWORD
BoxContainerSetConstIterator::~BoxContainerSetConstIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerSetConstIterator::operator * () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box* BoxContainerSetConstIterator::operator -> () const
{
   return *d_set_iter;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerSetConstIterator::operator () () const
{
   return **d_set_iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainerSetConstIterator::operator bool () const
{
   return d_set_iter != d_box_container->setEnd().d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerSetConstIterator::operator ! () const
{
   return d_set_iter == d_box_container->setEnd().d_set_iter;
}

SAMRAI_INLINE_KEYWORD
void BoxContainerSetConstIterator::operator ++ (
   int)
{
   d_set_iter++;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerSetConstIterator& BoxContainerSetConstIterator::operator ++ ()
{
   ++d_set_iter;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerSetConstIterator::operator == (
   const BoxContainerSetConstIterator& other) const
{
   return d_set_iter == other.d_set_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerSetConstIterator::operator != (
   const BoxContainerSetConstIterator& other) const
{
   return d_set_iter != other.d_set_iter;
}

}
}
