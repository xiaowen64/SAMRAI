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
#include "SAMRAI/hier/BoxContainerIterator.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainer& container,
   bool from_start):
   d_list_iter(from_start ? container.d_list.begin() :
               container.d_list.end()),
   d_set_iter(from_start ? container.d_set.begin() :
              container.d_set.end()),
   d_ordered(container.d_ordered)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainerConstIterator& other):
   d_list_iter(other.d_list_iter),
   d_set_iter(other.d_set_iter),
   d_ordered(other.d_ordered)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainerIterator& other):
   d_list_iter(other.d_list_iter),
   d_set_iter(other.d_set_iter),
   d_ordered(other.d_ordered)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator::~BoxContainerConstIterator()
{
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator&
BoxContainerConstIterator::operator = (
   const BoxContainerConstIterator& rhs)
{
   if (this != &rhs) {
      d_list_iter = rhs.d_list_iter;
      d_set_iter = rhs.d_set_iter;
      d_ordered = rhs.d_ordered;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerConstIterator::operator * () const
{
   return d_ordered ? **d_set_iter : *d_list_iter;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainerConstIterator::operator () () const
{
   return d_ordered ? **d_set_iter : *d_list_iter;
}

SAMRAI_INLINE_KEYWORD
const Box* BoxContainerConstIterator::operator -> () const
{
   return d_ordered ? *d_set_iter : &(*d_list_iter);
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator
BoxContainerConstIterator::operator ++ ( int )
{
   BoxContainerConstIterator return_iter(*this);
   if (d_ordered) {
      ++d_set_iter;
   } else {
      ++d_list_iter;
   }
   return return_iter;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerConstIterator& BoxContainerConstIterator::operator ++ ()
{
   if (d_ordered) {
      ++d_set_iter;
   } else {
      ++d_list_iter;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxContainerConstIterator
BoxContainerConstIterator::operator -- ( int )
{
   BoxContainerConstIterator return_iter(*this);
   if (d_ordered) {
      --d_set_iter;
   } else {
      --d_list_iter;
   }
   return return_iter;
}

SAMRAI_INLINE_KEYWORD
const BoxContainerConstIterator& BoxContainerConstIterator::operator -- ()
{
   if (d_ordered) {
      --d_set_iter;
   } else {
      --d_list_iter;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerConstIterator::operator == (
   const BoxContainerConstIterator& other) const
{
   return d_ordered ? d_set_iter == other.d_set_iter :
                      d_list_iter == other.d_list_iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainerConstIterator::operator != (
   const BoxContainerConstIterator& other) const
{
   return d_ordered ? d_set_iter != other.d_set_iter :
                      d_list_iter != other.d_list_iter;
}

}
}
