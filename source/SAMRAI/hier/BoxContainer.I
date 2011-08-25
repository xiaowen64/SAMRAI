/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const tbox::Dimension& dim,
   const int n):
   d_dim(dim),
   d_list(n, Box(dim))
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const Box& box):
   d_dim(box.getDim()),
   d_list(1, box)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const BoxContainer& other):
   d_dim(other.d_dim),
   d_list(other.d_list)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::~BoxContainer()
{
}

SAMRAI_INLINE_KEYWORD
int BoxContainer::size() const
{
   return static_cast<int>(d_list.size());
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::isEmpty() const
{
   return size() == 0;
}

SAMRAI_INLINE_KEYWORD
BoxContainer::ConstIterator BoxContainer::begin() const
{
   return BoxContainerConstIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::ConstIterator BoxContainer::end() const
{
   return BoxContainerConstIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::begin()
{
   return BoxContainerIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::end()
{
   return BoxContainerIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
Box& BoxContainer::front()
{
   return d_list.front();
}

SAMRAI_INLINE_KEYWORD
Box& BoxContainer::back()
{
   return d_list.back();
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainer::front() const
{
   return d_list.front();
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainer::back() const
{
   return d_list.back();
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::pushFront(const Box& item)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, item);
   d_list.push_front(item);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::pushBack(const Box& item)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, item);
   d_list.push_back(item);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::insertBefore(Iterator iter, const Box& item)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, item);
   d_list.insert(iter.d_list_iter, item);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::spliceFront(BoxContainer& boxes)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, boxes);
   d_list.splice(begin().d_list_iter, boxes.d_list);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::spliceBack(BoxContainer& boxes)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, boxes);
   boxes.spliceFront(*this);
   d_list.swap(boxes.d_list);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popFront()
{
   d_list.pop_front();
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popBack()
{
   d_list.pop_back();
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::erase(Iterator iter)
{
   d_list.erase(iter.d_list_iter);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::erase(Iterator first, Iterator last)
{
   d_list.erase(first.d_list_iter, last.d_list_iter);
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::clear()
{
   d_list.clear();
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::operator == (const BoxContainer& rhs) const
{
   bool isEqual = (d_dim == rhs.d_dim);
   if (isEqual) {
      isEqual = std::equal(d_list.begin(), d_list.end(),
            rhs.d_list.begin(), Box::box_equality());
   }
   return isEqual;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::operator != (const BoxContainer& rhs) const
{
   return !(*this == rhs);
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& BoxContainer::getDim() const
{
   return d_dim;
}

}
}
