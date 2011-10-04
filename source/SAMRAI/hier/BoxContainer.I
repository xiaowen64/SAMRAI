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
#include "SAMRAI/hier/BoxContainerOrderedConstIterator.h"
#include "SAMRAI/hier/BoxContainerOrderedConstReverseIterator.h"
#include "SAMRAI/hier/BoxContainerOrderedIterator.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const tbox::Dimension& dim,
   const int n):
   d_dim(dim),
   d_list(n, Box(dim)),
   d_set_created(false)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const Box& box):
   d_dim(box.getDim()),
   d_list(1, box),
   d_set_created(false)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const BoxContainer& other):
   d_dim(other.d_dim),
   d_list(other.d_list),
   d_set_created(false) 
{
   if (other.d_set_created) {
      makeOrdered();
   }
}

SAMRAI_INLINE_KEYWORD
BoxContainer::~BoxContainer()
{
}

SAMRAI_INLINE_KEYWORD
int BoxContainer::size() const
{
   if (!d_set_created) {
      return static_cast<int>(d_list.size());
   } else {
      return static_cast<int>(d_set.size());
   }
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
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::pushBack(const Box& item)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, item);
   d_list.push_back(item);
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::insertBefore(Iterator iter, const Box& item)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, item);
   d_list.insert(iter.d_list_iter, item);
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::spliceFront(BoxContainer& boxes)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, boxes);
   d_list.splice(begin().d_list_iter, boxes.d_list);
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::spliceBack(BoxContainer& boxes)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, boxes);
   boxes.spliceFront(*this);
   d_list.swap(boxes.d_list);
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popFront()
{
   d_list.pop_front();
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popBack()
{
   d_list.pop_back();
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::erase(Iterator iter)
{
   d_list.erase(iter.d_list_iter);
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::erase(Iterator first, Iterator last)
{
   d_list.erase(first.d_list_iter, last.d_list_iter);
   d_set.clear();
   d_set_created = false;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::clear()
{
   d_list.clear();
   d_set.clear();
   d_set_created = false; 
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::operator == (const BoxContainer& rhs) const
{
   bool isEqual = (d_dim == rhs.d_dim);
   if (isEqual) {
      if (!d_set_created || !rhs.d_set_created) { 
         isEqual = std::equal(d_list.begin(), d_list.end(),
                      rhs.d_list.begin(), Box::box_equality());
      } else {
         isEqual = std::equal(d_set.begin(), d_set.end(),
                      rhs.d_set.begin(), Box::id_equal());
      }
   }
   return isEqual;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::operator != (const BoxContainer& rhs) const
{
   return !(*this == rhs);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::OrderedConstIterator BoxContainer::orderedBegin() const
{
   TBOX_ASSERT(d_set_created || size() == 0);
   return OrderedConstIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::OrderedConstIterator BoxContainer::orderedEnd() const
{
   TBOX_ASSERT(d_set_created || size() == 0);
   return OrderedConstIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::OrderedConstReverseIterator BoxContainer::orderedRBegin() const
{
   TBOX_ASSERT(d_set_created || size() == 0);
   return OrderedConstReverseIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::OrderedConstReverseIterator BoxContainer::orderedREnd() const
{
   TBOX_ASSERT(d_set_created || size() == 0);
   return OrderedConstReverseIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::OrderedIterator BoxContainer::orderedBegin()
{
   if (!d_set_created) {
      makeOrdered();
   }
   return OrderedIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::OrderedIterator BoxContainer::orderedEnd()
{
   if (!d_set_created) {
      makeOrdered();
   }
   return OrderedIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& BoxContainer::getDim() const
{
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::separatePeriodicImages(
   std::vector<Box>& real_mapped_box_vector,
   std::vector<Box>& periodic_image_mapped_box_vector) const
{
   if (!isEmpty()) {

      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      real_mapped_box_vector.reserve(real_mapped_box_vector.size() + size());
      for (OrderedConstIterator ni = orderedBegin(); ni != orderedEnd(); ++ni) {
         const Box& mapped_box = *ni;
         if (mapped_box.getPeriodicId() == zero_shift_number) {
            real_mapped_box_vector.push_back(*ni);
         } else {
            periodic_image_mapped_box_vector.push_back(*ni);
         }
      }
   }
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::isOrdered() const
{
   return d_set_created;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::insert(const Box& box)
{
   if (!d_set_created) {
      makeOrdered();
   }

   const std::list<Box>::iterator& iter = d_list.insert(d_list.end(), box);
   Box* box_ptr(&(*iter));
   if (d_set.insert(box_ptr).second) {
      return true;
   } else {
      d_list.erase(iter);
      return false;
   }
}




}
}
