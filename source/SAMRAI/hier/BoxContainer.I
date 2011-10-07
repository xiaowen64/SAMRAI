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
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer() :
   d_set_created(false)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(bool ordered) :
   d_set_created(ordered)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const Box& box,
   bool ordered):
   d_list(1, box),
   d_set_created(ordered)
{
   if (d_set_created) {
      makeOrdered();
   }
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const BoxContainer& other):
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
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), item);
   }
#endif
   d_list.push_front(item);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::pushBack(const Box& item)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), item);
   }
#endif
   d_list.push_back(item);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::insertBefore(Iterator iter, const Box& item)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), item);
   }
#endif
   d_list.insert(iter.d_list_iter, item);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::spliceFront(BoxContainer& boxes)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0 && boxes.size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), boxes.front());
   }
#endif
   d_list.splice(begin().d_list_iter, boxes.d_list);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::spliceBack(BoxContainer& boxes)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0 && boxes.size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), boxes.front());
   }
#endif
   boxes.spliceFront(*this);
   d_list.swap(boxes.d_list);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popFront()
{
   d_list.pop_front();
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popBack()
{
   d_list.pop_back();
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::erase(Iterator iter)
{
   d_list.erase(iter.d_list_iter);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::erase(Iterator first, Iterator last)
{
   d_list.erase(first.d_list_iter, last.d_list_iter);
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }
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
   bool isEqual;
   if (!d_set_created || !rhs.d_set_created) { 
      isEqual = std::equal(d_list.begin(), d_list.end(),
                   rhs.d_list.begin(), Box::box_equality());
   } else {
      isEqual = std::equal(d_set.begin(), d_set.end(),
                   rhs.d_set.begin(), Box::id_equal());
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
void BoxContainer::separatePeriodicImages(
   std::vector<Box>& real_mapped_box_vector,
   std::vector<Box>& periodic_image_mapped_box_vector) const
{
   TBOX_ASSERT(d_set_created || size() == 0);

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
   return (d_set_created && size() != 0);
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::isUnordered() const
{
   return (!d_set_created && size() != 0);
}


SAMRAI_INLINE_KEYWORD
bool BoxContainer::insert(const Box& box)
{
   TBOX_ASSERT(box.getLocalId() != LocalId::getInvalidId());
   TBOX_ASSERT(box.getOwnerRank() != tbox::SAMRAI_MPI::getInvalidRank());
   TBOX_ASSERT(box.getBlockId() != BlockId::invalidId());

   if (!d_set_created) {
      TBOX_ASSERT(size() == 0);
      d_set_created = true;
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
