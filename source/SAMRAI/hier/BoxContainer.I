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
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer() :
   d_ordered(false)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(bool ordered) :
   d_ordered(ordered)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const Box& box,
   bool ordered):
   d_list(1, box),
   d_ordered(false)
{
   if (ordered) {
      order();
   }
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const BoxContainer& other):
   d_list(other.d_list),
   d_ordered(false) 
{
   if (other.d_ordered) {
      order();
   }
}

SAMRAI_INLINE_KEYWORD
BoxContainer::~BoxContainer()
{
}

SAMRAI_INLINE_KEYWORD
int BoxContainer::size() const
{
   if (!d_ordered) {
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
   return ConstIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::ConstIterator BoxContainer::end() const
{
   return ConstIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::begin()
{
   return Iterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::end()
{
   return Iterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainer::front() const
{
   return d_ordered ? **(d_set.begin()) : d_list.front();
}

SAMRAI_INLINE_KEYWORD
const Box& BoxContainer::back() const
{
   return d_ordered ? **(d_set.rbegin()) : d_list.back();
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::pushFront(const Box& item)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), item);
   }
#endif
   if (!d_ordered) {
      d_list.push_front(item);
   } else {
      TBOX_ERROR("pushFront on an ordered BoxContainer");
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
   if (!d_ordered) {
      d_list.push_back(item);
   } else {
      TBOX_ERROR("pushBack on an ordered BoxContainer");
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
   if (!d_ordered) {
      d_list.insert(iter.d_list_iter, item);
   } else {
      TBOX_ERROR("insertBefore on an ordered BoxContainer");
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
   if (!d_ordered) {
      d_list.splice(begin().d_list_iter, boxes.d_list);
   } else {
      TBOX_ERROR("spliceFront on an ordered BoxContainer");
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
   if (!d_ordered) {
      boxes.spliceFront(*this);
      d_list.swap(boxes.d_list);
   } else {
      TBOX_ERROR("spliceBack on an ordered BoxContainer");
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popFront()
{
   if (!d_ordered) {
      d_list.pop_front();
   } else {
      TBOX_ERROR("popFront on an ordered BoxContainer");
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popBack()
{
   if (!d_ordered) {
      d_list.pop_back();
   } else {
      TBOX_ERROR("popBack on an ordered BoxContainer");
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::clear()
{
   d_list.clear();
   d_set.clear();
   d_ordered = false; 
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::operator == (const BoxContainer& rhs) const
{
   bool isEqual = (d_ordered && rhs.d_ordered); 
   if (isEqual) {
      if (!d_ordered) { 
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
void BoxContainer::separatePeriodicImages(
   std::vector<Box>& real_mapped_box_vector,
   std::vector<Box>& periodic_image_mapped_box_vector) const
{
   if (!d_ordered) {
      TBOX_ERROR("separatePeriodicImages called on unordered container.");
   }

   if (!isEmpty()) {

      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      real_mapped_box_vector.reserve(real_mapped_box_vector.size() + size());
      for (ConstIterator ni = begin(); ni != end(); ++ni) {
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
   return (d_ordered);
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::insert(const Box& box)
{
   TBOX_ASSERT(box.getId().isValid());

   if (!d_ordered && size() == 0) {
      order();
   }

   if (!d_ordered) {
      TBOX_ERROR("BoxContainer::insert called on unordered container.");
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
