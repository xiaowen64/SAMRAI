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

/*
 *************************************************************************
 * Constructors and destructor
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer() :
   d_ordered(false)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(const bool ordered) :
   d_ordered(ordered)
{
}

SAMRAI_INLINE_KEYWORD
BoxContainer::BoxContainer(
   const Box& box,
   const bool ordered):
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

/*
 *************************************************************************
 * Size queries
 *************************************************************************
 */

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

/*
 *************************************************************************
 * Access methods
 *************************************************************************
 */

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

/*
 *************************************************************************
 * List-like semantics to manipulate unordered containers
 *************************************************************************
 */

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
      TBOX_ERROR("Attempted pushFront on an ordered BoxContainer");
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
      TBOX_ERROR("Attempted pushBack on an ordered BoxContainer");
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
      TBOX_ERROR("Attempted insertBefore on an ordered BoxContainer");
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::insertAfter(Iterator iter, const Box& item)
{
   if (!d_ordered) {
      Iterator tmp = iter;
      ++tmp;
      if (tmp == end()) {
         pushBack(item);
      } else {
         insertBefore(tmp, item);
      }
   } else {
      TBOX_ERROR("Attempted insertAfter called on ordered BoxContainer.");
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
      TBOX_ERROR("Attempted spliceFront on an ordered BoxContainer");
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
      TBOX_ERROR("Attempted spliceBack on an ordered BoxContainer");
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popFront()
{
   if (!d_ordered) {
      d_list.pop_front();
   } else {
      TBOX_ERROR("Attempted popFront on an ordered BoxContainer");
   }
}

SAMRAI_INLINE_KEYWORD
void BoxContainer::popBack()
{
   if (!d_ordered) {
      d_list.pop_back();
   } else {
      TBOX_ERROR("Attempted popBack on an ordered BoxContainer");
   }
}

/*
 *************************************************************************
 * Clear the container
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void BoxContainer::clear()
{
   d_list.clear();
   d_set.clear();
   d_ordered = false; 
}

/*
 *************************************************************************
 * Query the ordered/unordered state
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
bool BoxContainer::isOrdered() const
{
   return (d_ordered);
}


/*
 *************************************************************************
 * Swap all data
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void BoxContainer::swap(BoxContainer& other)
{
   d_list.swap(other.d_list);
   d_set.swap(other.d_set);

   bool other_set_created = other.d_ordered;
   other.d_ordered = d_ordered;
   d_ordered = other_set_created;
}


/*
 *************************************************************************
 * Set-like semantics for ordered containers
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::find(const Box& box) const
{
   if (!d_ordered) {
      TBOX_ERROR("find attempted on unordered BoxContainer.");
   }

   Iterator iter;
   iter.d_set_iter = d_set.find(const_cast<Box*>(&box));
   iter.d_ordered = true;

   return iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::lowerBound(const Box& box) const
{
   if (!d_ordered) {
      TBOX_ERROR("lowerBound attempted on unordered BoxContainer.");
   }

   Iterator iter;
   iter.d_set_iter = d_set.lower_bound(const_cast<Box*>(&box));
   iter.d_ordered = true;

   return iter;
}

SAMRAI_INLINE_KEYWORD
BoxContainer::Iterator BoxContainer::upperBound(const Box& box) const
{
   if (!d_ordered) {
      TBOX_ERROR("upperBound attempted on unordered BoxContainer.");
   }

   Iterator iter;
   iter.d_set_iter = d_set.upper_bound(const_cast<Box*>(&box));
   iter.d_ordered = true;

   return iter;
}

SAMRAI_INLINE_KEYWORD
bool BoxContainer::insert(const Box& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(box.getId().isValid());
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), box);
   }
#endif

   if (!d_ordered && size() == 0) {
      order();
   }

   if (!d_ordered) {
      TBOX_ERROR("BoxContainer::insert called on unordered BoxContainer.");
   }

   const std::list<Box>::iterator& iter = d_list.insert(d_list.end(), box);
   Box* box_ptr(&(*iter));
   if (d_set.insert(box_ptr).second) {
      box_ptr->lockId();
      return true;
   } else {
      d_list.erase(iter);
      return false;
   }

}

/*
 *************************************************************************
 * Equality/inequality
 *************************************************************************
 */

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




}
}
