/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2015 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/
#include "SAMRAI/hier/UncoveredBoxIterator.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/FlattenedHierarchy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

namespace SAMRAI {
namespace hier {

UncoveredBoxIterator::UncoveredBoxIterator(
   const PatchHierarchy* hierarchy,
   bool begin):
   d_hierarchy(hierarchy),
   d_uncovered_boxes_itr(BoxContainer().begin()),
   d_uncovered_boxes_itr_end(BoxContainer().end()),
   d_item(0)
{
   TBOX_ASSERT(hierarchy);

   d_finest_level_num = d_hierarchy->getFinestLevelNumber();
   if (begin) {
      // The iterator begins with the coarsest uncovered boxes so start looking
      // for uncovered boxes from level 0.
      d_level_num = -1;
      d_flattened_hierarchy = boost::make_shared<FlattenedHierarchy>(*d_hierarchy, 0, d_finest_level_num);
      findNextFinestUncoveredBoxes();
   } else {
      // The iterator ends at the end of the boxes on the finest level which,
      // by definition, are uncovered.
      d_level_num = d_finest_level_num + 1;
      d_flattened_hierarchy.reset();
   }
}

UncoveredBoxIterator::UncoveredBoxIterator(
   const UncoveredBoxIterator& other):
   d_hierarchy(other.d_hierarchy),
   d_level_num(other.d_level_num),
   d_uncovered_boxes_itr(other.d_uncovered_boxes_itr),
   d_uncovered_boxes_itr_end(other.d_uncovered_boxes_itr_end),
   d_item(new std::pair<boost::shared_ptr<Patch>, Box>(*other.d_item)),
   d_finest_level_num(other.d_finest_level_num),
   d_flattened_hierarchy(other.d_flattened_hierarchy)
{
}

UncoveredBoxIterator::~UncoveredBoxIterator()
{
   if (d_item) {
      delete d_item;
   }
}

UncoveredBoxIterator&
UncoveredBoxIterator::operator = (
   const UncoveredBoxIterator& rhs)
{
   if (this != &rhs) {
      d_hierarchy = rhs.d_hierarchy;
      d_level_num = rhs.d_level_num;
      d_uncovered_boxes_itr = rhs.d_uncovered_boxes_itr;
      d_uncovered_boxes_itr_end = rhs.d_uncovered_boxes_itr_end;
      d_item->first = rhs.d_item->first;
      d_item->second = rhs.d_item->second;
      d_finest_level_num = rhs.d_finest_level_num;
      d_flattened_hierarchy = rhs.d_flattened_hierarchy;
   }
   return *this;
}

const std::pair<boost::shared_ptr<Patch>, Box>&
UncoveredBoxIterator::operator * () const
{
   return *d_item;
}

const std::pair<boost::shared_ptr<Patch>, Box> *
UncoveredBoxIterator::operator -> () const
{
   return d_item;
}

bool
UncoveredBoxIterator::operator == (
   const UncoveredBoxIterator& rhs) const
{
   // Frist check and see if the iterators are working on the same hierarchies
   // and levels.  If not then they are not equal.
   bool result = d_hierarchy == rhs.d_hierarchy &&
      d_level_num == rhs.d_level_num;

   if (d_flattened_hierarchy == 0 && rhs.d_flattened_hierarchy != 0) {
      result = false;
   }
   if (d_flattened_hierarchy != 0 && rhs.d_flattened_hierarchy == 0) {
      result = false;
   }
   if (d_item == 0 && rhs.d_item != 0) {
      result = false;
   }
   if (d_item != 0 && rhs.d_item == 0) {
      result = false;
   }
   if (result) {
      if (d_item == 0 && rhs.d_item == 0) {
         result = true;
      }
      if (d_item && rhs.d_item) {
         if (d_item->second.isIdEqual(rhs.d_item->second) &&
             d_item->second.isSpatiallyEqual(rhs.d_item->second) &&
             d_item->first->getBox().isIdEqual(rhs.d_item->first->getBox()) &&
             d_item->first->getBox().isSpatiallyEqual(rhs.d_item->first->getBox())) {
            result = true;
         } else {
            result = false;
      }
   }

#if 0
      // Now check if the iterators are at the same point in the level.  If
      // they are both at the end of the level, then they are equal.  If they
      // are both in the middle of the level and on the same box, then they are
      // equal.  Otherwise one is at the end and the other is in the middle so
      // they are not equal.
      if (d_uncovered_boxes_itr == d_uncovered_boxes_itr_end &&
          rhs.d_uncovered_boxes_itr == rhs.d_uncovered_boxes_itr_end) {
      } else if (d_uncovered_boxes_itr != d_uncovered_boxes_itr_end &&
                 rhs.d_uncovered_boxes_itr != rhs.d_uncovered_boxes_itr_end) {
         result =
            (d_uncovered_boxes_itr->isSpatiallyEqual(
                *rhs.d_uncovered_boxes_itr));
      } else {
         result = false;
      }
#endif
   }
   return result;
}

bool
UncoveredBoxIterator::operator != (
   const UncoveredBoxIterator& rhs) const
{
   return !(*this == rhs);
}

UncoveredBoxIterator&
UncoveredBoxIterator::operator ++ ()
{
   incrementIterator();
   return *this;
}

UncoveredBoxIterator
UncoveredBoxIterator::operator ++ (
   int)
{
   // Save the state of the iterator.
   UncoveredBoxIterator tmp(*this);

   incrementIterator();

   // Return iterator in original state.
   return tmp;
}

void
UncoveredBoxIterator::incrementIterator()
{
   // Move to the next originating patch associated with the current box.  If
   // there are any left then make the iterator point to that patch.  If not
   // then move to the next uncovered box on the current level.  If there are
   // no more left, then check if this finest level.  If so the iteration is at
   // its end.  Otherwise, look for the next finer level with uncovered boxes.
   // If there are uncovered boxes left on the current level then set the
   // iterator return value to point to this next box and the first patch that
   // it overlaps.
   ++d_uncovered_boxes_itr;
   if (d_uncovered_boxes_itr != d_uncovered_boxes_itr_end) {
      d_current_patch_id = d_item->first->getBox().getBoxId();
      setIteratorItem();
   } else {

      bool id_found = false;

      bool new_level = false;
      while (d_level_num <= d_finest_level_num) {

         boost::shared_ptr<PatchLevel> this_level =
            d_hierarchy->getPatchLevel(d_level_num);

         const BoxContainer& this_level_boxes =
            this_level->getBoxLevel()->getBoxes();

         for (RealBoxConstIterator this_itr = this_level_boxes.realBegin();
              this_itr != this_level_boxes.realEnd(); ++this_itr) {

            if (!new_level &&
                this_itr->getBoxId() <= d_item->first->getBox().getBoxId()) {
               continue;
            }
            const BoxContainer& uncovered_boxes =
               d_flattened_hierarchy->getVisibleBoxes(*this_itr, d_level_num);

            if (!uncovered_boxes.empty()) {
               d_uncovered_boxes_itr = uncovered_boxes.begin();
               d_uncovered_boxes_itr_end = uncovered_boxes.end();
               d_current_patch_id = this_itr->getBoxId();
               id_found = true;
               break;
            }
         }
         if (id_found) {
            break;
         } else {
            ++d_level_num;
            new_level = true;
         }
      }

      if (id_found) {
         setIteratorItem();
      } else {
         d_flattened_hierarchy.reset();
         if (d_item) {
            delete d_item;
            d_item = 0; 
         }
      }
   }
}

void
UncoveredBoxIterator::findNextFinestUncoveredBoxes()
{
   ++d_level_num;
   boost::shared_ptr<PatchLevel> this_level =
      d_hierarchy->getPatchLevel(d_level_num);

   const BoxContainer& this_level_boxes = this_level->getBoxLevel()->getBoxes();

   bool id_found = false;

   while (d_level_num <= d_finest_level_num) { 

      boost::shared_ptr<PatchLevel> this_level =
         d_hierarchy->getPatchLevel(d_level_num);

      const BoxContainer& this_level_boxes =
         this_level->getBoxLevel()->getBoxes();

      for (BoxContainer::const_iterator this_itr = this_level_boxes.begin();
           this_itr != this_level_boxes.end(); ++this_itr) {
         const BoxContainer& uncovered_boxes =
            d_flattened_hierarchy->getVisibleBoxes(*this_itr, d_level_num);

         if (!uncovered_boxes.empty()) {
            d_uncovered_boxes_itr = uncovered_boxes.begin();
            d_uncovered_boxes_itr_end = uncovered_boxes.end();
            d_current_patch_id = this_itr->getBoxId(); 
            id_found = true;
            break;
         }
      }
      if (id_found) {
         break;
      } else {
         ++d_level_num;
      }
   }

   if (id_found) {
      setIteratorItem();
   } else {
      if (d_item) {
         delete d_item;
         d_item = 0;
      }
      d_flattened_hierarchy.reset();
   }
}

void
UncoveredBoxIterator::setIteratorItem()
{
   // Get rid of the overlapping boxes from the previous uncovered box.

   // Get the current uncovered box and the boxes from the level it came from
   // and find the overlap.
   const Box& cur_box = *d_uncovered_boxes_itr;
   boost::shared_ptr<PatchLevel> this_level =
      d_hierarchy->getPatchLevel(d_level_num);

   // Update item with the current originating patch and the current box.
   if (d_item) {
      d_item->first =
         this_level->getPatch(d_current_patch_id);
      d_item->second = cur_box;
   } else {
      d_item =
         new std::pair<boost::shared_ptr<Patch>, Box>(
            this_level->getPatch(d_current_patch_id),
            cur_box);
   }
}

}
}
