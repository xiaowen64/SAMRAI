/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_UncoveredBoxIterator_C
#define included_hier_UncoveredBoxIterator_C

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/UncoveredBoxIterator.h"

namespace SAMRAI {
namespace hier {

UncoveredBoxIterator::UncoveredBoxIterator(
   const PatchHierarchy* hierarchy,
   bool begin) :
   d_hierarchy(hierarchy),
   d_uncovered_boxes_itr(d_uncovered_boxes.begin()),
   d_uncovered_boxes_itr_end(d_uncovered_boxes.end()),
   d_item(0),
   d_cur_overlapping_level_box(d_overlapping_level_boxes.begin()),
   d_end_overlapping_level_boxes(d_overlapping_level_boxes.end())
{
   TBOX_ASSERT(hierarchy);

   d_finest_level_num = d_hierarchy->getFinestLevelNumber();
   if (begin) {
      // The iterator begins with the coarsest uncovered boxes so start looking
      // for uncovered boxes from level 0.
      d_level_num = -1;
      findNextFinestUncoveredBoxes();
   }
   else {
      // The iterator ends at the end of the boxes on the finest level which,
      // by definition, are uncovered.
      d_level_num = d_finest_level_num;
      d_uncovered_boxes =
         d_hierarchy->getPatchLevel(d_level_num)->getBoxLevel()->getBoxes();
      d_uncovered_boxes_itr = d_uncovered_boxes.end();
      d_uncovered_boxes_itr_end = d_uncovered_boxes.end();
   }
}

UncoveredBoxIterator::UncoveredBoxIterator(
   const UncoveredBoxIterator& other) :
   d_hierarchy(other.d_hierarchy),
   d_level_num(other.d_level_num),
   d_uncovered_boxes(other.d_uncovered_boxes),
   d_uncovered_boxes_itr(other.d_uncovered_boxes_itr),
   d_uncovered_boxes_itr_end(other.d_uncovered_boxes_itr_end),
   d_item(new std::pair<boost::shared_ptr<Patch>, Box>(*other.d_item)),
   d_finest_level_num(other.d_finest_level_num),
   d_overlapping_level_boxes(other.d_overlapping_level_boxes),
   d_cur_overlapping_level_box(other.d_cur_overlapping_level_box),
   d_end_overlapping_level_boxes(other.d_end_overlapping_level_boxes)
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
      d_uncovered_boxes = rhs.d_uncovered_boxes;
      d_uncovered_boxes_itr = rhs.d_uncovered_boxes_itr;
      d_uncovered_boxes_itr_end = rhs.d_uncovered_boxes_itr_end;
      d_item->first = rhs.d_item->first;
      d_item->second = rhs.d_item->second;
      d_finest_level_num = rhs.d_finest_level_num;
      d_overlapping_level_boxes = rhs.d_overlapping_level_boxes;
      d_cur_overlapping_level_box = d_cur_overlapping_level_box;
      d_end_overlapping_level_boxes = d_end_overlapping_level_boxes;
   }
   return *this;
}

const std::pair<boost::shared_ptr<Patch>, Box>&
UncoveredBoxIterator::operator * () const
{
   return *d_item;
}

const std::pair<boost::shared_ptr<Patch>, Box>*
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
   if (result) {
      // Now check if the iterators are at the same point in the level.  If
      // they are both at the end of the level, then they are equal.  If they
      // are both in the middle of the level and on the same box, then they are
      // equal.  Otherwise one is at the end and the other is in the middle so
      // they are not equal.
      if (d_uncovered_boxes_itr == d_uncovered_boxes_itr_end &&
          rhs.d_uncovered_boxes_itr == rhs.d_uncovered_boxes_itr_end) {
      }
      else if (d_uncovered_boxes_itr != d_uncovered_boxes_itr_end &&
               rhs.d_uncovered_boxes_itr != rhs.d_uncovered_boxes_itr_end) {
         result =
            (d_cur_overlapping_level_box->isSpatiallyEqual(
               *rhs.d_cur_overlapping_level_box)) &&
            (d_uncovered_boxes_itr->isSpatiallyEqual(
               *rhs.d_uncovered_boxes_itr));
      }
      else {
         result = false;
      }
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
   ++d_cur_overlapping_level_box;
   if (d_cur_overlapping_level_box != d_end_overlapping_level_boxes) {
      d_item->first =
         d_hierarchy->getPatchLevel(d_level_num)->getPatch(
            d_cur_overlapping_level_box->getBoxId());
   }
   else {
      ++d_uncovered_boxes_itr;
      if (d_uncovered_boxes_itr == d_uncovered_boxes_itr_end) {
         if (d_level_num != d_finest_level_num) {
            findNextFinestUncoveredBoxes();
         }
      }
      else {
         findOverlappedPatch();
      }
   }
}

void
UncoveredBoxIterator::findNextFinestUncoveredBoxes()
{
   // We're through with the current uncovered boxes.
   d_uncovered_boxes.clear();

   // Continue through this loop until a level with uncovered boxes if found or
   // the finest level is reached.
   while (d_uncovered_boxes.isEmpty()) {
      // Move to the next level and get its boxes.  We'll check and see if any
      // of them are not covered by the next finer level.
      ++d_level_num;
      boost::shared_ptr<PatchLevel> this_level =
         d_hierarchy->getPatchLevel(d_level_num);
      d_uncovered_boxes = this_level->getBoxLevel()->getBoxes();

      // If this is the finest level then the boxes are all uncovered so quit.
      if (d_level_num == d_finest_level_num) {
         break;
      }

      // Get the next level and the connector between the current level and the
      // next.  Use a width of 0 to get actual overlaps.
      IntVector width(d_hierarchy->getDim(), 0);
      boost::shared_ptr<PatchLevel> next_level =
         d_hierarchy->getPatchLevel(d_level_num+1);
      const Connector& this_to_next =
         this_level->findConnector(
            *next_level,
            width,
            CONNECTOR_CREATE,
            true);
      d_uncovered_boxes.unorder();

      // Gather all overlapping boxes into a BoxContainer.
      BoxContainer overlaps;
      for (BoxContainer::iterator ln_itr(d_uncovered_boxes.begin());
           ln_itr != d_uncovered_boxes.end(); ++ln_itr) {
         const BoxId& ln_box_id = ln_itr->getBoxId();
         if (this_to_next.hasNeighborSet(ln_box_id)) {
            Connector::ConstNeighborhoodIterator ln_nbrhd(
               this_to_next.findLocal(ln_box_id));
            for (Connector::ConstNeighborIterator overlap_itr(this_to_next.begin(ln_nbrhd));
                 overlap_itr != this_to_next.end(ln_nbrhd); ++overlap_itr) {
               overlaps.pushBack(*overlap_itr);
            }
         }
      }

      // Coarsen the overlapping Boxes to put them into the index space of the
      // current level.  Then remove the intersections of these overlaps from
      // the current level's boxes.
      overlaps.coarsen(next_level->getRatioToCoarserLevel());
      d_uncovered_boxes.removeIntersections(overlaps);

      // Clear overlaps for next pass through loop.
      overlaps.clear();
   }

   // Set the internal iterator over the uncovered boxes and find the Patch
   // overlapped by the box this iterator points to.
   d_uncovered_boxes_itr = d_uncovered_boxes.begin();
   d_uncovered_boxes_itr_end = d_uncovered_boxes.end();
   d_level_boxes =
      d_hierarchy->getPatchLevel(d_level_num)->getBoxLevel()->getBoxes();
   d_level_boxes.makeTree(&(*d_hierarchy->getGridGeometry()));
   findOverlappedPatch();
}

void
UncoveredBoxIterator::findOverlappedPatch()
{
   // Get rid of the overlapping boxes from the previous uncovered box.
   d_overlapping_level_boxes.clear();

   // Get the current uncovered box and the boxes from the level it came from
   // and find the overlap.
   Box cur_box = *d_uncovered_boxes_itr;
   boost::shared_ptr<PatchLevel> this_level =
      d_hierarchy->getPatchLevel(d_level_num);
   d_level_boxes.findOverlapBoxes(d_overlapping_level_boxes, cur_box);

   // At least one of the level boxes must overlap the current uncovered box.
   size_t num_overlapping_level_boxes = d_overlapping_level_boxes.size();
   TBOX_ASSERT(num_overlapping_level_boxes >= 1);

   // If only one of the level boxes overlapped the current box then we're
   // done.  Otherwise we have overlapping patches and need to do more work.
   if (num_overlapping_level_boxes > 1) {
      // We want the patches associated with any overlapping level box that
      // contains the current box.
      BoxContainer::iterator first(d_overlapping_level_boxes.begin());
      for (BoxContainer::iterator itr(d_overlapping_level_boxes.begin());
           itr != d_overlapping_level_boxes.end(); ++itr) {
         if (itr->contains(cur_box)) {
            if (first != itr) {
               d_overlapping_level_boxes.erase(first, itr);
               first = itr;
            }
            ++first;
         }
      }
      if (first != d_overlapping_level_boxes.end()) {
         d_overlapping_level_boxes.erase(first,
            d_overlapping_level_boxes.end());
      }

      // We will defer reporting overlaps with cur_box.  To do this find any
      // uncovered box yet to be processed that intersects or equals cur_box.
      // For each of these boxes that intersects cur_box remove the
      // intersection with cur_box.  If there is nothing left then cur_box
      // should be skipped.  Otherwise make cur_box the first box left when
      // the intersection is removed and place any other boxes left at the end
      // of d_uncovered_boxes.  If the box equals cur_box skip cur_box.
      bool skip_box = false;
      BoxContainer::const_iterator future = d_uncovered_boxes_itr;
      ++future;
      for (; future != d_uncovered_boxes_itr_end; ++future) {
         if (future->isSpatiallyEqual(cur_box)) {
            skip_box = true;
            break;
         }
         else if (future->intersects(cur_box)) {
            BoxContainer cur(cur_box);
            BoxContainer fut(*future);
            cur.removeIntersections(fut);
            if (cur.isEmpty()) {
               skip_box = true;
               break;
            }
            else if (cur.size() == 1) {
               cur_box = cur.front();
            }
            else {
               cur_box = cur.front();
               cur.popFront();
               for (BoxContainer::const_iterator cur_itr(cur.begin());
                    cur_itr != cur.end(); ++cur_itr) {
                  d_uncovered_boxes.pushBack(*cur_itr);
               }
            }
         }
      }

      // If the box should be skipped increment the iterator to move to the
      // next uncovered box.  Otherwise, find the overlapping level boxes for
      // what's left of cur_box.
      if (skip_box) {
         // Move d_cur_overlapping_level_box to just before the end so that
         // when incrementIterator bumps it up it is at the end and we skip
         // cur_box.
         BoxContainer::const_iterator next(d_cur_overlapping_level_box);
         ++next;
         while (true) {
            if (next == d_end_overlapping_level_boxes) {
               break;
            }
            else {
               ++d_cur_overlapping_level_box;
               next = d_cur_overlapping_level_box;
               ++next;
            }
         }
         incrementIterator();
         return;
      }
      else {
         d_overlapping_level_boxes.clear();
         d_level_boxes.findOverlapBoxes(d_overlapping_level_boxes, cur_box);
      }
   }
   d_cur_overlapping_level_box = d_overlapping_level_boxes.begin();
   d_end_overlapping_level_boxes = d_overlapping_level_boxes.end();

   // Update item with the current originating patch and the current box.
   if (d_item) {
      d_item->first =
         this_level->getPatch(d_cur_overlapping_level_box->getBoxId());
      d_item->second = cur_box;
   }
   else {
      d_item =
         new std::pair<boost::shared_ptr<Patch>, Box>(
            this_level->getPatch(d_cur_overlapping_level_box->getBoxId()),
            cur_box);
   }
}

}
}

#endif
