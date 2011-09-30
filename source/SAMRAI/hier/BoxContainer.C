/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_BoxContainer_C
#define included_hier_BoxContainer_C

#include "SAMRAI/hier/BoxContainer.h"

#include "SAMRAI/hier/BoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/Index.h"
#ifdef MB_BOXTREE_EXISTS
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/MultiblockBoxTree.h"
#else
#include "SAMRAI/hier/BoxTree.h"
#endif

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainer.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int BoxContainer::HIER_BOX_CONTAINER_VERSION = 0;


BoxContainer::BoxContainer(
   Iterator first,
   Iterator last):
   d_dim(first().getDim()),
   d_list(),
   d_set_created(false)
{
   while (first != last) {
      d_list.push_back(first());
      ++first;
   }
}

BoxContainer& BoxContainer::operator = (
   const BoxContainer& rhs)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, rhs);
   clear();

   if (this != &rhs) {
      d_list = rhs.d_list;
      if (rhs.d_set_created) { 
         makeSet();
      } else {
         d_set_created = false;
      }
   }
   return *this;
}

BoxContainer& BoxContainer::operator = (
   const tbox::Array<tbox::DatabaseBox>& rhs)
{
   clear();

   const int n = rhs.size();
   for (int j = 0; j < n; j++) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, rhs[j]);
      pushBack(Box(rhs[j]));
   }
   d_set_created = false;

   return *this;
}

BoxContainer::BoxContainer(
   const tbox::Array<tbox::DatabaseBox>& other):
   d_dim(other.size() == 0 ? tbox::Dimension::getInvalidDimension() :
         other[0].getDim()),
   d_set_created(false)
{
   const int n = other.size();
   for (int j = 0; j < n; j++) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, other[j]);
      pushBack(Box(other[j]));
   }
}

void BoxContainer::insertAfter(Iterator iter, const Box& item)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, item);
   Iterator tmp = iter;
   ++tmp;
   if (tmp == end()) {
      pushBack(item);
   } else {
      insertBefore(tmp, item);
   }
}

/*************************************************************************
 *
 * Function simplify() takes the complicated container of boxes and
 * coalesces regions together where possible.
 *
 * The canonical ordering for boxes is defined such that boxes which
 * lie next to each other in higher dimensions are coalesced together
 * before boxes which lie next to each other in lower dimensions.
 * Thus, we try to coalesce two boxes together on the higher
 * dimensions first.
 *
 * Assuming that two boxes a and b of dimension DIM are in canonical
 * order for dimensions d+1, ..., D, we can coalesce them together on
 * dimension d if:
 *
 *      (1) the lower and upper bounds for a and b agree for all
 *          dimensions greater than d
 *      (2) boxes a and b overlap or are next to each other in
 *          dimension d
 *      (3) boxes a and b overlap for all dimensions less than d
 *
 * If these conditions hold, then we break up the two boxes and put
 * them into the container of non-canonical boxes.
 *
 *************************************************************************
 */
void BoxContainer::simplify()
{
   // Start coalescing on the highest dimension of the containers and work down
   // While there are non-canonical boxes, pick somebody out of the container.

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   if (!isEmpty()) {

      BoxContainer notCanonical(d_dim);
      for (int d = d_dim.getValue() - 1; d >= 0; d--) {
         notCanonical.spliceBack(*this);
         while (!notCanonical.isEmpty()) {
            Box tryMe = notCanonical.front();
            notCanonical.popFront();

            // Pick somebody off of the canonical container and compare
            // against tryMe.

            if (!tryMe.empty()) {
               bool combineDaPuppies = false;
               Iterator l(*this);
               for ( ; l; l++) {
                  const Box andMe = l();

                  const Index& al = andMe.lower();
                  const Index& ah = andMe.upper();
                  const Index& bl = tryMe.lower();
                  const Index& bh = tryMe.upper();

                  combineDaPuppies = true;
                  for (int du = d + 1; du < d_dim.getValue(); du++) {
                     if ((al(du) != bl(du)) || (ah(du) != bh(du))) {
                        combineDaPuppies = false;
                        break;
                     }
                  }
                  if (combineDaPuppies) {
                     if ((bl(d) > ah(d) + 1) || (bh(d) < al(d) - 1)) {
                        combineDaPuppies = false;
                     } else {
                        for (int dl = 0; dl < d; dl++) {
                           if ((bl(dl) > ah(dl)) || (bh(dl) < al(dl))) {
                              combineDaPuppies = false;
                              break;
                           }
                        }
                     }
                  }
                  if (combineDaPuppies) {
                     break;
                  }
               }

               // If we are at the end of the canonical container, then just
               // add.  Otherwise, burst tryMe and andMe and put on
               // noncanonical.

               if (!combineDaPuppies) {
                  pushBack(tryMe);
               } else {
                  Box andMe = l();
                  erase(l);
                  const Index& bl = tryMe.lower();
                  const Index& bh = tryMe.upper();
                  Index il = andMe.lower();
                  Index ih = andMe.upper();
                  for (int dl = 0; dl < d; dl++) {
                     if (il(dl) < bl(dl)) {
                        il(dl) = bl(dl);
                     }
                     if (ih(dl) > bh(dl)) {
                        ih(dl) = bh(dl);
                     }
                  }
                  if (bl(d) < il(d)) {
                     il(d) = bl(d);
                  }
                  if (bh(d) > ih(d)) {
                     ih(d) = bh(d);
                  }
                  Box intersection(il, ih);
                  notCanonical.pushFront(intersection);
                  if (d > 0) {
                     notCanonical.burstBoxes(tryMe, intersection, d);
                     notCanonical.burstBoxes(andMe, intersection, d);
                  }
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Coalesce boxes in the container where possible.  The resulting box
 * container will contain a non-overlapping set of boxes covering the
 * identical region of index space covered by the original container.
 * Two boxes may be coalesced if their union is a box (recall that union
 * is not closed over boxes), and they have a non-empty intersection or
 * they are adjacent to each other in index space.  Empty boxes in the
 * container are removed during this process.  Also, the boxes are
 * coalesced in the order in which they appear in the container.  No
 * attempt is made to coalesce boxes in any particular way (e.g., to
 * achieve the smallest number of boxes).
 *
 *************************************************************************
 */
void BoxContainer::coalesce()
{
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   Iterator tb(*this);
   while (tb) {

      bool found_match = false;

      Iterator tb2 = tb;
      tb2++;

      while (!found_match && tb2) {

         if (tb2().coalesceWith(tb())) {
            found_match = true;
            erase(tb);
         }

         tb2++;
      }

      if (found_match) {
         tb = begin();
      } else {
         tb++;
      }
   }
}

void BoxContainer::grow(
   const IntVector& ghosts)
{
   for (Iterator i(*this); i; i++) {
      i().grow(ghosts);
   }
}

void BoxContainer::shift(
   const IntVector& offset)
{
   for (Iterator i(*this); i; i++) {
      i().shift(offset);
   }
}

void BoxContainer::refine(
   const IntVector& ratio)
{
   for (Iterator i(*this); i; i++) {
      i().refine(ratio);
   }
}

void BoxContainer::coarsen(
   const IntVector& ratio)
{
   for (Iterator i(*this); i; i++) {
      i().coarsen(ratio);
   }
}

void BoxContainer::makeSet()
{
   d_set.clear();
   for (Iterator i(*this); i; i++) {
      d_set.insert(&(i()));
   }
   d_set_created = true;
}

BoxContainerSetIterator
BoxContainer::insert(BoxContainerSetIterator position,
       const Box& box)
{
   if (!d_set_created) {
      makeSet();
   }
   const std::list<Box>::iterator& iter = d_list.insert(d_list.end(), box);

   BoxContainerSetIterator insert_iter(*this);
   int old_size = d_set.size();
   insert_iter.d_set_iter = d_set.insert(position.d_set_iter, &(*iter));
   if (d_set.size() == old_size) {
      d_list.erase(iter);
   }
   return insert_iter;
}

void BoxContainer::insert ( BoxContainerSetConstIterator first,
                            BoxContainerSetConstIterator last )
{
   if (!d_set_created) {
      makeSet();
   }

   for (std::set<const Box*>::const_iterator set_iter = first.d_set_iter;
        set_iter != last.d_set_iter; ++set_iter) {
      const std::list<Box>::iterator& iter = d_list.insert(d_list.end(), **set_iter);
      if (!d_set.insert(&(*iter)).second) {
         d_list.erase(iter); 
      }
   }

}

BoxContainerSetIterator BoxContainer::find(const Box& box) const
{
   BoxContainerSetIterator iter(*this);
   iter.d_set_iter = d_set.find(&box);

   return iter;
}

void
BoxContainer::swap(BoxContainer& other)
{
   d_list.swap(other.d_list);
   d_set.swap(other.d_set);

   bool other_set_created = other.d_set_created;
   other.d_set_created = d_set_created;
   d_set_created = other_set_created;
}

void
BoxContainer::removePeriodicImageBoxes()
{
   for (SetIterator na = setBegin(); na != setEnd(); ) {
      if (na->isPeriodicImage()) {
         erase(na++);
      }
      else {
         ++na;
      }
   }
}


void BoxContainer::erase(BoxContainerSetIterator iter)
{
   const Box& box = **(iter.d_set_iter);
   d_set.erase(iter.d_set_iter);

   for (std::list<Box>::iterator bi = d_list.begin(); bi != d_list.end(); ++bi) {
      if (bi->getId() == box.getId() && bi->isSpatiallyEqual(box)) {
         d_list.erase(bi);
         break;
      }
   }

}

int BoxContainer::erase(const Box& box)
{
   bool ret = d_set.erase(&box);
   for (std::list<Box>::iterator bi = d_list.begin(); bi != d_list.end(); ++bi) {
      if (bi->getId() == box.getId() && bi->isSpatiallyEqual(box)) {
         d_list.erase(bi++);
         break;
      }
   }

   return ret;
}

#if 0
void BoxContainer::erase(BoxContainerSetIterator first,
                         BoxContainerSetIterator last)
{
   d_set.erase(first.d_set_iter, last.d_set_iter);
   //TODO:: Linear search to erase from d_list? 
}
#endif
BoxContainerSetIterator BoxContainer::lower_bound(const Box& box) const
{
   BoxContainerSetIterator iter(*this);

   iter.d_set_iter = d_set.lower_bound(&box);

   return iter;
}

BoxContainerSetIterator BoxContainer::upper_bound(const Box& box) const
{
   BoxContainerSetIterator iter(*this);

   iter.d_set_iter = d_set.upper_bound(&box);

   return iter;
}

void BoxContainer::rotate(
   const Transformation::RotationIdentifier rotation_ident)
{
   if (getDim().getValue() == 2 || getDim().getValue() == 3) {
      for (Iterator i(*this); i; i++) {
         i().rotate(rotation_ident);
      }
   } else {
      NULL_USE(rotation_ident);

      TBOX_ERROR("BoxContainer::rotate() error ..."
         << "\n   Rotation only implemented for 2D and 3D " << std::endl);
   }
}

int BoxContainer::getTotalSizeOfBoxes() const
{
   int size = 0;
   if (!d_set_created) {
      for (ConstIterator i(*this); i; i++) {
         size += i().size();
      }
   } else {
      for (SetConstIterator i(*this); i; i++) {
         size += i().size(); 
      }
   } 
   return size;
}

bool BoxContainer::contains(
   const Index& idx) const
{
   for (ConstIterator i(*this); i; i++) {
      if (i().contains(idx)) {
         return true;
      }
   }
   return false;
}

/*
 *************************************************************************
 *
 * Return the bounding box for all boxes in the BoxContainer.
 *
 *************************************************************************
 */
Box BoxContainer::getBoundingBox() const
{
   if (isEmpty()) {
      TBOX_WARNING("Bounding box container is empty");
      const tbox::Dimension dim(tbox::Dimension::getInvalidDimension());
      Box empty(dim);
      return empty;
   } else {
      Box bbox(d_dim);
      for (ConstIterator i(*this); i; i++) {
         bbox += i();
      }

      return bbox;
   }
}

/*
 *************************************************************************
 *
 * Test the box container for intersections among its boxes.
 *
 *************************************************************************
 */
bool BoxContainer::boxesIntersect() const
{
   bool intersections = false;

   ConstIterator tryMe(*this);
   ConstIterator whatAboutMe(*this);
   whatAboutMe++;
   while (!intersections && tryMe) {
      while (!intersections && whatAboutMe) {
         if (!((tryMe() * whatAboutMe()).size() == 0)) {
            intersections = true;
         }
         whatAboutMe++;
      }
      tryMe++;
      whatAboutMe = tryMe;
      whatAboutMe++;
   }
   return intersections;
}

/*
 *************************************************************************
 *
 * Return the current container without the portions that intersect
 * takeaway.
 *
 *************************************************************************
 */
void BoxContainer::removeIntersections(
   const Box& takeaway)
{
   if (isEmpty()) {
      return;
   }

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   const unsigned short dim = takeaway.getDim().getValue();
   Iterator insertion_pt(*this);
   while (insertion_pt) {
      Box& tryme = *insertion_pt;
      if (!tryme.intersects(takeaway)) {
         insertion_pt++;
      } else {
         Iterator tmp = insertion_pt;
         burstBoxes(tryme, takeaway, dim, insertion_pt);
         insertion_pt++;
         erase(tmp);
      }
   }
}

void BoxContainer::removeIntersections(
   const BoxContainer& takeaway)
{
   if (isEmpty()) {
      return;
   }

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   for (ConstIterator remove(takeaway); remove; remove++) {
      const Box& byebye = remove();
      removeIntersections(byebye);
   }
}

void BoxContainer::removeIntersections(
   const BoxTree& takeaway)
{
   if (isEmpty()) {
      return;
   }

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   std::vector<const Box *> overlap_mapped_boxes;
   Iterator itr(*this);
   while (itr) {
      const Box& tryme = *itr;
      takeaway.findOverlapBoxes(overlap_mapped_boxes, tryme);
      if (overlap_mapped_boxes.empty()) {
         itr++;
      } else {
         Iterator sublist_start = itr;
         Iterator sublist_end = sublist_start;
         sublist_end++;
         for (size_t i = 0;
              i < overlap_mapped_boxes.size() && sublist_start != sublist_end;
              ++i) {
            Iterator insertion_pt = sublist_start;
            removeIntersectionsFromSublist(
               *overlap_mapped_boxes[i],
               sublist_start,
               sublist_end,
               insertion_pt);
         }
         overlap_mapped_boxes.clear();
         itr = sublist_end;
      }
   }
}

#ifdef MB_BOXTREE_EXISTS
void BoxContainer::removeIntersections(
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   const MultiblockBoxTree& takeaway,
   bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   const tbox::ConstPointer<hier::GridGeometry>
   & grid_geometry(takeaway.getGridGeometry());

   std::vector<const Box *> overlap_mapped_boxes;
   Iterator itr(*this);
   while (itr) {
      const Box& tryme = *itr;
      takeaway.findOverlapBoxes(overlap_mapped_boxes,
         tryme,
         block_id,
         refinement_ratio,
         include_singularity_block_neighbors);
      if (overlap_mapped_boxes.empty()) {
         itr++;
      } else {
         Iterator sublist_start = itr;
         Iterator sublist_end = sublist_start;
         sublist_end++;
         for (size_t i = 0;
              i < overlap_mapped_boxes.size() && sublist_start != sublist_end;
              ++i) {
            Iterator insertion_pt = sublist_start;
            const BlockId& overlap_box_block_id =
               overlap_mapped_boxes[i]->getBlockId();
            if (overlap_box_block_id != block_id) {
               Box overlap_box = *overlap_mapped_boxes[i];
               grid_geometry->transformBox(overlap_box,
                  refinement_ratio,
                  block_id,
                  overlap_box_block_id);
               removeIntersectionsFromSublist(
                  overlap_box,
                  sublist_start,
                  sublist_end,
                  insertion_pt);
            } else {
               removeIntersectionsFromSublist(
                  *overlap_mapped_boxes[i],
                  sublist_start,
                  sublist_end,
                  insertion_pt);
            }
         }
         overlap_mapped_boxes.clear();
         itr = sublist_end;
      }
   }
}
#endif

void BoxContainer::removeIntersections(
   const Box& box,
   const Box& takeaway)
{
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   /*
    * The box container MUST be empty to use this function (see comments
    * in header file for discussion of why). If the two boxes intersect,
    * form a BoxContainer that contains the boxes resulting from removing
    * the intersection of box with takeaway.  If the two boxes do not
    * intersect, simply add box to the box container (no intersection removed).
    */
   TBOX_ASSERT(isEmpty());

   if (box.intersects(takeaway)) {
      burstBoxes(box, takeaway, box.getDim().getValue());
   } else {
      pushBack(box);
   }

}

void BoxContainer::removeIntersectionsFromSublist(
   const Box& takeaway,
   Iterator& sublist_start,
   Iterator& sublist_end,
   Iterator& insertion_pt)
{
   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   const unsigned short dim = takeaway.getDim().getValue();
   Iterator itr = sublist_start;
   while (itr != sublist_end) {
      Box& tryme = *itr;
      if (!tryme.intersects(takeaway)) {
         itr++;
      } else {
         burstBoxes(tryme, takeaway, dim, insertion_pt);
         Iterator tmp = itr;
         itr++;
         if (tmp == sublist_start) {
            sublist_start++;
         }
         erase(tmp);
      }
      insertion_pt = itr;
   }
}

/*
 *************************************************************************
 *
 * Return the boxes in the current container that intersect the index
 * space of the argument.
 *
 *************************************************************************
 */
void BoxContainer::intersectBoxes(
   const Box& keep)
{
   if (isEmpty()) {
      return;
   }

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   Iterator i(*this);
   Box overlap(i().getDim());
   while (i) {
      Box& tryMe = *i;
      tryMe.intersect(keep, overlap);
      if (!overlap.empty()) {
         tryMe = overlap;
         i++;
      } else {
         Iterator tmp = i;
         i++;
         erase(tmp);
      }
   }
}

void BoxContainer::intersectBoxes(
   const BoxContainer& keep)
{
   if (isEmpty()) {
      return;
   }

   if (d_set_created) {
      d_set.clear();
      d_set_created = false;
   }

   Iterator insertion_pt(*this);
   Box overlap(insertion_pt().getDim());
   while (insertion_pt) {
      Iterator tmp = insertion_pt;
      const Box& tryme = *insertion_pt;
      for (ConstIterator i(keep); i; i++) {
         tryme.intersect(i(), overlap);
         if (!overlap.empty()) {
            insertAfter(insertion_pt, overlap);
            insertion_pt++;
         }
      }
      insertion_pt++;
      erase(tmp);
   }
}

void BoxContainer::intersectBoxes(
   const BoxTree& keep)
{
   if (isEmpty()) {
      return;
   }

   std::vector<const Box *> overlap_mapped_boxes;
   Box overlap(front().getDim());
   Iterator itr(*this);
   Iterator insertion_pt = itr;
   while (itr) {
      const Box& tryme = *itr;
      keep.findOverlapBoxes(overlap_mapped_boxes, tryme);
      for (size_t i = 0; i < overlap_mapped_boxes.size(); ++i) {
         tryme.intersect(*overlap_mapped_boxes[i], overlap);
         if (!overlap.empty()) {
            insertAfter(insertion_pt, overlap);
            insertion_pt++;
         }
      }
      overlap_mapped_boxes.clear();
      Iterator tmp = itr;
      insertion_pt++;
      itr = insertion_pt;
      erase(tmp);
   }
}

#ifdef MB_BOXTREE_EXISTS
void BoxContainer::intersectBoxes(
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   const MultiblockBoxTree& keep,
   bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   const tbox::ConstPointer<hier::GridGeometry>
   & grid_geometry(keep.getGridGeometry());

   std::vector<const Box *> overlap_mapped_boxes;
   Box overlap(front().getDim());
   Iterator itr(*this);
   Iterator insertion_pt = itr;
   while (itr) {
      const Box& tryme = *itr;
      keep.findOverlapBoxes(overlap_mapped_boxes,
         tryme,
         block_id,
         refinement_ratio,
         include_singularity_block_neighbors);
      for (size_t i = 0; i < overlap_mapped_boxes.size(); ++i) {
         const BlockId& overlap_box_block_id =
            overlap_mapped_boxes[i]->getBlockId();
         if (overlap_box_block_id != block_id) {
            Box overlap_box = *overlap_mapped_boxes[i];
            grid_geometry->transformBox(overlap_box,
               refinement_ratio,
               block_id,
               overlap_box_block_id);
            tryme.intersect(overlap_box, overlap);
            if (!overlap.empty()) {
               insertAfter(insertion_pt, overlap);
               insertion_pt++;
            }
         } else {
            tryme.intersect(*overlap_mapped_boxes[i], overlap);
            if (!overlap.empty()) {
               insertAfter(insertion_pt, overlap);
               insertion_pt++;
            }
         }
      }
      overlap_mapped_boxes.clear();
      Iterator tmp = itr;
      insertion_pt++;
      itr = insertion_pt;
      erase(tmp);
   }
}
#endif

/*
 *************************************************************************
 *
 * Type conversion from a BoxContainer to an Array of
 * tbox::DatabaseBoxes.
 *
 *************************************************************************
 */
BoxContainer::operator tbox::Array<tbox::DatabaseBox>() const
{
   tbox::Array<tbox::DatabaseBox> new_Array(size());

   int j = 0;
   for (ConstIterator i(*this); i; i++) {
      new_Array[j++] = (tbox::DatabaseBox)(*i);
   }

   return new_Array;
}

/*
 *************************************************************************
 *
 * Print boxes in container.
 *
 *************************************************************************
 */
void BoxContainer::print(
   std::ostream& os) const
{
   int i = 0;
   for (ConstIterator b(*this); b; b++) {
      os << "Box # " << i << ":  " << b() << std::endl;
      i++;
   }
}

/*
 *************************************************************************
 *
 * Break up box bursty against box solid and adds the pieces to
 * container.  The bursting is done on dimensions 0 through dimension-1,
 * starting with lowest dimensions first to try to maintain the
 * canonical representation for the bursted domains.
 *
 *************************************************************************
 */
void BoxContainer::burstBoxes(
   const Box& bursty,
   const Box& solid,
   const int dimension)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(bursty, solid);
   TBOX_ASSERT(dimension <= bursty.getDim().getValue());

   // Set up the lower and upper bounds of the regions for ease of access

   Index burstl = bursty.lower();
   Index bursth = bursty.upper();
   const Index& solidl = solid.lower();
   const Index& solidh = solid.upper();

   // Break bursty region against solid region along low dimensions first

   for (int d = 0; d < dimension; d++) {
      if (bursth(d) > solidh(d)) {
         Index newl = burstl;
         newl(d) = solidh(d) + 1;
         pushBack(Box(newl, bursth));
         bursth(d) = solidh(d);
      }
      if (burstl(d) < solidl(d)) {
         Index newh = bursth;
         newh(d) = solidl(d) - 1;
         pushBack(Box(burstl, newh));
         burstl(d) = solidl(d);
      }
   }
}

/*
 *************************************************************************
 *
 * Break up box bursty against box solid and adds the pieces to
 * container starting at the location pointed to by the supplied
 * iterator.  The bursting is done on dimensions 0 through dimension-1,
 * starting with lowest dimensions first to try to maintain the
 * canonical representation for the bursted domains.
 *
 *************************************************************************
 */
void BoxContainer::burstBoxes(
   const Box& bursty,
   const Box& solid,
   const int dimension,
   Iterator& insertion_pt)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(bursty, solid);
   TBOX_ASSERT(dimension <= bursty.getDim().getValue());

   // Set up the lower and upper bounds of the regions for ease of access

   Index burstl = bursty.lower();
   Index bursth = bursty.upper();
   const Index& solidl = solid.lower();
   const Index& solidh = solid.upper();

   // Break bursty region against solid region along low dimensions first

   for (int d = 0; d < dimension; d++) {
      if (bursth(d) > solidh(d)) {
         Index newl = burstl;
         newl(d) = solidh(d) + 1;
         insertAfter(insertion_pt, Box(newl, bursth));
         bursth(d) = solidh(d);
         insertion_pt++;
      }
      if (burstl(d) < solidl(d)) {
         Index newh = bursth;
         newh(d) = solidl(d) - 1;
         insertAfter(insertion_pt, Box(burstl, newh));
         burstl(d) = solidl(d);
         insertion_pt++;
      }
   }
}

/*
 ***********************************************************************
 * Construct a BoxContainer consisting of the Boxes in this BoxContainer
 * in the requested block.
 ***********************************************************************
 */
tbox::Pointer<BoxContainer>
BoxContainer::getSingleBlockBoxContainer(
   const tbox::Dimension& dim,
   const BlockId& which_block) const
{
   BoxSetSingleBlockIterator itr(*this, which_block);
   BoxContainer* boxes_in_block = new BoxContainer(dim);
   while (itr.isValid()) {
      const Box& mapped_box = *itr;
      TBOX_ASSERT(dim == mapped_box.getDim());
      boxes_in_block->pushBack(mapped_box);
      ++itr;
   }
   return tbox::Pointer<BoxContainer>(boxes_in_block);
}

/*
 ***********************************************************************
 * Insert Box owners into a single set container.
 ***********************************************************************
 */
void BoxContainer::getOwners(
   std::set<int>& owners) const
{
   for (BoxContainerSetConstIterator i_nabr = setBegin();
        i_nabr != setEnd(); ++i_nabr) {
      const int owner = (*i_nabr).getOwnerRank();
      owners.insert(owner);
   }
}

/*
 ***********************************************************************
 * Unshift periodic image Boxes from a BoxContainer.
 ***********************************************************************
 */
void BoxContainer::unshiftPeriodicImageBoxes(
   BoxContainer& output_mapped_boxes,
   const IntVector& refinement_ratio) const
{
   BoxContainerSetIterator hint = output_mapped_boxes.setBegin();

   if (!isEmpty()) {
      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      for (BoxContainerSetConstIterator na = setBegin(); na != setEnd(); ++na) {
         if (na->isPeriodicImage()) {
            const Box unshifted_mapped_box(
               *na, zero_shift_number, refinement_ratio);
            hint = output_mapped_boxes.insert(hint, unshifted_mapped_box);
         } else {
            hint = output_mapped_boxes.insert(hint, *na);
         }
      }
   }
}

/*
 ***********************************************************************
 * Write the BoxContainer to a database.
 ***********************************************************************
 */
void BoxContainer::putToDatabase(
   tbox::Database& database) const
{
   database.putInteger(
      "HIER_BOX_CONTAINER_VERSION", HIER_BOX_CONTAINER_VERSION);

   const int mbs_size = size();
   database.putInteger("mapped_box_set_size", mbs_size);
   if (mbs_size > 0) {

      std::vector<int> local_ids;
      std::vector<int> ranks;
      std::vector<int> block_ids;
      std::vector<int> periodic_ids;
      local_ids.reserve(mbs_size);
      ranks.reserve(mbs_size);
      block_ids.reserve(mbs_size);
      periodic_ids.reserve(mbs_size);

      tbox::Array<tbox::DatabaseBox> db_box_array(mbs_size);

      int counter = -1;
      for (BoxContainer::SetConstIterator ni = setBegin();
           ni != setEnd(); ++ni) {
         local_ids.push_back(ni->getLocalId().getValue());
         ranks.push_back(ni->getOwnerRank());
         block_ids.push_back(ni->getBlockId().getBlockValue());
         periodic_ids.push_back(ni->getPeriodicId().getPeriodicValue());
         db_box_array[++counter] = *ni;
      }

      database.putIntegerArray(
         "local_indices", &local_ids[0], mbs_size);
      database.putIntegerArray(
         "ranks", &ranks[0], mbs_size);
      database.putIntegerArray(
         "block_ids", &block_ids[0], mbs_size);
      database.putIntegerArray(
         "periodic_ids", &periodic_ids[0], mbs_size);
      database.putDatabaseBoxArray(
         "boxes", &db_box_array[0], mbs_size);
   }
}


/*
 ***********************************************************************
 * Read the BoxContainer from a database.
 ***********************************************************************
 */
void BoxContainer::getFromDatabase(
   tbox::Database& database)
{
   const unsigned int mbs_size = database.getInteger("mapped_box_set_size");
   if (mbs_size > 0) {
      std::vector<int> local_ids(mbs_size);
      std::vector<int> ranks(mbs_size);
      std::vector<int> block_ids(mbs_size);
      std::vector<int> periodic_ids(mbs_size);
      tbox::Array<tbox::DatabaseBox> db_box_array(mbs_size);

      database.getIntegerArray(
         "local_indices", &local_ids[0], mbs_size);
      database.getIntegerArray(
         "ranks", &ranks[0], mbs_size);
      database.getIntegerArray(
         "block_ids", &block_ids[0], mbs_size);
      database.getIntegerArray(
         "periodic_ids", &periodic_ids[0], mbs_size);
      database.getDatabaseBoxArray(
         "boxes", &db_box_array[0], mbs_size);

      for (unsigned int i = 0; i < mbs_size; ++i) {
         Box box(db_box_array[i]);
         Box mapped_box(
            box,
            LocalId(local_ids[i]),
            ranks[i],
            BlockId(block_ids[i]),
            PeriodicId(periodic_ids[i]));
         insert(setEnd(), mapped_box);
      }
   }
}


/*
 ***********************************************************************
 * Avoid communication in this method.  It is often used for debugging.
 ***********************************************************************
 */
void BoxContainer::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   NULL_USE(detail_depth);
   for (ConstIterator bi = begin(); bi != end(); ++bi) {
      Box mapped_box(*bi);
      co << border << "    "
         << mapped_box << "   "
         << mapped_box.numberCells() << '\n';
   }
}

/*
 ***********************************************************************
 * Construct a BoxContainer Outputter with formatting parameters.
 ***********************************************************************
 */

BoxContainer::Outputter::Outputter(
   const BoxContainer& mapped_box_set,
   const std::string& border,
   int detail_depth):
   d_set(mapped_box_set),
   d_border(border),
   d_detail_depth(detail_depth)
{
}

/*
 ***********************************************************************
 * Print out a BoxContainer according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream& operator << (
   std::ostream& s,
   const BoxContainer::Outputter& format)
{
   format.d_set.recursivePrint(s, format.d_border, format.d_detail_depth);
   return s;
}


/*
 ***********************************************************************
 * Return a Outputter that can dump the BoxContainer to a stream.
 ***********************************************************************
 */

BoxContainer::Outputter BoxContainer::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
}




}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
