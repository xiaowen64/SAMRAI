/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A list of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_BoxList_C
#define included_hier_BoxList_C

#include "SAMRAI/hier/BoxList.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/BoxTree.h"
#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxList.I"
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
#if 0
/*
 *************************************************************************
 *                                                                       *
 * Implement the various constructors and the assignment operator for    *
 * the list of boxes.  Note that none of these routines modify their     *
 * arguments.                                                            *
 *                                                                       *
 *************************************************************************
 */

BoxList::BoxList(
   const Box& box):
   tbox::List<Box>(),
   d_dim(box.getDim())
{
   addItem(box);
}

BoxList::BoxList(
   const BoxList& list):
   tbox::List<Box>(),
   d_dim(list.getDim())
{
   copyItems(list);
}

BoxList::BoxList(
   const tbox::Array<tbox::DatabaseBox>& array):
   d_dim(array[0].getDim())
{
   const int n = array.getSize();

   for (int j = 0; j < n; j++) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, array[j]);
      const Box this_box(array[j]);
      appendItem(this_box);
   }
}

BoxList::operator tbox::Array<tbox::DatabaseBox>() const
{
   int number_boxes = getNumberOfBoxes();
   tbox::Array<tbox::DatabaseBox> new_Array(number_boxes);

   int j = 0;
   for (Iterator i(*this); i; i++) {
      new_Array[j++] = (tbox::DatabaseBox)(*i);
   }

   return new_Array;
}

BoxList& BoxList::operator = (
   const tbox::Array<tbox::DatabaseBox>& array)
{
   const int n = array.getSize();
   clearItems();

   for (int j = 0; j < n; j++) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, array[j]);
      const Box this_box(array[j]);
      appendItem(this_box);
   }
   return *this;
}

BoxList& BoxList::operator = (
   const BoxList& list)
{
//   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, list);
   if (this != &list) {
      clearItems();
      copyItems(list);
   }
   return *this;
}

/*************************************************************************
 *                                                                       *
 * Function simplifyBoxes() takes the complicated list of boxes and      *
 * coalesces regions together where possible.                            *
 *                                                                       *
 * The canonical ordering for boxes is defined such that boxes which     *
 * 1* lie next to each other in higher dimensions are coalesced together *
 * before boxes which lie next to each other in lower dimensions.        *
 * Thus, we try to coalesce two boxes together on the higher             *
 * dimensions first.                                                     *
 *                                                                       *
 * Assuming that two boxes a and b of dimension DIM are in canonical     *
 * order for dimensions d+1, ..., D, we can coalesce them together on    *
 * dimension d if:                                                       *
 *                                                                       *
 *      (1) the lower and upper bounds for a and b agree for all         *
 *          dimensions greater than d                                    *
 *      (2) boxes a and b overlap or are next to each other in           *
 *          dimension d                                                  *
 *      (3) boxes a and b overlap for all dimensions less than d         *
 *                                                                       *
 * If these conditions hold, then we break up the two boxes and put      *
 * them onto the list of non-canonical boxes.                            *
 *                                                                       *
 *************************************************************************
 */

void BoxList::simplifyBoxes()
{
   // Start coalescing on the highest dimension of the lists and work down
   // While there are non-canonical boxes, pick somebody off of the list

   if (!this->isEmpty()) {

      tbox::Dimension dim(this->getFirstItem().getDim());

      BoxList notCanonical;
      for (int d = dim.getValue() - 1; d >= 0; d--) {
         notCanonical.catenateItems(*this);
         while (!notCanonical.isEmpty()) {
            Box tryMe = notCanonical.getFirstItem();
            notCanonical.removeFirstItem();

            // Pick somebody off of the canonical list and compare against tryMe

            if (!tryMe.empty()) {
               bool combineDaPuppies = false;
               Iterator l;
               for (l = this->listStart(); l; l++) {
                  const Box andMe = l();

                  const Index& al = andMe.lower();
                  const Index& ah = andMe.upper();
                  const Index& bl = tryMe.lower();
                  const Index& bh = tryMe.upper();

                  combineDaPuppies = true;
                  for (int du = d + 1; du < dim.getValue(); du++) {
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
                  if (combineDaPuppies) break;
               }

               // If we are at the end of the canonical list, then just add
               // Otherwise, burst tryMe and andMe and put on noncanonical

               if (!combineDaPuppies) {
                  appendItem(tryMe);
               } else {
                  Box andMe = l();
                  removeItem(l);
                  const Index& bl = tryMe.lower();
                  const Index& bh = tryMe.upper();
                  Index il = andMe.lower();
                  Index ih = andMe.upper();
                  for (int dl = 0; dl < d; dl++) {
                     if (il(dl) < bl(dl)) il(dl) = bl(dl);
                     if (ih(dl) > bh(dl)) ih(dl) = bh(dl);
                  }
                  if (bl(d) < il(d)) il(d) = bl(d);
                  if (bh(d) > ih(d)) ih(d) = bh(d);
                  Box intersection(il, ih);
                  notCanonical.addItem(intersection);
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
 *                                                                       *
 * Break up box bursty against box solid and adds the pieces to list.    *
 * The bursting is done on dimensions 0 through dimension-1, starting    *
 * with lowest dimensions first to try to maintain the canonical         *
 * representation for the bursted domains.                               *
 *                                                                       *
 *************************************************************************
 */

void BoxList::burstBoxes(
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
         appendItem(Box(newl, bursth));
         bursth(d) = solidh(d);
      }
      if (burstl(d) < solidl(d)) {
         Index newh = bursth;
         newh(d) = solidl(d) - 1;
         appendItem(Box(burstl, newh));
         burstl(d) = solidl(d);
      }
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Break up box bursty against box solid and adds the pieces to list.    *
 * The bursting is done on dimensions 0 through dimension-1, starting    *
 * with lowest dimensions first to try to maintain the canonical         *
 * representation for the bursted domains.                               *
 *                                                                       *
 *************************************************************************
 */

void BoxList::burstBoxes(
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
         addItemAfter(insertion_pt, Box(newl, bursth));
         bursth(d) = solidh(d);
         insertion_pt++;
      }
      if (burstl(d) < solidl(d)) {
         Index newh = bursth;
         newh(d) = solidl(d) - 1;
         addItemAfter(insertion_pt, Box(burstl, newh));
         burstl(d) = solidl(d);
         insertion_pt++;
      }
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Return the current list without the portions that intersect takeaway. *
 *                                                                       *
 *************************************************************************
 */

void BoxList::removeIntersections(
   const Box& takeaway)
{
   if (isEmpty()) {
      return;
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
         removeItem(tmp);
      }
   }
}

void BoxList::removeIntersectionsFromSublist(
   const Box& takeaway,
   Iterator& sublist_start,
   Iterator& sublist_end,
   Iterator& insertion_pt)
{
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
         removeItem(tmp);
      }
      insertion_pt = itr;
   }
}

void BoxList::removeIntersections(
   const Box& box,
   const Box& takeaway)
{
   /*
    * The boxlist MUST be empty to use this function (see comments
    * in header file for discussion of why). If the two boxes intersect,
    * form a boxlist that contains the boxes resulting from removing
    * the intersection of box with takeaway.  If the two boxes do not
    * intersect, simply add box to the box list (no intersection removed).
    */
   TBOX_ASSERT(this->isEmpty());

   if (box.intersects(takeaway)) {
      burstBoxes(box, takeaway, box.getDim().getValue());
   } else {
      appendItem(box);
   }

}

void BoxList::removeIntersections(
   const BoxList& takeaway)
{
   if (isEmpty()) {
      return;
   }

   for (Iterator remove(takeaway); remove; remove++) {
      const Box& byebye = remove();
      removeIntersections(byebye);
   }
}

void BoxList::removeIntersections(
   const BoxTree& takeaway)
{
   if (isEmpty()) {
      return;
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

void BoxList::removeIntersections(
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   const MultiblockBoxTree& takeaway,
   bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
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

void BoxList::intersectBoxes(
   const BoxTree& boxes)
{
   if (isEmpty()) {
      return;
   }

   std::vector<const Box *> overlap_mapped_boxes;
   Box overlap(getFirstItem().getDim());
   Iterator itr(*this);
   Iterator insertion_pt = itr;
   while (itr) {
      const Box& tryme = *itr;
      boxes.findOverlapBoxes(overlap_mapped_boxes, tryme);
      for (size_t i = 0; i < overlap_mapped_boxes.size(); ++i) {
         tryme.intersect(*overlap_mapped_boxes[i], overlap);
         if (!overlap.empty()) {
            addItemAfter(insertion_pt, overlap);
            insertion_pt++;
         }
      }
      overlap_mapped_boxes.clear();
      Iterator tmp = itr;
      insertion_pt++;
      itr = insertion_pt;
      removeItem(tmp);
   }
}

void BoxList::intersectBoxes(
   const BlockId& block_id,
   const IntVector& refinement_ratio,
   const MultiblockBoxTree& boxes,
   bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   const tbox::ConstPointer<hier::GridGeometry>
   & grid_geometry(boxes.getGridGeometry());

   std::vector<const Box *> overlap_mapped_boxes;
   Box overlap(getFirstItem().getDim());
   Iterator itr(*this);
   Iterator insertion_pt = itr;
   while (itr) {
      const Box& tryme = *itr;
      boxes.findOverlapBoxes(overlap_mapped_boxes,
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
               addItemAfter(insertion_pt, overlap);
               insertion_pt++;
            }
         } else {
            tryme.intersect(*overlap_mapped_boxes[i], overlap);
            if (!overlap.empty()) {
               addItemAfter(insertion_pt, overlap);
               insertion_pt++;
            }
         }
      }
      overlap_mapped_boxes.clear();
      Iterator tmp = itr;
      insertion_pt++;
      itr = insertion_pt;
      removeItem(tmp);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Return the boxes in the current list that intersect the index space   *
 * of the argument.                                                      *
 *                                                                       *
 *************************************************************************
 */

void BoxList::intersectBoxes(
   const Box& box)
{
   if (isEmpty()) {
      return;
   }

   Iterator i(*this);
   Box overlap(i().getDim());
   while (i) {
      Box& tryme = *i;
      tryme.intersect(box, overlap);
      if (!overlap.empty()) {
         tryme = overlap;
         i++;
      } else {
         Iterator tmp = i;
         i++;
         removeItem(tmp);
      }
   }
}

void BoxList::intersectBoxes(
   const BoxList& boxes)
{
   if (isEmpty()) {
      return;
   }

   Iterator insertion_pt(*this);
   Box overlap(insertion_pt().getDim());
   while (insertion_pt) {
      Iterator tmp = insertion_pt;
      const Box& tryme = *insertion_pt;
      for (Iterator boxes_itr(boxes); boxes_itr; boxes_itr++) {
         tryme.intersect(boxes_itr(), overlap);
         if (!overlap.empty()) {
            addItemAfter(insertion_pt, overlap);
            insertion_pt++;
         }
      }
      insertion_pt++;
      removeItem(tmp);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Coalesce boxes on the list where possible.  The resulting box list    *
 * will contain a non-overlapping set of boxes covering the identical    *
 * region of index space covered by the original list.  Two boxes may be *
 * coalesced if their union is a box (recall that union is not closed    *
 * over boxes), and they have a non-empty intersection or they are       *
 * adjacent to each other in index space.  Empty boxes on the list are   *
 * removed during this process.  Also, the boxes are coalesced in the    *
 * order in which they appear on the list.  No attempt is made to        *
 * coalesce boxes in any particular way (e.g., to achieve the smallest   *
 * number of boxes).                                                     *
 *                                                                       *
 *************************************************************************
 */

void BoxList::coalesceBoxes()
{
   Iterator tb = this->listStart();
   while (tb) {

      bool found_match = false;

      Iterator tb2 = tb;
      tb2++;

      while (!found_match && tb2) {

         if (tb2().coalesceWith(tb())) {
            found_match = true;
            removeItem(tb);
         }

         tb2++;
      }

      if (found_match) {
         tb = this->listStart();
      } else {
         tb++;
      }
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Sort boxes in list from largest to smallest in size with a heap sort. *
 *                                                                       *
 *************************************************************************
 */

void BoxList::heapify(
   Box** heap,
   const int i,
   const int j)
{
   const int l = 2 * i + 1;
   const int r = l + 1;
   int s = i;
   if ((l < j) && (heap[s]->size() > heap[l]->size())) s = l;
   if ((r < j) && (heap[s]->size() > heap[r]->size())) s = r;
   if (s != i) {
      Box* tmp = heap[s];
      heap[s] = heap[i];
      heap[i] = tmp;
      heapify(heap, s, j);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Count the total size of all of the boxes in the boxlist.              *
 *                                                                       *
 *************************************************************************
 */

int BoxList::getTotalSizeOfBoxes() const
{
   int size = 0;
   for (Iterator i(*this); i; i++) {
      size += i().size();
   }
   return size;
}

/*
 *************************************************************************
 *                                                                       *
 * Perform simple operations (contains, grow, shift, refine, coarsen)    *
 * on all elements in the box list.  These functions simply iterate over *
 * all list boxes and apply the operation to each box.                   *
 *                                                                       *
 *************************************************************************
 */

bool BoxList::contains(
   const Index& p) const
{
   for (Iterator i(*this); i; i++) {
      if (i().contains(p)) return true;
   }
   return false;
}

void BoxList::grow(
   const IntVector& ghosts)
{
   for (Iterator i(*this); i; i++) {
      i().grow(ghosts);
   }
}

void BoxList::shift(
   const IntVector& offset)
{
   for (Iterator i(*this); i; i++) {
      i().shift(offset);
   }
}

void BoxList::rotate(
   const Transformation::RotationIdentifier rotation_ident)
{
   if (getDim().getValue() == 2 || getDim().getValue() == 3) {
      for (Iterator i(*this); i; i++) {
         i().rotate(rotation_ident);
      }
   } else {
      NULL_USE(rotation_ident);

      TBOX_ERROR("BoxList::rotate() error ..."
         << "\n   Rotation only implemented for 2D and 3D " << std::endl);
   }
}

void BoxList::refine(
   const IntVector& ratio)
{
   for (Iterator i(*this); i; i++) {
      i().refine(ratio);
   }
}

void BoxList::coarsen(
   const IntVector& ratio)
{
   for (Iterator i(*this); i; i++) {
      i().coarsen(ratio);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Test the box list for intersections among its boxes.                  *
 *                                                                       *
 *************************************************************************
 */

bool BoxList::boxesIntersect() const
{
   bool intersections = false;

   Iterator tryMe(*this);
   Iterator whatAboutMe(*this);
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
 *                                                                       *
 * Return the bounding box for all boxes in the box list.                *
 *                                                                       *
 *************************************************************************
 */

Box BoxList::getBoundingBox() const
{
   if (this->isEmpty()) {
      TBOX_WARNING("Bounding box list is empty");
      const tbox::Dimension dim(tbox::Dimension::getInvalidDimension());
      Box empty(dim);
      return empty;
   } else {
      const tbox::Dimension& dim(this->getFirstItem().getDim());

      Box bbox(dim);

      for (Iterator i(*this); i; i++) {
         bbox += i();
      }

      return bbox;
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Print boxes in list.                                                  *
 *                                                                       *
 *************************************************************************
 */

void BoxList::print(
   std::ostream& os) const
{
   int i = 0;
   for (Iterator b(*this); b; b++) {
      os << "Box # " << i << ":  " << b() << std::endl;
      i++;
   }
}
#endif
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
