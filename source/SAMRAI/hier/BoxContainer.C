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

BoxContainer::BoxContainer(
   Iterator first,
   Iterator last):
   d_dim(first().getDim()),
   d_list()
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
   if (this != &rhs) {
      d_list = rhs.d_list;
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
   return *this;
}

BoxContainer::BoxContainer(
   const tbox::Array<tbox::DatabaseBox>& other):
   d_dim(other.size() == 0 ? tbox::Dimension::getInvalidDimension() :
                             other[0].getDim())
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
   }
   else {
      insertBefore(tmp, item);
   }
}

/*************************************************************************
 *                                                                       *
 * Function simplify() takes the complicated container of boxes and      *
 * coalesces regions together where possible.                            *
 *                                                                       *
 * The canonical ordering for boxes is defined such that boxes which     *
 * lie next to each other in higher dimensions are coalesced together    *
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
 * them into the container of non-canonical boxes.                       *
 *                                                                       *
 *************************************************************************
 */
void BoxContainer::simplify()
{
   // Start coalescing on the highest dimension of the containers and work down
   // While there are non-canonical boxes, pick somebody out of the container.

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
               for (; l; l++) {
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
                     }
                     else {
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
               }
               else {
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
 *                                                                       *
 * Coalesce boxes in the container where possible.  The resulting box    *
 * container will contain a non-overlapping set of boxes covering the    *
 * identical region of index space covered by the original container.    *
 * Two boxes may be coalesced if their union is a box (recall that union *
 * is not closed over boxes), and they have a non-empty intersection or  *
 * they are adjacent to each other in index space.  Empty boxes in the   *
 * container are removed during this process.  Also, the boxes are       *
 * coalesced in the order in which they appear in the container.  No     *
 * attempt is made to coalesce boxes in any particular way (e.g., to     *
 * achieve the smallest number of boxes).                                *
 *                                                                       *
 *************************************************************************
 */
void BoxContainer::coalesce()
{
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
      }
      else {
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
   for (ConstIterator i(*this); i; i++) {
      size += i().size();
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
 *                                                                       *
 * Return the bounding box for all boxes in the BoxContainer.            *
 *                                                                       *
 *************************************************************************
 */
Box BoxContainer::getBoundingBox() const
{
   if (isEmpty()) {
      TBOX_WARNING("Bounding box container is empty");
      const tbox::Dimension dim(tbox::Dimension::getInvalidDimension());
      Box empty(dim);
      return empty;
   }
   else {
      Box bbox(d_dim);
      for (ConstIterator i(*this); i; i++) {
         bbox += i();
      }

      return bbox;
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Test the box container for intersections among its boxes.             *
 *                                                                       *
 *************************************************************************
 */
bool BoxContainer::boxesIntersect() const
{
   bool intersections = false;

   ConstIterator tryMe(*this);
   ConstIterator whatAboutMe(*this);
   whatAboutMe++;
   while (!intersections && tryMe) {
      while(!intersections && whatAboutMe) {
         if (!((tryMe()*whatAboutMe()).size() == 0)) {
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
 * Return the current container without the portions that intersect      *
 * takeaway.                                                             *
 *                                                                       *
 *************************************************************************
 */
void BoxContainer::removeIntersections(
   const Box& takeaway)
{
   if (isEmpty()) {
      return;
   }

   const unsigned short dim = takeaway.getDim().getValue();
   Iterator insertion_pt(*this);
   while (insertion_pt) {
      Box &tryme = *insertion_pt;
      if (!tryme.intersects(takeaway)) {
         insertion_pt++;
      }
      else {
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

   std::vector<const Box*> overlap_mapped_boxes;
   Iterator itr(*this);
   while (itr) {
      const Box& tryme = *itr;
      takeaway.findOverlapBoxes(overlap_mapped_boxes, tryme);
      if (overlap_mapped_boxes.empty()) {
         itr++;
      }
      else {
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
   const BlockId &block_id,
   const IntVector &refinement_ratio,
   const MultiblockBoxTree& takeaway,
   bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   const tbox::ConstPointer<hier::GridGeometry>
      &grid_geometry(takeaway.getGridGeometry());

   std::vector<const Box*> overlap_mapped_boxes;
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
            const BlockId &overlap_box_block_id =
               overlap_mapped_boxes[i]->getBlockId();
            if ( overlap_box_block_id != block_id ) {
               Box overlap_box = *overlap_mapped_boxes[i];
               grid_geometry->transformBox( overlap_box,
                                            refinement_ratio,
                                            block_id,
                                            overlap_box_block_id );
               removeIntersectionsFromSublist(
                  overlap_box,
                  sublist_start,
                  sublist_end,
                  insertion_pt);
            }
            else {
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
   }
   else {
      pushBack(box);
   }

}

void BoxContainer::removeIntersectionsFromSublist(
   const Box& takeaway,
   Iterator& sublist_start,
   Iterator& sublist_end,
   Iterator& insertion_pt)
{
   const unsigned short dim = takeaway.getDim().getValue();
   Iterator itr = sublist_start;
   while (itr != sublist_end) {
      Box &tryme = *itr;
      if (!tryme.intersects(takeaway)) {
         itr++;
      }
      else {
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
 *                                                                       *
 * Return the boxes in the current container that intersect the index    *
 * space of the argument.                                                *
 *                                                                       *
 *************************************************************************
 */
void BoxContainer::intersectBoxes(
   const Box& keep)
{
   if (isEmpty()) {
      return;
   }

   Iterator i(*this);
   Box overlap(i().getDim());
   while (i) {
      Box &tryMe = *i;
      tryMe.intersect(keep, overlap);
      if (!overlap.empty()) {
         tryMe = overlap;
         i++;
      }
      else {
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

   Iterator insertion_pt(*this);
   Box overlap(insertion_pt().getDim());
   while (insertion_pt) {
      Iterator tmp = insertion_pt;
      const Box &tryme = *insertion_pt;
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

   std::vector<const Box*> overlap_mapped_boxes;
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
   const BlockId &block_id,
   const IntVector &refinement_ratio,
   const MultiblockBoxTree& keep,
   bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   const tbox::ConstPointer<hier::GridGeometry>
      &grid_geometry(keep.getGridGeometry());

   std::vector<const Box*> overlap_mapped_boxes;
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
         const BlockId &overlap_box_block_id =
            overlap_mapped_boxes[i]->getBlockId();
         if ( overlap_box_block_id != block_id ) {
            Box overlap_box = *overlap_mapped_boxes[i];
            grid_geometry->transformBox( overlap_box,
                                         refinement_ratio,
                                         block_id,
                                         overlap_box_block_id );
            tryme.intersect(overlap_box, overlap);
            if (!overlap.empty()) {
               insertAfter(insertion_pt, overlap);
               insertion_pt++;
            }
         }
         else {
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
 *                                                                       *
 * Type conversion from a BoxContainer to an Array of                    *
 * tbox::DatabaseBoxes.                                                  *
 *                                                                       *
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
 *                                                                       *
 * Print boxes in container.                                             *
 *                                                                       *
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
 *                                                                       *
 * Break up box bursty against box solid and adds the pieces to          *
 * container.  The bursting is done on dimensions 0 through dimension-1, *
 * starting with lowest dimensions first to try to maintain the          *
 * canonical representation for the bursted domains.                     *
 *                                                                       *
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
 *                                                                       *
 * Break up box bursty against box solid and adds the pieces to          *
 * container starting at the location pointed to by the supplied         *
 * iterator.  The bursting is done on dimensions 0 through dimension-1,  *
 * starting with lowest dimensions first to try to maintain the          *
 * canonical representation for the bursted domains.                     *
 *                                                                       *
 *************************************************************************
 */
void BoxContainer::burstBoxes(
   const Box& bursty,
   const Box& solid,
   const int dimension,
   Iterator &insertion_pt)
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
