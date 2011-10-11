/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Binary tree of Boxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_BoxTree_C
#define included_hier_BoxTree_C

#include "SAMRAI/hier/BoxTree.h"

#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

tbox::Pointer<tbox::Timer> BoxTree::t_build_tree[tbox::Dimension::
                                                 MAXIMUM_DIMENSION_VALUE];
tbox::Pointer<tbox::Timer> BoxTree::t_search[tbox::Dimension::
                                             MAXIMUM_DIMENSION_VALUE];
unsigned int BoxTree::s_num_build[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] =
{ 0 };
unsigned int BoxTree::s_num_generate[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int BoxTree::s_num_duplicate[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int BoxTree::s_num_search[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] =
{ 0 };
unsigned int BoxTree::s_num_sorted_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { 0 };
unsigned int BoxTree::s_num_found_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int BoxTree::s_max_sorted_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { 0 };
unsigned int BoxTree::s_max_found_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int BoxTree::s_max_lin_search[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { 0 };

tbox::StartupShutdownManager::Handler
BoxTree::s_initialize_finalize_handler(
   BoxTree::initializeCallback,
   0,
   0,
   BoxTree::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *************************************************************************
 */

BoxTree::BoxTree(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_bounding_box(dim),
   d_block_id(BlockId::invalidId()),
   d_partition_dim(0)
{
}

BoxTree::BoxTree(
   const tbox::Dimension& dim,
   const BoxSet& mapped_boxes,
   size_t min_number):
   d_dim(dim),
   d_bounding_box(dim),
   d_block_id(BlockId::invalidId())
{
   ++s_num_build[d_dim.getValue() - 1];
   s_num_sorted_box[d_dim.getValue() - 1] +=
      static_cast<int>(mapped_boxes.size());
   s_max_sorted_box[d_dim.getValue() - 1] = tbox::MathUtilities<int>::Max(
         s_max_sorted_box[d_dim.getValue() - 1],
         static_cast<int>(mapped_boxes.size()));
   t_build_tree[d_dim.getValue() - 1]->start();
   min_number = (min_number < 1) ? 1 : min_number;

#ifdef DEBUG_CHECK_ASSERTIONS
   // Catch empty boxes so sorting logic does not have to.
   for (BoxSet::ConstIterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end();
        ++ni) {
      TBOX_ASSERT(!ni->empty());
   }
#endif

   /*
    * Implementation note: We can simply copy mapped_boxes into
    * d_mapped_boxes and call privateGenerateTree using:
    *
    *   d_mapped_boxes.insert(mapped_boxes.begin(),
    *                         mapped_boxes.end());
    *   privateGenerateTree(d_mapped_boxes, min_number);
    *
    * However, this extra copy slows things down about 30%.
    * So we live with the repetitious code to do the same thing
    * that privateGenerateTree, except with a BoxSet instead of a
    * std::vector<Box >.
    */

   /*
    * Compute the bounding box for the set of mapped boxes.  Also get
    * BlockId from the given mapped_boxes.
    */
   if (!mapped_boxes.isEmpty()) {
      TBOX_ASSERT(mapped_boxes.begin()->getBlockId() != BlockId::invalidId());
      d_block_id = mapped_boxes.begin()->getBlockId();
   }
   for (BoxSet::ConstIterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end(); ++ni) {
      d_bounding_box += (*ni);
      TBOX_ASSERT(ni->getBlockId() == d_block_id);
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if (mapped_boxes.size() <= min_number) {
      d_mapped_boxes = mapped_boxes; 
   } else {

      /*
       * Partition the boxes into three sets, using the midpoint of
       * the longest dimension of the bounding box:
       *
       * - those that belong to me (intersects the midpoint plane).  Put
       * these in d_mapped_boxes.
       *
       * - those that belong to my left child (lower than the midpoint
       * plane)
       *
       * - those that belong to my right child (higher than the midpoint
       * plane)
       */

      const IntVector bbsize = d_bounding_box.numberCells();
      d_partition_dim = 0;
      for (int d = 1; d < d_dim.getValue(); ++d) {
         if (bbsize(d_partition_dim) < bbsize(d)) {
            d_partition_dim = d;
         }
      }

      int midpoint =
         (d_bounding_box.lower(d_partition_dim)
          + d_bounding_box.upper(d_partition_dim)) / 2;

      BoxSet left_mapped_boxes, right_mapped_boxes;
      for (BoxSet::ConstIterator ni = mapped_boxes.begin();
           ni != mapped_boxes.end(); ++ni) {
         const Box& mapped_box = *ni;
         if (mapped_box.upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
         } else if (mapped_box.lower(d_partition_dim) > midpoint) {
            right_mapped_boxes.insert(right_mapped_boxes.end(), mapped_box);
         } else {
            d_mapped_boxes.insert(d_mapped_boxes.end(), mapped_box);
         }
      }

      setupChildren(min_number, left_mapped_boxes, right_mapped_boxes);

   }

   if (s_max_lin_search[d_dim.getValue() - 1] < d_mapped_boxes.size()) {
      s_max_lin_search[d_dim.getValue() - 1] =
         static_cast<int>(d_mapped_boxes.size());
   }

   t_build_tree[d_dim.getValue() - 1]->stop();
}

/*
 *************************************************************************
 * Constructs a BoxTree that represents the physical
 * domain specified by boxes.
 *************************************************************************
 */

BoxTree::BoxTree(
   const tbox::Dimension& dim,
   const hier::BoxList& boxes,
   const BlockId& block_id,
   size_t min_number):
   d_dim(dim),
   d_bounding_box(d_dim),
   d_block_id(block_id)
{
   t_build_tree[d_dim.getValue() - 1]->start();
   ++s_num_build[d_dim.getValue() - 1];
   s_num_sorted_box[d_dim.getValue() - 1] += boxes.size();
   s_max_sorted_box[d_dim.getValue() - 1] = tbox::MathUtilities<int>::Max(
         s_max_sorted_box[d_dim.getValue() - 1],
         boxes.size());
   min_number = (min_number < 1) ? 1 : min_number;

#ifdef DEBUG_CHECK_ASSERTIONS
   // Catch empty boxes so sorting logic does not have to.
   // Make sure that all boxes have d_block_id.
   for (BoxList::ConstIterator ni(boxes); ni != boxes.end(); ++ni) {
      TBOX_ASSERT(!(*ni).empty());
   }
#endif

   /*
    * Compute this mapped_box's domain, which is the bounding box
    * for the list of boxes.
    */
   for (BoxList::ConstIterator li(boxes); li != boxes.end(); ++li) {
      d_bounding_box += *li;
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if ((size_t)boxes.size() <= min_number) {
      if (boxes.isOrdered()) {
         d_mapped_boxes = boxes;
      } else {
         d_mapped_boxes.order();
         LocalId count(-1);
         for (BoxList::ConstIterator li(boxes); li != boxes.end(); ++li) {
            const Box n(*li, ++count, 0, d_block_id);
            d_mapped_boxes.insert(d_mapped_boxes.end(), n);
         }
      }
   } else {

      /*
       * Partition the boxes into three sets, using the midpoint of
       * the longest dimension of the bounding box:
       *
       * - those that belong to me (intersects the midpoint plane).  Put
       * these in d_mapped_boxes.
       *
       * - those that belong to my left child (lower than the midpoint
       * plane)
       *
       * - those that belong to my right child (higher than the midpoint
       * plane)
       */

      const IntVector bbsize = d_bounding_box.numberCells();
      d_partition_dim = 0;
      for (int d = 1; d < dim.getValue(); ++d) {
         if (bbsize(d_partition_dim) < bbsize(d)) {
            d_partition_dim = d;
         }
      }

      int midpoint =
         (d_bounding_box.lower(d_partition_dim)
          + d_bounding_box.upper(d_partition_dim)) / 2;

      BoxSet left_mapped_boxes, right_mapped_boxes;
      LocalId count(-1);
      Box mapped_box(d_dim);
      for (BoxList::ConstIterator li(boxes); li != boxes.end(); ++li) {
         if (boxes.isOrdered()) {
            mapped_box = *li;
         } else {
            mapped_box.initialize(*li, ++count, 0, d_block_id);
         }
         if (mapped_box.upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
         } else if (mapped_box.lower(d_partition_dim) > midpoint) {
            right_mapped_boxes.insert(right_mapped_boxes.end(), mapped_box);
         } else {
            d_mapped_boxes.insert(d_mapped_boxes.end(), mapped_box);
         }
      }

      setupChildren(min_number, left_mapped_boxes, right_mapped_boxes);

   }

   if (s_max_lin_search[d_dim.getValue() - 1] < d_mapped_boxes.size()) {
      s_max_lin_search[d_dim.getValue() - 1] =
         static_cast<int>(d_mapped_boxes.size());
   }

   t_build_tree[d_dim.getValue() - 1]->stop();
}

/*
 *************************************************************************
 * dtor
 *************************************************************************
 */

BoxTree::~BoxTree()
{
}

/*
 *************************************************************************
 * Assignment operator.
 *
 * We share the children with the reference BoxTree.  This is
 * safe because the trees are never changed once they are set up.  If
 * one BoxTree changes, it simply sets its children pointers to
 * NULL instead of trying to change the children.  Other
 * BoxTrees sharing the children do not see any changes to the
 * children.
 *************************************************************************
 */

BoxTree& BoxTree::operator = (
   const BoxTree& r)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, r);

   d_bounding_box = r.d_bounding_box;
   d_left_child = r.d_left_child;
   d_right_child = r.d_right_child;
   d_mapped_boxes = r.d_mapped_boxes;
   d_partition_dim = r.d_partition_dim;
   d_center_child = r.d_center_child;
   return *this;
}

/*
 *************************************************************************
 * Generate the tree from a given mutable set of mapped_boxes.
 *************************************************************************
 */
void BoxTree::generateTree(
   BoxSet& mapped_boxes,
   size_t min_number)
{
   t_build_tree[d_dim.getValue() - 1]->start();
   ++s_num_build[d_dim.getValue() - 1];
   s_num_sorted_box[d_dim.getValue() - 1] +=
      static_cast<int>(mapped_boxes.size());
   s_max_sorted_box[d_dim.getValue() - 1] = tbox::MathUtilities<int>::Max(
         s_max_sorted_box[d_dim.getValue() - 1],
         static_cast<int>(mapped_boxes.size()));

#ifdef DEBUG_CHECK_ASSERTIONS
   // Catch empty boxes so sorting logic does not have to.
   // Ensure all Boxes are all in the same block.
   for (BoxSet::ConstIterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end();
        ++ni) {
      TBOX_ASSERT(!ni->empty());
      TBOX_ASSERT(ni->getBlockId() == mapped_boxes.begin()->getBlockId());
   }
#endif

   clear();
   mapped_boxes.swap(d_mapped_boxes);
   privateGenerateTree(min_number);
   t_build_tree[d_dim.getValue() - 1]->stop();
}

/*
 *************************************************************************
 * Generate the tree from a given mutable vector of Boxes.
 * The vector will be changed and its output state is undefined.
 *
 * Methods taking various input containers of Boxes could
 * simply copy the input Boxes into a vector, then call this
 * method.  However, we don't do that for efficiency reasons.  The
 * extra copy turns out to be significant.  Therefore, the
 * constructors have code similar to privateGenerateTree to split
 * the incoming Boxes into three groups.  These groups
 * are turned into child branches by setupChildren.
 *
 * This method is not timed using the Timers.  Only the public
 * itnerfaces are timed.  Isolating the recursive code in
 * privateGenerateTree also helps in timing the methods, because timer
 * starts/stops can be removed from the recursive codes.
 *************************************************************************
 */
void BoxTree::privateGenerateTree(
   size_t min_number)
{
   ++s_num_generate[d_dim.getValue() - 1];

   if (d_mapped_boxes.size()) {
      d_block_id = d_mapped_boxes.begin()->getBlockId();
   }

   /*
    * Compute this tree's domain, which is the bounding box for the
    * constituent boxes.
    */
   for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
        ni != d_mapped_boxes.end(); ++ni) {
      d_bounding_box += *ni;
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if (d_mapped_boxes.size() > min_number) {
      /*
       * Partition the boxes into three sets, using the midpoint of
       * the longest dimension of the bounding box:
       *
       * - those that belong to me (intersects the midpoint plane).  Put
       * these in d_mapped_boxes.
       *
       * - those that belong to my left child (lower than the midpoint
       * plane)
       *
       * - those that belong to my right child (higher than the midpoint
       * plane)
       */

      const IntVector bbsize = d_bounding_box.numberCells();
      d_partition_dim = 0;
      for (int d = 1; d < d_dim.getValue(); ++d) {
         if (bbsize(d_partition_dim) < bbsize(d)) {
            d_partition_dim = d;
         }
      }

      int midpoint =
         (d_bounding_box.lower(d_partition_dim)
          + d_bounding_box.upper(d_partition_dim)) / 2;

      BoxSet left_mapped_boxes, right_mapped_boxes;
      for (BoxSet::Iterator ni = d_mapped_boxes.begin();
           ni != d_mapped_boxes.end(); ) {
         const Box& mapped_box = *ni;
         if (mapped_box.upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
            BoxSet::Iterator curr = ni;
            ++ni;
            d_mapped_boxes.erase(curr);
         } else if (mapped_box.lower(d_partition_dim) > midpoint) {
            right_mapped_boxes.insert(right_mapped_boxes.end(), mapped_box);
            BoxSet::Iterator curr = ni;
            ++ni;
            d_mapped_boxes.erase(curr);
         } else {
            ++ni;
         }
      }

      setupChildren(min_number, left_mapped_boxes, right_mapped_boxes);
   }

   if (s_max_lin_search[d_dim.getValue() - 1] < d_mapped_boxes.size()) {
      s_max_lin_search[d_dim.getValue() - 1] =
         static_cast<int>(d_mapped_boxes.size());
   }
}

/*
 **************************************************************************
 * This method finishes the tree generation by setting up the children
 * branches.  It expects the Boxes be have been split into
 * left_mapped_boxes, right_mapped_boxes, and d_mapped_boxes.  It will
 * generate the d_left_child and d_right_child.  If d_mapped_boxes is
 * big enough, it will generate d_center_child.
 *
 **************************************************************************
 */
void BoxTree::setupChildren(
   const size_t min_number,
   BoxSet& left_mapped_boxes,
   BoxSet& right_mapped_boxes)
{
   const size_t total_size =
      left_mapped_boxes.size() + right_mapped_boxes.size() + d_mapped_boxes.size();

   /*
    * If all Boxes are in a single child, the child is just as
    * big as its parent, so there is no point recursing.  Put
    * everything into d_mapped_boxes so the check below will prevent
    * recursion.
    */
   if (left_mapped_boxes.size() == total_size) {
      left_mapped_boxes.swap(d_mapped_boxes);
   } else if (right_mapped_boxes.size() == total_size) {
      right_mapped_boxes.swap(d_mapped_boxes);
   }

#if 0
   tbox::plog << "Split " << d_mapped_boxes.size() << "  " << d_bounding_box
              << " across " << d_partition_dim << " at " << mid << " into "
              << ' ' << left_mapped_boxes.size()
              << ' ' << cent_mapped_boxes.size()
              << ' ' << right_mapped_boxes.size()
              << std::endl;
#endif
   /*
    * If d_mapped_boxes is big enough, generate a center child for it.
    */
   if (d_mapped_boxes.size() > min_number /* recursion criterion */ &&
       d_mapped_boxes.size() < total_size /* avoid infinite recursion */) {
      d_center_child = new BoxTree(d_dim);
      d_mapped_boxes.swap(d_center_child->d_mapped_boxes);
      d_center_child->privateGenerateTree(min_number);
      d_mapped_boxes.clear();   // No longer needed for tree construction or search.
   }

   /*
    * Recurse to build this node's left and right children.
    */
   if (!left_mapped_boxes.isEmpty()) {
      d_left_child = new BoxTree(d_dim);
      left_mapped_boxes.swap(d_left_child->d_mapped_boxes);
      d_left_child->privateGenerateTree(min_number);
   }
   if (!right_mapped_boxes.isEmpty()) {
      d_right_child = new BoxTree(d_dim);
      right_mapped_boxes.swap(d_right_child->d_mapped_boxes);
      d_right_child->privateGenerateTree(min_number);
   }
}

bool BoxTree::hasOverlap(
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   bool has_overlap = false;
   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         has_overlap = d_center_child->hasOverlap(box);
      } else {
         for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const Box& mapped_box = *ni;
            if (box.intersects(mapped_box)) {
               has_overlap = true;
               break;
            }
         }
      }

      if (!has_overlap && d_left_child) {
         has_overlap = d_left_child->hasOverlap(box);
      }

      if (!has_overlap && d_right_child) {
         has_overlap = d_right_child->hasOverlap(box);
      }
   }
   return has_overlap;
}

void BoxTree::findOverlapBoxes(
   std::vector<Box>& overlap_mapped_boxes,
   const Box& box,
   bool recursive_call) const
{
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      num_found_box = static_cast<int>(overlap_mapped_boxes.size());
      t_search[d_dim.getValue() - 1]->start();
   }

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      } else {
         for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const Box& mapped_box = *ni;
            if (box.intersects(mapped_box)) {
               overlap_mapped_boxes.push_back(mapped_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      }
   }

   if (!recursive_call) {
      t_search[d_dim.getValue() - 1]->stop();
      num_found_box = static_cast<int>(overlap_mapped_boxes.size())
         - num_found_box;
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
}

void BoxTree::findOverlapBoxes(
   std::vector<const Box *>& overlap_mapped_boxes,
   const Box& box,
   bool recursive_call) const
{
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      num_found_box = static_cast<int>(overlap_mapped_boxes.size());
      t_search[d_dim.getValue() - 1]->start();
   }

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      } else {
         for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const Box& mapped_box = *ni;
            if (box.intersects(mapped_box)) {
               overlap_mapped_boxes.push_back(&mapped_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      }
   }

   if (!recursive_call) {
      t_search[d_dim.getValue() - 1]->stop();
      num_found_box = static_cast<int>(overlap_mapped_boxes.size())
         - num_found_box;
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
}

void BoxTree::findOverlapBoxes(
   BoxList& overlap_boxes,
   const Box& box,
   bool recursive_call) const
{
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      num_found_box = static_cast<int>(overlap_boxes.size());
      t_search[d_dim.getValue() - 1]->start();
   }

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_boxes, box, true);
      } else {
         for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const Box& this_box = *ni;
            if (box.intersects(this_box)) {
               overlap_boxes.insert(this_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_boxes, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_boxes, box, true);
      }
   }


   if (!recursive_call) {
      t_search[d_dim.getValue() - 1]->stop();
      num_found_box = static_cast<int>(overlap_boxes.size()) - num_found_box;
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
}

void BoxTree::clear()
{
   d_bounding_box.setEmpty();
   d_left_child.setNull();
   d_right_child.setNull();
   d_mapped_boxes.clear();
   d_center_child.setNull();
}

bool BoxTree::isInitialized() const
{
   return !d_bounding_box.empty();
}

const Box& BoxTree::getBoundingBox() const
{
   return d_bounding_box;
}

const tbox::Dimension& BoxTree::getDim() const
{
   return d_dim;
}
#if 0
void BoxTree::findOverlapBoxes(
   BoxSet& overlap_mapped_boxes,
   const Box& box,
   bool recursive_call) const
{
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      num_found_box = static_cast<int>(overlap_mapped_boxes.size());
      t_search[d_dim.getValue() - 1]->start();
   }

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      } else {
         for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const Box& mapped_box = *ni;
            if (box.intersects(mapped_box)) {
               overlap_mapped_boxes.insert(mapped_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_mapped_boxes, box, true);
      }
   }

   if (!recursive_call) {
      t_search[d_dim.getValue() - 1]->stop();
      num_found_box = static_cast<int>(overlap_mapped_boxes.size())
         - num_found_box;
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
}
#endif
void BoxTree::findOverlapBoxes(
   Connector& overlap_connector,
   const Box& box,
   bool recursive_call) const
{
   const BoxId& box_id = box.getId();
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      if (overlap_connector.hasNeighborSet(box_id)) {
         num_found_box = overlap_connector.numLocalNeighbors(box_id);
      }
      else {
         num_found_box = 0;
      }
      t_search[d_dim.getValue() - 1]->start();
   }

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_connector, box, true);
      } else {
         for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const Box& mapped_box = *ni;
            if (box.intersects(mapped_box)) {
	      overlap_connector.insertLocalNeighbor(mapped_box, box_id);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_connector, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_connector, box, true);
      }
   }

   if (!recursive_call) {
      t_search[d_dim.getValue() - 1]->stop();
      if (overlap_connector.hasNeighborSet(box_id)) {
         num_found_box =
            overlap_connector.numLocalNeighbors(box_id) - num_found_box;
      }
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
   return;
}

void BoxTree::getBoxes(
   std::vector<Box>& mapped_boxes) const
{
   if (d_center_child) {
      d_center_child->getBoxes(mapped_boxes);
   } else {
      for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
           ni != d_mapped_boxes.end(); ++ni) {
         mapped_boxes.push_back(*ni);
      }
   }

   if (d_left_child) {
      d_left_child->getBoxes(mapped_boxes);
   }

   if (d_right_child) {
      d_right_child->getBoxes(mapped_boxes);
   }
}

tbox::Pointer<BoxTree> BoxTree::createRefinedTree(
   const IntVector& ratio) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);
   TBOX_ASSERT(ratio >= IntVector::getOne(d_dim));

   BoxTree* rval = new BoxTree(d_dim);

   rval->d_partition_dim = d_dim.getValue();

   rval->d_bounding_box = d_bounding_box;
   rval->d_bounding_box.refine(ratio);

   for (BoxSet::ConstIterator ni = d_mapped_boxes.begin();
        ni != d_mapped_boxes.end(); ++ni) {
      Box refined_box = *ni;
      refined_box.refine(ratio);
      rval->d_mapped_boxes.insert(refined_box);
   }

   if (!d_center_child.isNull()) {
      rval->d_center_child = d_center_child->createRefinedTree(ratio);
   }
   if (!d_left_child.isNull()) {
      rval->d_left_child = d_left_child->createRefinedTree(ratio);
   }
   if (!d_right_child.isNull()) {
      rval->d_right_child = d_right_child->createRefinedTree(ratio);
   }

   return tbox::Pointer<BoxTree>(rval);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void BoxTree::initializeCallback()
{
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++i) {
      const std::string dim_str(tbox::Utilities::intToString(i + 1));
      t_build_tree[i] = tbox::TimerManager::getManager()->
         getTimer(std::string("hier::BoxTree::build_tree[") + dim_str + "]");
      t_search[i] = tbox::TimerManager::getManager()->
         getTimer(std::string("hier::BoxTree::search[") + dim_str + "]");
   }
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */
void BoxTree::finalizeCallback()
{
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++i) {
      t_build_tree[i].setNull();
      t_search[i].setNull();
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */
void BoxTree::resetStatistics(
   const tbox::Dimension& dim)
{
   s_num_build[dim.getValue() - 1] = 0;
   s_num_generate[dim.getValue() - 1] = 0;
   s_num_duplicate[dim.getValue() - 1] = 0;
   s_num_search[dim.getValue() - 1] = 0;
   s_num_sorted_box[dim.getValue() - 1] = 0;
   s_num_found_box[dim.getValue() - 1] = 0;
   s_max_sorted_box[dim.getValue() - 1] = 0;
   s_max_found_box[dim.getValue() - 1] = 0;
   s_max_lin_search[dim.getValue() - 1] = 0;
}

/*
 ***************************************************************************
 ***************************************************************************
 */
void BoxTree::printStatistics(
   const tbox::Dimension& dim)
{
   tbox::plog << "BoxTree local stats:"
              << "  build=" << s_num_build[dim.getValue() - 1]
              << "  generate=" << s_num_generate[dim.getValue() - 1]
              << "  duplicate=" << s_num_duplicate[dim.getValue() - 1]
              << "  search=" << s_num_search[dim.getValue() - 1]
              << "  sorted_box=" << s_num_sorted_box[dim.getValue() - 1]
              << "  found_box=" << s_num_found_box[dim.getValue() - 1]
              << "  max_sorted_box=" << s_max_sorted_box[dim.getValue() - 1]
              << "  max_found_box=" << s_max_found_box[dim.getValue() - 1]
              << "  max_lin_search=" << s_max_lin_search[dim.getValue() - 1]
              << std::endl;

   tbox::Statistician* st = tbox::Statistician::getStatistician();
   tbox::Pointer<tbox::Statistic> bdstat = st->getStatistic("num_build",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> gnstat = st->getStatistic("num_generate",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> dpstat = st->getStatistic("num_duplicate",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> srstat = st->getStatistic("num_search",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> sbstat = st->getStatistic("num_sorted_box",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> fbstat = st->getStatistic("num_found_box",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> msbstat = st->getStatistic("max_sorted_box",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> mfbstat = st->getStatistic("max_found_box",
         "PROC_STAT");
   tbox::Pointer<tbox::Statistic> lsstat = st->getStatistic("max_lin_search",
         "PROC_STAT");

   static int seq_num = 0;
   bdstat->recordProcStat(s_num_build[dim.getValue() - 1], seq_num);
   gnstat->recordProcStat(s_num_generate[dim.getValue() - 1], seq_num);
   dpstat->recordProcStat(s_num_duplicate[dim.getValue() - 1], seq_num);
   srstat->recordProcStat(s_num_search[dim.getValue() - 1], seq_num);
   sbstat->recordProcStat(s_num_sorted_box[dim.getValue() - 1], seq_num);
   fbstat->recordProcStat(s_num_found_box[dim.getValue() - 1], seq_num);
   lsstat->recordProcStat(s_max_lin_search[dim.getValue() - 1], seq_num);

   st->finalize(false);
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int nproc = mpi.getSize();

   double avg;
   double min, max;
   int rmin(0), rmax(0);

   int doublewidth = 6;
   int intwidth = 6;
   int namewidth = 20;

   min = max = avg = s_num_build[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << bdstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_generate[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << gnstat->getName()
              << "  " << std::setw(namewidth) << std::setw(doublewidth)
              << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_duplicate[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << dpstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_search[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << srstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_sorted_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << sbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_found_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << fbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_max_sorted_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << msbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_max_found_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << mfbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_max_lin_search[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << bdstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   ++seq_num;
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
