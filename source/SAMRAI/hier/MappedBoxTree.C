/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Binary tree of MappedBoxes for overlap searches. 
 *
 ************************************************************************/

#ifndef included_hier_MappedBoxTree_C
#define included_hier_MappedBoxTree_C

#include "SAMRAI/hier/MappedBoxTree.h"
#include "SAMRAI/hier/BoxList.h"
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

tbox::Pointer<tbox::Timer> MappedBoxTree::t_build_tree[tbox::Dimension::
                                                       MAXIMUM_DIMENSION_VALUE];
tbox::Pointer<tbox::Timer> MappedBoxTree::t_search[tbox::Dimension::
                                                   MAXIMUM_DIMENSION_VALUE];
unsigned int MappedBoxTree::s_num_build[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] =
{ 0 };
unsigned int MappedBoxTree::s_num_generate[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int MappedBoxTree::s_num_duplicate[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int MappedBoxTree::s_num_search[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] =
{ 0 };
unsigned int MappedBoxTree::s_num_sorted_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { 0 };
unsigned int MappedBoxTree::s_num_found_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int MappedBoxTree::s_max_sorted_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { 0 };
unsigned int MappedBoxTree::s_max_found_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]
   =
   { 0 };
unsigned int MappedBoxTree::s_max_lin_search[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { 0 };

tbox::StartupShutdownManager::Handler
MappedBoxTree::s_initialize_finalize_handler(
   MappedBoxTree::initializeCallback,
   0,
   0,
   MappedBoxTree::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *************************************************************************
 */

MappedBoxTree::MappedBoxTree():
   d_dim(tbox::Dimension::getInvalidDimension()),
   d_bounding_box(tbox::Dimension::getInvalidDimension()),
   d_block_id(BlockId::invalidId()),
   d_partition_dim(0)
{
   TBOX_ERROR("Using forbidden MappedBoxTree constructor.\n"
              <<"This constructor should never be invoked.");
}

MappedBoxTree::MappedBoxTree(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_bounding_box(dim),
   d_block_id(BlockId::invalidId()),
   d_partition_dim(0)
{
}

MappedBoxTree::MappedBoxTree(
   const tbox::Dimension& dim,
   const MappedBoxSet& mapped_boxes,
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
   for (MappedBoxSet::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end();
        ++ni) {
      TBOX_ASSERT(!ni->getBox().empty());
   }
#endif

   /*
    * Implementation note: We can simply copy mapped_boxes into
    * d_mapped_boxes and call privateGenerateTree using:
    *
    *   d_mapped_boxes.insert(d_mapped_boxes.end(),
    *                         mapped_boxes.begin(),
    *                         mapped_boxes.end());
    *   privateGenerateTree(d_mapped_boxes, d_partition_dim, min_number);
    *
    * However, this extra copy slows things down about 30%.
    * So we live with the repetitious code to do the same thing
    * that privateGenerateTree, except with a MappedBoxSet instead of a
    * std::vector<MappedBox >.
    */

   /*
    * Compute the bounding box for the set of mapped boxes.  Also get
    * BlockId from the given mapped_boxes.
    */
   if ( !mapped_boxes.empty() ) {
      TBOX_ASSERT(mapped_boxes.begin()->getBlockId() != BlockId::invalidId());
      d_block_id = mapped_boxes.begin()->getBlockId();
   }
   for (MappedBoxSet::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end(); ++ni) {
      d_bounding_box += (*ni).getBox();
      TBOX_ASSERT(ni->getBlockId() == d_block_id);
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if (mapped_boxes.size() <= min_number) {
      d_mapped_boxes.insert(d_mapped_boxes.end(), mapped_boxes.begin(), mapped_boxes.end());
   }
   else {

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

      std::vector<MappedBox> left_mapped_boxes, right_mapped_boxes;
      for (MappedBoxSet::const_iterator ni = mapped_boxes.begin();
           ni != mapped_boxes.end(); ++ni) {
         const MappedBox& mapped_box = *ni;
         if (mapped_box.getBox().upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
         } else if (mapped_box.getBox().lower(d_partition_dim) > midpoint) {
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
 * Constructs a MappedBoxTree that represents the physical
 * domain specified by box.
 *************************************************************************
 */

MappedBoxTree::MappedBoxTree(
   const tbox::Dimension& dim,
   const std::vector<MappedBox>& mapped_boxes,
   size_t min_number):
   d_dim(dim),
   d_bounding_box(dim),
   d_block_id(BlockId::invalidId())
{
   t_build_tree[d_dim.getValue() - 1]->start();
   ++s_num_build[d_dim.getValue() - 1];
   s_num_sorted_box[d_dim.getValue() - 1] +=
      static_cast<int>(mapped_boxes.size());
   s_max_sorted_box[d_dim.getValue() - 1] = tbox::MathUtilities<int>::Max(
         s_max_sorted_box[d_dim.getValue() - 1],
         static_cast<int>(mapped_boxes.size()));
   min_number = (min_number < 1) ? 1 : min_number;

#ifdef DEBUG_CHECK_ASSERTIONS
   // Catch empty boxes so sorting logic does not have to.
   for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end();
        ++ni) {
      TBOX_ASSERT(!ni->getBox().empty());
   }
#endif

   /*
    * Compute the bounding box for the vector of mapped boxes.  Also get
    * BlockId from the mapped boxes.
    */
   if ( !mapped_boxes.empty() ) {
      TBOX_ASSERT(mapped_boxes.begin()->getBlockId() != BlockId::invalidId());
      d_block_id = mapped_boxes.begin()->getBlockId();
   }
   for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end(); ++ni) {
      d_bounding_box += (*ni).getBox();
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
   }
   else {

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

      std::vector<MappedBox> left_mapped_boxes, right_mapped_boxes;
      for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
           ni != mapped_boxes.end(); ++ni) {
         const MappedBox& mapped_box = *ni;
         if (mapped_box.getBox().upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
         } else if (mapped_box.getBox().lower(d_partition_dim) > midpoint) {
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
 * Constructs a MappedBoxTree that represents the physical
 * domain specified by box.
 *************************************************************************
 */

MappedBoxTree::MappedBoxTree(
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
   for (hier::BoxList::Iterator ni(boxes); ni; ni++) {
      TBOX_ASSERT(!(*ni).empty());
   }
#endif

   /*
    * Compute this mapped_box's domain, which is the bounding box
    * for the list of boxes.
    */
   for (hier::BoxList::Iterator li(boxes); li; li++) {
      d_bounding_box += *li;
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if ((size_t)boxes.size() <= min_number) {
      d_mapped_boxes.reserve(boxes.size());
      LocalId count(-1);
      for (hier::BoxList::Iterator li(boxes); li; li++) {
         const MappedBox n(*li, ++count, 0, d_block_id);
         d_mapped_boxes.insert(d_mapped_boxes.end(), n);
      }
   }
   else {

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

      std::vector<MappedBox> left_mapped_boxes, right_mapped_boxes;
      LocalId count(-1);
      for (hier::BoxList::Iterator li(boxes); li; li++) {
         const MappedBox mapped_box(*li, ++count, 0);
         if (mapped_box.getBox().upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
         } else if (mapped_box.getBox().lower(d_partition_dim) > midpoint) {
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

MappedBoxTree::~MappedBoxTree()
{
}

/*
 *************************************************************************
 * Assignment operator.
 *
 * We share the children with the reference MappedBoxTree.  This is
 * safe because the trees are never changed once they are set up.  If
 * one MappedBoxTree changes, it simply sets its children pointers to
 * NULL instead of trying to change the children.  Other
 * MappedBoxTrees sharing the children do not see any changes to the
 * children.
 *************************************************************************
 */

MappedBoxTree& MappedBoxTree::operator = (
   const MappedBoxTree& r)
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
 * Generate the tree from a given mutable vector of mapped_boxes.
 * The vector will be changed and its output state is undefined.
 *************************************************************************
 */
void MappedBoxTree::generateTree(
   std::vector<MappedBox>& mapped_boxes,
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
   // Ensure all MappedBoxes are all in the same block.
   for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end();
        ++ni) {
      TBOX_ASSERT(!ni->getBox().empty());
      TBOX_ASSERT(ni->getBlockId() == mapped_boxes.begin()->getBlockId());
   }
#endif

   clear();
   privateGenerateTree(mapped_boxes, min_number);
   t_build_tree[d_dim.getValue() - 1]->stop();
}

/*
 *************************************************************************
 * Generate the tree from a given mutable vector of MappedBoxes.
 * The vector will be changed and its output state is undefined.
 *
 * Methods taking various input containers of MappedBoxes could
 * simply copy the input MappedBoxes into a vector, then call this
 * method.  However, we don't do that for efficiency reasons.  The
 * extra copy turns out to be significant.  Therefore, the
 * constructors have code similar to privateGenerateTree to split
 * the incoming MappedBoxes into three groups.  These groups
 * are turned into child branches by setupChildren.
 *
 * This method is not timed using the Timers.  Only the public
 * itnerfaces are timed.  Isolating the recursive code in
 * privateGenerateTree also helps in timing the methods, because timer
 * starts/stops can be removed from the recursive codes.
 *************************************************************************
 */
void MappedBoxTree::privateGenerateTree(
   std::vector<MappedBox>& mapped_boxes,
   size_t min_number)
{
   ++s_num_generate[d_dim.getValue() - 1];

   if (mapped_boxes.size()) {
      d_block_id = mapped_boxes[0].getBlockId();
   }

   /*
    * Compute this tree's domain, which is the bounding box for the
    * constituent boxes.
    */
   for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
        ni != mapped_boxes.end(); ++ni) {
      d_bounding_box += (*ni).getBox();
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if (mapped_boxes.size() <= min_number) {
      d_mapped_boxes.swap(mapped_boxes);
   }
   else {

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

      std::vector<MappedBox> left_mapped_boxes, right_mapped_boxes;
      for (std::vector<MappedBox>::const_iterator ni = mapped_boxes.begin();
           ni != mapped_boxes.end(); ++ni) {
         const MappedBox& mapped_box = *ni;
         if (mapped_box.getBox().upper(d_partition_dim) <= midpoint) {
            left_mapped_boxes.insert(left_mapped_boxes.end(), mapped_box);
         } else if (mapped_box.getBox().lower(d_partition_dim) > midpoint) {
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

   return;
}


/*
**************************************************************************
* This method finishes the tree generation by setting up the children
* branches.  It expects the MappedBoxes be have been split into
* left_mapped_boxes, right_mapped_boxes, and d_mapped_boxes.  It will
* generate the d_left_child and d_right_child.  If d_mapped_boxes is
* big enough, it will generate d_center_child.
*
**************************************************************************
*/
void MappedBoxTree::setupChildren(
   const size_t min_number,
   std::vector<MappedBox> &left_mapped_boxes,
   std::vector<MappedBox> &right_mapped_boxes )
{
   const size_t total_size =
      left_mapped_boxes.size() + right_mapped_boxes.size() + d_mapped_boxes.size();

   /*
    * If all MappedBoxes are in a single child, the child is just as
    * big as its parent, so there is no point recursing.  Put
    * everything into d_mapped_boxes so the check below will prevent
    * recursion.
    */
   if (left_mapped_boxes.size() == total_size) {
      swap(left_mapped_boxes, d_mapped_boxes);
   } else if (right_mapped_boxes.size() == total_size) {
      swap(right_mapped_boxes, d_mapped_boxes);
   }

   /*
    * If d_mapped_boxes is big enough, generate a center child for it.
    */
   if ( d_mapped_boxes.size() > min_number /* recursion criterion */ &&
        d_mapped_boxes.size() < total_size /* avoid infinite recursion */ ) {
      d_center_child = new MappedBoxTree(d_dim);
      d_center_child->privateGenerateTree(d_mapped_boxes, min_number);
      d_mapped_boxes.clear();   // No longer needed for tree construction or search.
   }

   /*
    * Recurse to build this node's left and right children.
    */
   if (!left_mapped_boxes.empty()) {
      d_left_child = new MappedBoxTree(d_dim);
      d_left_child->privateGenerateTree(left_mapped_boxes, min_number);
   }
   if (!right_mapped_boxes.empty()) {
      d_right_child = new MappedBoxTree(d_dim);
      d_right_child->privateGenerateTree(right_mapped_boxes, min_number);
   }

   return;
}




bool MappedBoxTree::hasOverlap(
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   return privateHasOverlap(box);
}

void MappedBoxTree::findOverlapMappedBoxes(
   std::vector<MappedBox>& overlap_mapped_boxes,
   const Box& box) const
{
   ++s_num_search[d_dim.getValue() - 1];
   int num_found_box = static_cast<int>(overlap_mapped_boxes.size());
   t_search[d_dim.getValue() - 1]->start();
   privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
   t_search[d_dim.getValue() - 1]->stop();
   num_found_box = static_cast<int>(overlap_mapped_boxes.size()) -
                   num_found_box;
   s_max_found_box[d_dim.getValue()
                   - 1] =
      tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
         num_found_box);
   s_num_found_box[d_dim.getValue() - 1] += num_found_box;
}

void MappedBoxTree::findOverlapMappedBoxes(
   hier::BoxList& overlap_mapped_boxes,
   const Box& box) const
{
   ++s_num_search[d_dim.getValue() - 1];
   int num_found_box = static_cast<int>(overlap_mapped_boxes.size());
   t_search[d_dim.getValue() - 1]->start();
   privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
   t_search[d_dim.getValue() - 1]->stop();
   num_found_box = static_cast<int>(overlap_mapped_boxes.size()) -
                   num_found_box;
   s_max_found_box[d_dim.getValue()
                   - 1] =
      tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
         num_found_box);
   s_num_found_box[d_dim.getValue() - 1] += num_found_box;
}

void MappedBoxTree::clear()
{
   d_bounding_box.setEmpty();
   d_left_child.setNull();
   d_right_child.setNull();
   d_mapped_boxes.clear();
   d_center_child.setNull();
}

bool MappedBoxTree::isInitialized() const
{
   return !d_bounding_box.empty();
}

const Box& MappedBoxTree::getBoundingBox() const
{
   return d_bounding_box;
}

const tbox::Dimension& MappedBoxTree::getDim() const
{
   return d_dim;
}

void MappedBoxTree::findOverlapMappedBoxes(
   MappedBoxSet& overlap_mapped_boxes,
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   ++s_num_search[d_dim.getValue() - 1];
   int num_found_box = static_cast<int>(overlap_mapped_boxes.size());
   t_search[d_dim.getValue() - 1]->start();
   privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
   t_search[d_dim.getValue() - 1]->stop();
   num_found_box = static_cast<int>(overlap_mapped_boxes.size()) -
                   num_found_box;
   s_max_found_box[d_dim.getValue()
                   - 1] =
      tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
         num_found_box);
   s_num_found_box[d_dim.getValue() - 1] += num_found_box;
}

bool MappedBoxTree::privateHasOverlap(
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   bool has_overlap = false;
   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         has_overlap = d_center_child->privateHasOverlap(box);
      } else {
         for (std::vector<MappedBox>::const_iterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const MappedBox& mapped_box = *ni;
            if (box.intersects(mapped_box.getBox())) {
               has_overlap = true;
               break;
            }
         }
      }

      if (!has_overlap && d_left_child) {
         has_overlap = d_left_child->privateHasOverlap(box);
      }

      if (!has_overlap && d_right_child) {
         has_overlap = d_right_child->privateHasOverlap(box);
      }
   }
   return has_overlap;
}

void MappedBoxTree::privateFindOverlapMappedBoxes(
   MappedBoxSet& overlap_mapped_boxes,
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      } else {
         for (std::vector<MappedBox>::const_iterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const MappedBox& mapped_box = *ni;
            if (box.intersects(mapped_box.getBox())) {
               overlap_mapped_boxes.insert(mapped_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      }

      if (d_right_child) {
         d_right_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      }
   }
}

void MappedBoxTree::privateFindOverlapMappedBoxes(
   std::vector<MappedBox>& overlap_mapped_boxes,
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      } else {
         for (std::vector<MappedBox>::const_iterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const MappedBox& mapped_box = *ni;
            if (box.intersects(mapped_box.getBox())) {
               overlap_mapped_boxes.insert(
                  overlap_mapped_boxes.end(), mapped_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      }

      if (d_right_child) {
         d_right_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      }
   }
}

void MappedBoxTree::privateFindOverlapMappedBoxes(
   hier::BoxList& overlap_mapped_boxes,
   const Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      } else {
         for (std::vector<MappedBox>::const_iterator ni = d_mapped_boxes.begin();
              ni != d_mapped_boxes.end(); ++ni) {
            const MappedBox& mapped_box = *ni;
            if (box.intersects(mapped_box.getBox())) {
               overlap_mapped_boxes.appendItem(mapped_box.getBox());
            }
         }
      }

      if (d_left_child) {
         d_left_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      }

      if (d_right_child) {
         d_right_child->privateFindOverlapMappedBoxes(overlap_mapped_boxes, box);
      }
   }
}

void MappedBoxTree::getMappedBoxes(
   std::vector<MappedBox>& mapped_boxes) const
{
   if (d_center_child) {
      d_center_child->getMappedBoxes(mapped_boxes);
   } else {
      mapped_boxes.insert(
         mapped_boxes.end(),
         d_mapped_boxes.begin(), d_mapped_boxes.end());
   }

   if (d_left_child) {
      d_left_child->getMappedBoxes(mapped_boxes);
   }

   if (d_right_child) {
      d_right_child->getMappedBoxes(mapped_boxes);
   }
}

tbox::Pointer<MappedBoxTree> MappedBoxTree::createRefinedTree(
   const IntVector& ratio) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);
   TBOX_ASSERT(ratio >= IntVector::getOne(d_dim));

   MappedBoxTree* rval = new MappedBoxTree(d_dim);

   rval->d_partition_dim = d_dim.getValue();

   rval->d_bounding_box = d_bounding_box;
   rval->d_bounding_box.refine(ratio);

   rval->d_mapped_boxes = d_mapped_boxes;
   for (std::vector<MappedBox>::iterator ni = rval->d_mapped_boxes.begin();
        ni != rval->d_mapped_boxes.end(); ++ni) {
      (*ni).getBox().refine(ratio);
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

   return tbox::Pointer<MappedBoxTree>(rval);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappedBoxTree::initializeCallback()
{
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++i) {
      const std::string dim_str( tbox::Utilities::intToString(i+1) );
      t_build_tree[i] = tbox::TimerManager::getManager()->
         getTimer(std::string("hier::MappedBoxTree::build_tree[") + dim_str + "]");
      t_search[i] = tbox::TimerManager::getManager()->
                  getTimer(std::string("hier::MappedBoxTree::search[") + dim_str + "]");
   }
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure  *
 * memory for timers does not leak.                                        *
 ***************************************************************************
 */
void MappedBoxTree::finalizeCallback()
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
void MappedBoxTree::resetStatistics(
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
void MappedBoxTree::printStatistics(
   const tbox::Dimension& dim)
{
   tbox::plog << "MappedBoxTree local stats:"
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
