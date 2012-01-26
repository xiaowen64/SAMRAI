/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Binary tree of Boxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_BoxTree
#define included_hier_BoxTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/tbox/Timer.h"

#include <boost/shared_ptr.hpp>
#include <vector>
#include <map>

namespace SAMRAI {
namespace hier {

class Connector;

/*!
 * @brief Utility sorting Boxes into tree-like form for finding
 * box overlaps.
 *
 * This class recursively splits a set of Boxes into tree-like
 * form and stores them for fast searches.  The recursive
 * splitting stops when the number of boxes in a leaf node of the tree
 * is less than a minimum number specified in the constructor.
 *
 * All mapped boxes in a BoxTree must exist in the same index space.
 * This means that the must all have the same BlockId value.
 *
 * Overlap searches are done by
 * - hasOverlap()
 * - findOverlapBoxes()
 *
 * Information about the boxes in the tree are given by
 * - getBoundingBox()
 * - getBoxes()
 */

class BoxTree:public tbox::DescribedClass
{
public:

   /*!
    * @brief Constructs a BoxTree from set of Boxes.
    *
    * @param[in] dim
    *
    * @param[in] mapped_boxes.  No empty boxes are allowed.  An assertion
    *                           failure will occur if the mapped boxes in this
    *                           input set do not all have the same BlockId.
    *
    * @param[in] min_number Split up sets of boxes while the number of
    * boxes in a subset is greater than this value.  Setting to a
    * larger value tends to make tree building faster but tree
    * searching slower, and vice versa.  @b Default: 10
    */
   BoxTree(
      const tbox::Dimension& dim,
      const BoxContainer& mapped_boxes,
      int min_number = 10);

   /*!
    * @brief Destructor.
    */
   ~BoxTree();

   /*!
    * @brief Reset to uninitialized state.
    *
    * The dimension of boxes in the tree cannot be changed.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   void
   clear();

   /*!
    * @brief Check whether the tree has been initialized.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   bool
   isInitialized() const;

   //@{

   //! @name Access to state data

   /*!
    * @brief Return the dimension of the boxes in the tree.
    */
   const tbox::Dimension& getDim() const
   {
      return d_dim;
   }

   //@}

   //@{

   //! @name Overlap checks

   /*!
    * @brief Whether the given box has an overlap with Boxes in the
    * tree.
    *
    * @param[in] box The box is assumed to be in same index space as
    * those in the tree.
    */
   bool
   hasOverlap(
      const Box& box) const;

   /*!
    * @brief Find all boxes that overlap the given \b box.
    *
    * Analogous to findOverlapBoxes returning a BoxContainer 
    * but avoids the copies.  If the returned overlap_mapped_boxes are used
    * in a context in which the BoxTree is constant there is no point
    * in incurring the cost of copying the tree's Boxes.  Just return
    * a vector of their addresses.
    *
    * @param[out] overlap_mapped_boxes Pointers to Boxes that overlap
    * with box.
    *
    * @param[in] box the specified box whose overlaps are requested.
    * An assertion failure will occur if the box does not have the same
    * BlockId as the tree.
    */
   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_mapped_boxes,
      const Box& box,
      bool recursive_call = false) const;

   /*!
    * @brief Find all boxes that overlap the given \b box.
    *
    * To avoid unneeded work, the output @b overlap_boxes container
    * is not emptied.  Overlapping Boxes are simply added.
    *
    * @param[out] overlap_boxes Boxes that overlap with box.  The ordered/
    * unordered state of this container is not changed from its state at
    * entry of this method
    *
    * @param[in] box the specified box whose overlaps are requested.
    * An assertion failure will occur if the box does not have the same
    * BlockId as the tree.
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box,
      bool recursive_call = false) const;

   //@}

   /*!
    * @brief Create a similar tree with the boxes refined by a given
    * ratio.
    *
    * @param[in] ratio The boxes are refined by this ratio.
    *
    * Note that there is no coresponding version to create a coarsened
    * tree.  Coarsened trees cannot be trivially generated like
    * refined trees can.  To create a coarsened tree, you must
    * manually get the boxes, coarsen them and use them to build a new
    * tree.
    */
   boost::shared_ptr<BoxTree>
   createRefinedTree(
      const IntVector& ratio) const;

   /*!
    * @brief Assignment operator.
    *
    * @param[in] r
    */
   BoxTree&
   operator = (
      const BoxTree& r);

   /*!
    * @brief Print statistics on number of constructor calls, tree
    * builds, tree searches, etc.
    *
    * This method is for developers to analyze performance.
    */
   static void
   printStatistics(
      const tbox::Dimension& dim);

   /*!
    * @brief Reset statistics on number of constructor calls, tree
    * builds, tree searches, etc.
    *
    * This method is for developers to analyze performance.
    */
   static void
   resetStatistics(
      const tbox::Dimension& dim);

private:
   /*!
    * @brief Default constructor is private to disallow user access.
    * Objects are normally constructed with at least a dimension.
    */
   BoxTree();

   /*!
    * @brief Constructor building an uninitialized object.
    *
    * Private as it is used only internally to create child trees.
    * The object can be initialized using generateTree().
    *
    * @param[in] dim
    */
   explicit BoxTree(
      const tbox::Dimension& dim);

   /*!
    * @brief Private recursive function for generating the search tree.
    *
    * d_mapped_boxes is changed in the process (for efficiency reasons).
    * On output it will contain any boxes that are not assigned to a child
    * tree.
    *
    * The object is not cleared in this method.  If the object has
    * been initialized, it should be cleared before calling this
    * method.  @see clear().
    *
    * @param min_number.  @b Default: 10
    */
   void
   privateGenerateTree(
      int min_number = 10);

   /*!
    * @brief Set up the child branches.
    *
    * This method is called after splitting the Boxes into the
    * left_mapped_boxes and right_mapped_boxes, with boxes straddling
    * the divider stored in d_mapped_boxes.  It generates
    * d_left_child, d_right_child and, if needed, d_center_child.
    *
    * @param[in] min_number
    *
    * @param[in,out] left_mapped_boxes
    *
    * @param[in,out] right_mapped_boxes
    */
   void
   setupChildren(
      const int min_number,
      std::list<const Box*>& left_mapped_boxes,
      std::list<const Box*>& right_mapped_boxes);

   /*!
    * @brief Set up static class members.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief Dimension corresponds to the dimension of boxes in the
    * tree.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Bounding box of all the Boxes in this tree.
    */
   Box d_bounding_box;

   /*!
    * @brief BlockId
    */
   BlockId d_block_id;

   /*!
    * boost::shared_ptrs to familial mapped_boxes.
    */
   boost::shared_ptr<BoxTree> d_left_child;
   boost::shared_ptr<BoxTree> d_right_child;

   /*!
    * @brief A tree for Boxes that are not given to the left or
    * right children.
    */
   boost::shared_ptr<BoxTree> d_center_child;

   /*!
    * @brief Boxes that are contained within the physical domain
    * that this tree represents.  When we have a small number of boxes
    * that do not warant the overhead of a child tree, the boxes go here.
    */
   std::list<const Box*> d_mapped_boxes;

   /*!
    * @brief Dimension along which the input box triples are
    * partitioned.
    */
   int d_partition_dim;

   /*
    * Timers are static to keep the objects light-weight.
    */
   static boost::shared_ptr<tbox::Timer> t_build_tree[tbox::Dimension::
                                                      MAXIMUM_DIMENSION_VALUE];
   static boost::shared_ptr<tbox::Timer> t_search[tbox::Dimension::
                                                  MAXIMUM_DIMENSION_VALUE];

   static unsigned int s_num_build[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_num_generate[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_num_duplicate[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_num_search[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_num_sorted_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_num_found_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_max_sorted_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_max_found_box[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static unsigned int s_max_lin_search[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
