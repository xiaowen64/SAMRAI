/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_BoxContainer
#define included_hier_BoxContainer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"

#include <iostream>
#include <list>
#include <set>


namespace SAMRAI {
namespace hier {

class BoxContainerIterator;
class BoxContainerConstIterator;
class BoxTree;
class MultiblockBoxTree;

/*!
 * @brief A generic container for Boxes.
 *
 * This container makes use of the semantics of a list which implies ordering.
 * The ordering of the Boxes in the container is determined by the user of the
 * container.  Specifically, the ordering is explicitly determined by how the
 * user inserts Boxes into the container.
 *
 * @see hier::Box
 */
class BoxContainer
{
friend class BoxContainerIterator;
friend class BoxContainerConstIterator;

public:
   // Typedefs.

   /*!
    * @brief The iterator for class BoxContainer.
    */
   typedef BoxContainerIterator Iterator;

   /*!
    * @brief The const iterator for class BoxContainer.
    */
   typedef BoxContainerConstIterator ConstIterator;

   // Constructors.

   /*!
    * @brief Default constructor creates empty container in unordered state.
    */
   explicit BoxContainer();

   /*!
    * @brief Creates empty container in state determined by boolean
    *
    * param[in] ordered   Container will be ordered if true, unordered if false.
    */
   explicit BoxContainer(const bool ordered);

   /*!
    * @brief Create container containing members from another container.
    *
    * Members in the range [first, last) are copied to new container.
    *
    * @param[in] first
    * @param[in] last
    */
   explicit BoxContainer(
      Iterator first,
      Iterator last,
      bool ordered = false);

   /*!
    * @brief Create a container with 1 box.
    *
    * @param[in] box Box to copy into new container.
    */
   explicit BoxContainer(
      const Box& box,
      bool ordered = false);

   /*!
    * @brief Copy constructor from another BoxContainer.
    *
    * @param[in] other
    */
   BoxContainer(
      const BoxContainer& other);

   /*!
    * @brief Assignment from other BoxContainer.
    *
    * @param[in] rhs
    */
   BoxContainer&
   operator = (
      const BoxContainer& rhs);

   /*!
    * @brief Assignment from an array of tbox::DatabaseBox objects.
    *
    * @param[in] rhs
    */
   BoxContainer&
   operator = (
      const tbox::Array<tbox::DatabaseBox>& rhs);

   /*!
    * @brief Copy constructor from an array of tbox::DatabaseBox objects.
    *
    * @param[in] other
    */
   explicit BoxContainer(
      const tbox::Array<tbox::DatabaseBox>& other);

   // Destructor.

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainer();

   // Size.

   /*!
    * @brief Return the number of boxes in the container.
    *
    * @return The number of boxes in the container.
    */
   int
   size() const;

   /*!
    * @brief Returns true if size() == 0.
    *
    * @return True if the container is empty.
    */
   bool
   isEmpty() const;

   // Iteration.

   /*!
    * @brief Return a ConstIterator pointing to the start of the container.
    *
    * @return An immutable iterator pointing to the first box.
    */
   ConstIterator
   begin() const;

   /*!
    * @brief Return a ConstIterator pointing to the end of the container.
    *
    * @return An immutable iterator pointing to the last box.
    */
   ConstIterator
   end() const;

   /*!
    * @brief Return an Iterator pointing to the start of the container.
    *
    * @return A mutable iterator pointing to the first box.
    */
   Iterator
   begin();

   /*!
    * @brief Return an Iterator pointing to the end of the container.
    *
    * @return A mutable iterator pointing to the last box.
    */
   Iterator
   end();

   // Access.

   /*!
    * @brief Returns the first element in the container.
    *
    * @return An immutable reference to the first Box in the container.
    */
   const Box&
   front() const;

   /*!
    * @brief Returns the first element in the container.
    *
    * @return An immutable reference to the last Box in the container.
    */
   const Box&
   back() const;

   // Insertion.

   /*!
    * @brief Adds "item" to the "front" of the container.
    *
    * Makes "item" the member of the container pointed to by begin().
    *
    * @param[in] item
    */
   void
   pushFront(
      const Box& item);

   /*!
    * @brief Adds "item" to the "end" of the container.
    *
    * Makes "item" the member of the container pointed to by end().
    *
    * @param[in] item
    */
   void
   pushBack(
      const Box& item);

   /*!
    * @brief Add "item" to specific place in the container.
    *
    * Places "item" immediately before the member of the container pointed
    * to by "iter".
    *
    * @param[in] iter Location to add item before.
    * @param[in] item Box to add to container.
    */
   void
   insertBefore(
      Iterator iter,
      const Box& item);

   /*!
    * @brief Add "item" to specific place in the container.
    *
    * Places "item" immediately after the member of the container pointed
    * to by "iter".
    *
    * @param[in] iter Location to add item after.
    * @param[in] item Box to add to container.
    */
   void
   insertAfter(
      Iterator iter,
      const Box& item);

   /*!
    * @brief Prepends the Boxes in "boxes" to this BoxContainer.
    *
    * "Boxes" will be empty following this operation.
    *
    * @param[in] boxes
    */
   void
   spliceFront(
      BoxContainer& boxes);

   /*!
    * @brief Appends the Boxes in "boxes" to this BoxContainer.
    *
    * "Boxes" will be empty following this operation.
    *
    * @param[in] boxes
    */
   void
   spliceBack(
      BoxContainer& boxes);

   // Erasure.

   /*!
    * @brief Remove the first member of the container.
    */
   void
   popFront();

   /*!
    * @brief Remove the last member of the container.
    */
   void
   popBack();

   /*!
    * @brief Remove the member of the container pointed to by "iter".
    *
    * @param[in] iter
    */
   void
   erase(
      Iterator iter);

   /*!
    * @brief Remove the members of the container in the range [first, last).
    *
    * @param[in] first
    * @param[in] last
    */
   void
   erase(
      Iterator first,
      Iterator last);

   /*!
    * @brief Removes all the members of the container.
    */
   void
   clear();

/*******************************************************/
// set methods

   void order();
   void unorder();

   bool isOrdered() const;

   bool insert(const Box& box);

   Iterator insert ( Iterator position,
                     const Box& box );

   void insert ( ConstIterator first,
                 ConstIterator last );

   void insert ( Iterator first,
                 Iterator last );

   Iterator find(const Box& box) const;
   Iterator lower_bound(const Box& box) const;
   Iterator upper_bound(const Box& box) const;

   int erase(const Box& box);

   void
   swap(BoxContainer& other);

   /*!
    * @brief Insert Box owners into a single set container.
    *
    * @param[out] owners
    */
   void
   getOwners(
      std::set<int>& owners) const;

   /*!
    * @brief Split a BoxContainer into two vector<Box>
    * objects, one containing real Boxes and one containing their
    * periodic images.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] real_mapped_box_vector
    *
    * @param[out] periodic_image_mapped_box_vector
    */
   void
   separatePeriodicImages(
      std::vector<Box>& real_mapped_box_vector,
      std::vector<Box>& periodic_image_mapped_box_vector) const;

   /*!
    * @brief Remove periodic image Boxes.
    */
   void
   removePeriodicImageBoxes();

   /*!
    * @brief Unshift periodic image Boxes
    *
    * Change periodic image Boxes to their unshifted position.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] output_mapped_boxes
    *
    * @param[in] refinement_ratio Refinement ratio where the boxes
    * live.
    */
   void
   unshiftPeriodicImageBoxes(
      BoxContainer& output_mapped_boxes,
      const IntVector& refinement_ratio) const;

   /*!
    * @brief Write the BoxSet to a database.
    */
   void
   putToDatabase(
      tbox::Database& database) const;

   /*!
    * @brief Read the BoxSet from a database.
    */
   void
   getFromDatabase(
      tbox::Database& database);


   /*!
    * @brief Intermediary between BoxContainer and output streams,
    * adding ability to control the output.  See
    * BoxContainer::format().
    */
   class Outputter
   {

      friend std::ostream&
      operator << (
         std::ostream& s,
         const Outputter& f);

private:
      friend class BoxContainer;

      /*!
       * @brief Construct the Outputter with a BoxContainer and the
       * parameters needed to output the BoxContainer to a stream.
       */
      Outputter(
         const BoxContainer& mapped_box_set,
         const std::string& border,
         int detail_depth = 0);

      void
      operator = (
         const Outputter& rhs);               // Unimplemented private.

      const BoxContainer& d_set;

      const std::string d_border;

      const int d_detail_depth;
   };

   /*!
    * @brief Return a object to that can format the BoxContainer for
    * inserting into output streams.
    *
    * Usage example (printing with a tab indentation):
    * @verbatim
    *    cout << "my mapped_boxes:\n" << mapped_boxes.format("\t") << endl;
    * @endverbatim
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 0) const;

   /*!
    * @brief Print the contents of the object recursively.
    *
    * @param[in] output_stream
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
/*
   void
   recursivePrint(
      std::ostream& output_stream,
      const std::string& left_border,
      int detail_depth) const;
*/
 

// end set methods
/**********************************************************/

   // Equivalence.

   /*!
    * @brief Returns true if contents of rhs are identical.
    *
    * @return true if this and rhs are identical.
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const BoxContainer& rhs) const;

   /*!
    * @brief Returns true if contents of rhs not are identical.
    *
    * @return true if this and rhs are not identical.
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const BoxContainer& rhs) const;

   // Container manipulation.

   /*!
    * @brief Place the boxes in the container into a canonical ordering.
    *
    * The canonical ordering for boxes is defined such that boxes that lie
    * next to each other in higher dimensions are coalesced together before
    * boxes that lie next to each other in lower dimensions.  This ordering
    * provides a standard representation that can be used to compare box
    * containers.  The canonical ordering also does not allow any overlap
    * between the boxes in the container.  This routine is potentially
    * expensive, since the running time is \f$O(N^2)\f$ for N boxes.  None
    * of the domain calculus routines call simplify(); all calls to simplify
    * the boxes must be explicit.  Note that this routine is distinct from
    * coalesce(), which is not guaranteed to produce a canonical ordering.
    */
   void
   simplify();

   /*!
    * @brief Combine any boxes in the container which may be coalesced.
    *
    * Two boxes may be coalesced if their union is a box (recall that boxes
    * are not closed under index set unions).  Empty boxes in the container
    * are removed during this process.  Note that this is potentially an
    * expensive calculation (e.g., it will require \f$(N-1)!\f$ box
    * comparisons for a box container with \f$N\f$ boxes in the worst
    * possible case).  So this routine should be used sparingly.  Also note
    * that this routine is different than simplify() since it does not
    * produce a canonical ordering.  In particular, this routine processes
    * the boxes in the order in which they appear in the container, rather
    * than attempting to coalesce boxes along specific coordinate directions
    * before others.
    */
   void
   coalesce();

   // Box calculus.

   /*!
    * @brief Grow boxes in the container by the specified ghost cell width.
    *
    * @param[in] ghosts
    */
   void
   grow(
      const IntVector& ghosts);

   /*!
    * @brief Shift boxes in the container by the specified offset.
    *
    * @param[in] offset
    */
   void
   shift(
      const IntVector& offset);

   /*!
    * @brief Refine boxes in container by the specified refinement ratio.
    *
    * @param[in] ratio
    */
   void
   refine(
      const IntVector& ratio);

   /*!
    * @brief Coarsen boxes in container by the specified coarsening ratio.
    *
    * @param[in] ratio
    */
   void
   coarsen(
      const IntVector& ratio);

   /*!
    * @brief Rotate boxes in container according to rotation_ident.
    *
    * @note Currently works only in 2D.
    *
    * @param[in] rotation_ident
    */
   void
   rotate(
      const Transformation::RotationIdentifier rotation_ident);

   /*!
    * @brief Count total number of indices in the boxes in the container.
    *
    * @return Total number of indices of all boxes in the container.
    */
   int
   getTotalSizeOfBoxes() const;

   /*!
    * @brief Determine if "idx" lies within bounds of boxes in container.
    *
    * @return true if idx lies within bounds of boxes in container.
    *
    * @param[in] idx
    */
   bool
   contains(
      const Index& idx) const;

   /*!
    * @brief Returns the bounding box for all the boxes in the container.
    *
    * @return The bounding box for all the boxes in the container.
    */
   Box
   getBoundingBox() const;

   /*!
    * @brief Check for non-empty intersection among boxes in container.
    *
    * @return Returns true if there exists any non-empty intersection among
    * the boxes in the container.
    */
   bool
   boxesIntersect() const;

   /*!
    * @brief Remove from each box the portions that intersect takeaway.
    *
    * This operation can be thought of as a set difference defined over the
    * abstract AMR box index space.  Performing the set difference will
    * require \f$O(N)\f$ time for a container with \f$N\f$ boxes.  For each
    * box, b, in this container this operation computes b-(b^takeaway) where
    * '^' indicates intersection.
    *
    * @param[in] takeaway What to exclude from each box in the container.
    */
   void
   removeIntersections(
      const Box& takeaway);

   /*!
    * @brief Remove from each box portions intersecting boxes in takeaway.
    *
    * For each box, b, in this container and for each box, t, in takeaway
    * this operation computes b-(b^t) where '^' indicates intersection.
    *
    * @param[in] takeaway What to exclude from each box in the container.
    */
   void
   removeIntersections(
      const BoxContainer& takeaway);

   /*!
    * @brief Remove from each box portions intersecting boxes in takeaway.
    *
    * BoxTree has an efficient overlap search method so this
    * version of removeIntersection is relatively fast.
    * For each box, b, in this container and for each box, t, in takeaway
    * this operation computes b-(b^t) where '^' indicates intersection.
    *
    * @param[in] takeaway What to exclude from each box in the container.
    */
   void
   removeIntersections(
      const BoxTree& takeaway);

   /*!
    * @brief Remove from each box portions intersecting boxes in takeaway.
    *
    * Use extra data to provide needed information in a multiblock setting.
    * MultiblockBoxTree has an efficient overlap search method so this
    * version of removeIntersection is relatively fast.
    *
    * @param[in] block_id Assume all boxes in this BoxContainer belong in
    * the index space specified this BlockId.
    *
    * @param[in] refinement_ratio Assume all boxes in this BoxContainer
    * belong in this refinement ratio.
    *
    * @param[in] takeaway The boxes to take away from this BoxContainer.
    */
   void
   removeIntersections(
      const BlockId& block_id,
      const IntVector& refinement_ratio,
      const MultiblockBoxTree& takeaway,
      bool include_singularity_block_neighbors = false);

   /*!
    * @brief Remove from box the portions intersecting takeaway.
    *
    * This is special version for the case where the container is empty
    * initially.  Upon completion this container contains the result of the
    * removal from box of the intersection of box with takeaway.  If the
    * boxes do not intersect, box is simply added to this container.  This
    * routine is primarily suited for applications which are looking only
    * for the intersection of two boxes.  This operation computes
    * box-(box^takeaway) where '^' indicates intersection.
    *
    * @param[in] box
    * @param[in] takaway
    */
   void
   removeIntersections(
      const Box& box,
      const Box& takeaway);

   /*!
    * @brief Keep the intersection of the container's boxes and keep.
    *
    * Performing the intersection will require \f$O(N)\f$ time for a
    * container with \f$N\f$ boxes.  The complement of removeIntersections.
    *
    * @param[in] keep
    */
   void
   intersectBoxes(
      const Box& keep);

   /*!
    * @brief Keep the intersection of the container's boxes and keep's boxes
    *
    * Intersect the boxes in the current container against the boxes in the
    * specified container.  The intersection calculation will require
    * \f$O(N^2)\f$ time for containers with \f$N\f$ boxes.  The complement
    * of removeIntersections.
    *
    * @param[in] keep
    */
   void
   intersectBoxes(
      const BoxContainer& keep);

   /*!
    * @brief Keep the intersection of the container's boxes and keep's boxes
    *
    * BoxTree has an efficient overlap search method so this
    * version of intersectBoxes is relatively fast.  The complement of
    * removeIntersections.
    *
    * @param[in] keep
    */
   void
   intersectBoxes(
      const BoxTree& keep);

   /*!
    * @brief Keep the intersection of the container's boxes and keep's boxes
    *
    * Use extra data to provide needed information in a multiblock setting.
    * MultiblockBoxTree has an efficient overlap search method so this
    * version of intersectBoxes is relatively fast.  The complement of
    * removeIntersection.
    *
    * @param[in] block_id Assume all boxes in this BoxContainer belong in
    * the index space specified this BlockId.
    *
    * @param[in] refinement_ratio Assume all boxes in this BoxContainer
    * belong in this refefinement ratio.
    *
    * @param[in] boxes The boxes to intersect with this BoxContainer.
    */
   void
   intersectBoxes(
      const BlockId& block_id,
      const IntVector& refinement_ratio,
      const MultiblockBoxTree& keep,
      bool include_singularity_block_neighbors = false);

   /*!
    * @brief Returns a BoxContainer containing the Boxes from this container 
    * in the requested block.
    */
   void
   getSingleBlockBoxContainer(
      BoxContainer& container,
      const BlockId& which_block) const;

   // Database I/O.

   /*!
    * @brief Conversion from BoxContainer to tbox::Array<tbox::DatabaseBox>.
    */
   operator tbox::Array<tbox::DatabaseBox>() const;

   // Debug output.

   /*!
    * @brief Print each box in the container to the specified output stream.
    *
    * @param[in] os
    */
   void
   print(
      std::ostream& os = tbox::plog) const;

private:

   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_BOX_CONTAINER_VERSION;

   /*!
    * @brief Break up bursty against solid and adds the pieces to container.
    *
    * The bursting is done on dimensions 0 through dimension-1, starting
    * with lowest dimensions first to try to maintain the canonical
    * representation for the bursted domains.
    *
    * @param[in] bursty
    * @param[in] solid
    * @param[in] dimension
    */
   void
   burstBoxes(
      const Box& bursty,
      const Box& solid,
      const int dimension);

   /*!
    * @brief Break up bursty against solid and adds the pieces to container
    * starting at location pointed to by itr.
    *
    * The bursting is done on dimensions 0 through dimension-1, starting
    * with lowest dimensions first to try to maintain the canonical
    * representation for the bursted domains.
    *
    * @param[in] bursty
    * @param[in] solid
    * @param[in] dimension
    * @param[in] itr
    */
   void
   burstBoxes(
      const Box& bursty,
      const Box& solid,
      const int dimension,
      Iterator& itr);

   /*!
    * @brief Remove from each box in the sublist of this container defined
    * by sublist_start and sublist_end portions intersecting takeaway.
    *
    * @param[in] takeaway
    * @param[in] sublist_start
    * @param[in] sublist_end
    * @param[in] insertion_pt Where to put new boxes created by this
    * operation.
    */
   void
   removeIntersectionsFromSublist(
      const Box& takeaway,
      Iterator& sublist_start,
      Iterator& sublist_end,
      Iterator& insertion_pt);

   /*
    * The underlying container representation.  This class is a wrapper.
    */
   std::list<Box> d_list;

   std::set<Box*, Box::id_less> d_set;

   bool d_ordered;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainer.I"
#endif

#endif
