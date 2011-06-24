/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
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

#define MB_MAPPEDBOXTREE_EXISTS

namespace SAMRAI {
namespace hier {

class BoxContainerIterator;
class BoxContainerConstIterator;
class MappedBoxTree;
#ifdef MB_MAPPEDBOXTREE_EXISTS
class MultiblockMappedBoxTree;
#endif

/*!
 * @brief A generic container for Boxes.
 *
 * Boxes (unlike MappedBoxes) do not have any intrinsic notion of ordering.
 * This container makes use of the semantics of a list which implies ordering.
 * The ordering of the Boxes in the container is determined by the user of the
 * container.  Specifically, the ordering is explicitly determined by how the
 * user inserts Boxes into the container.
 *
 * @see hier::Box
 */
class BoxContainer
{
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
       * @brief Create a container with n empty boxes each with dimension dim.
       *
       * @param[in] dim
       * @param[in] n
       */
      explicit BoxContainer(
         const tbox::Dimension& dim,
         const int n = 0);

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
         Iterator last);

      /*!
       * @brief Create a container with 1 box.
       *
       * @param[in] box Box to copy into new container.
       */
      explicit BoxContainer(
         const Box& box);

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
       * @brief Copy constructor from another BoxContainer.
       *
       * @param[in] other
       */
      BoxContainer(
         const BoxContainer& other);

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
       * @return A mutable reference to the first Box in the container.
       */
      Box&
      front();

      /*!
       * @brief Returns the last element in the container.
       *
       * @return A mutable reference to the last Box in the container.
       */
      Box&
      back();

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



      // Equivalence.

      /*!
       * @brief Returns true if contents of rhs are identical.
       *
       * @return true if this and rhs are identical.
       *
       * @param[in] rhs
       */
      bool
      operator == (const BoxContainer& rhs) const;

      /*!
       * @brief Returns true if contents of rhs not are identical.
       *
       * @return true if this and rhs are not identical.
       *
       * @param[in] rhs
       */
      bool
      operator != (const BoxContainer& rhs) const;



      // Box properties.

      /*!
       * @brief Dimension of Box contents.
       *
       * @return The dimension of each Box in the container.
       */
      const tbox::Dimension&
      getDim() const;



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
       * MappedBoxTree has an efficient overlap search method so this
       * version of removeIntersection is relatively fast.
       * For each box, b, in this container and for each box, t, in takeaway
       * this operation computes b-(b^t) where '^' indicates intersection.
       *
       * @param[in] takeaway What to exclude from each box in the container.
       */
      void
      removeIntersections(
         const MappedBoxTree& takeaway);

#ifdef MB_MAPPEDBOXTREE_EXISTS
      /*!
       * @brief Remove from each box portions intersecting boxes in takeaway.
       *
       * Use extra data to provide needed information in a multiblock setting.
       * MultiblockMappedBoxTree has an efficient overlap search method so this
       * version of removeIntersection is relatively fast.
       *
       * @param[in] block_id Assume all boxes in this BoxContainer belong in
       * the index space specified this BlockId.
       *
       * @param[in] refinement_ratio Assume all boxes in this BoxContainer
       * belong in this refefinement ratio.
       *
       * @param[in] takeaway The boxes to take away from this BoxContainer.
       */
      void
      removeIntersections(
         const BlockId &block_id,
         const IntVector &refinement_ratio,
         const MultiblockMappedBoxTree& takeaway);
#endif

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
       * MappedBoxTree has an efficient overlap search method so this
       * version of intersectBoxes is relatively fast.  The complement of
       * removeIntersections.
       *
       * @param[in] keep
       */
      void
      intersectBoxes(
         const MappedBoxTree& keep);

#ifdef MB_MAPPEDBOXTREE_EXISTS
      /*!
       * @brief Keep the intersection of the container's boxes and keep's boxes
       *
       * Use extra data to provide needed information in a multiblock setting.
       * MultiblockMappedBoxTree has an efficient overlap search method so this
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
         const BlockId &block_id,
         const IntVector &refinement_ratio,
         const MultiblockMappedBoxTree& keep);
#endif



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
       * Default constructor just to be clear that there is none.
       */
      BoxContainer();

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

      /*
       * The dimemsion of each box in the container.
       */
      tbox::Dimension d_dim;

      /*
       * The underlying container representation.  This class is a wrapper.
       */
      std::list<Box> d_list;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainer.I"
#endif

#endif
