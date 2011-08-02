/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A list of boxes with basic domain calculus operations 
 *
 ************************************************************************/

#ifndef included_hier_BoxList
#define included_hier_BoxList

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MultiblockMappedBoxTree.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/PIO.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

class MappedBoxTree;

/**
 * Class BoxList represents a linked list of boxes.  It defines
 * basic box calculus operations such as set intersection, union, and
 * difference.
 *
 * @see hier::Box
 */

class BoxList:public tbox::List<Box>
{
public:
   /**
    * The iterator for class BoxList.  This is a convenient
    * alias to the list iterator tbox::List< Box >::Iterator.
    */
   typedef tbox::List<Box>::Iterator Iterator;

   /*
    * Default constructor.
    */
   BoxList();

   /**
    * Create an empty box list with no boxes.
    */
   BoxList(const tbox::Dimension& dim);

   /**
    * Create a box list with one box in it.
    */
   explicit BoxList(
      const Box& box);

   /**
    * Create a box list and copy the boxes from the argument list.
    */
   BoxList(
      const BoxList& list);

   /**
    * Create a regular box list from an array of tbox::DatabaseBox objects.
    */
   explicit BoxList(
      const tbox::Array<tbox::DatabaseBox>& array);

   /**
    * Type conversion from BoxList to tbox::Array<tbox::DatabaseBox>.
    */
   operator tbox::Array<tbox::DatabaseBox>() const;

   /**
    * Create a box list using boxes in tbox::Array<tbox::DatabaseBox> for
    * the data.
    */
   BoxList&
   operator = (
      const tbox::Array<tbox::DatabaseBox>& array);

   /**
    * Copy boxes from the argument list.
    */
   BoxList&
   operator = (
      const BoxList& list);

   /**
    * The destructor releases all list storage.
    */
   ~BoxList();

   /**
    * Return integer number of boxes in the box list.  Note that this
    * function merely calls the getNumberOfItems() function in the tbox::List
    * base class.
    */
   int
   getNumberOfBoxes() const;

   /**
    * Return the number of boxes in the list.  Identical to getNumberOfBoxes(),
    * but this method is common to several container classes.
    */
   int
   size() const;

   /**
    * Place the boxes on the list into a canonical ordering.  The canonical
    * ordering for boxes is defined such that boxes that lie next to each
    * other in higher dimensions are coalesced together before boxes that
    * lie next to each other in lower dimensions.  This ordering provides
    * a standard representation that can be used to compare box lists.
    * The canonical ordering also does not allow any overlap between the
    * boxes on the list.  This routine is potentially expensive, since the
    * running time is \f$O(N^2)\f$ for N boxes.  None of the domain calculus
    * routines call simplifyBoxes(); all calls to simplify the boxes must
    * be explicit.  Note that this routine is distinct from coalesceBoxes(),
    * which is not guaranteed to produce a canonical ordering.
    */
   void
   simplifyBoxes();

   /**
    * Add the box to the list of boxes.  Note that this routine does not
    * simplify the box list.  Thus, the new box may overlap with boxes
    * that already reside on the list.
    */
   void
   unionBoxes(
      const Box& box);

   /**
    * Add the boxes to the list of boxes.  Note that this routine does not
    * simplify the box list.  Thus, the new boxes may overlap with boxes
    * that already reside on the list.
    */
   void
   unionBoxes(
      const BoxList& boxes);

   /**
    * Remove from the current boxlist the portions that intersect
    * the box takeaway.  This operation can be thought of as a set
    * difference defined over the abstract AMR box index space.
    * Performing the set difference will require \f$O(N)\f$ time for a
    * list with \f$N\f$ boxes.
    */
   void
   removeIntersections(
      const Box& takeaway);

   /**
    * Remove from the current boxlist the portions that intersect the
    * boxes in the BoxList takeaway.
    */
   void
   removeIntersections(
      const BoxList& takeaway);

   /**
    * Remove from the current boxlist the portions that intersect the
    * boxes in the MappedBoxTree takeaway.
    *
    * MappedBoxTree has an efficient overlap search method so this
    * version of removeIntersection is relatively fast.
    */
   void
   removeIntersections(
      const MappedBoxTree& takeaway);

   /**
    * Remove from the current boxlist the portions that intersect the
    * boxes in the MultiblockMappedBoxTree takeaway.  Use extra data
    * to provide needed information in a multiblock setting.
    *
    * MultiblockMappedBoxTree has an efficient overlap search method
    * so this version of removeIntersection is relatively fast.
    *
    * @param[in] block_id Assume all boxes in this BoxList belong in
    * the index space specified this BlockId.
    *
    * @param[in] refinement_ratio Assume all boxes in this BoxList
    * belong in this refefinement ratio.
    *
    * @param[in] takeaway The boxes to take away from this BoxList.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   void
   removeIntersections(
      const BlockId &block_id,
      const IntVector &refinement_ratio,
      const MultiblockMappedBoxTree& takeaway,
      bool include_singularity_block_neighbors = false);

   /**
    * A special version for the case where the BoxList is empty initially,
    * this routine builds the list of boxes that get formed when intersecting
    * box with takeaway.  If the boxes do not intersect, box is added to
    * the boxlist.  This routine is primarily suited for applications
    * which are looking only for the intersection of two boxes.
    */
   void
   removeIntersections(
      const Box& box,
      const Box& takeaway);

   /**
    * Intersect the current boxlist against the specified box.  Performing
    * the intersection will require \f$O(N)\f$ time for a list with \f$N\f$ boxes.
    */
   void
   intersectBoxes(
      const Box& box);

   /**
    * Intersect the current boxlist against the specified boxlist.
    * The intersection calculation will require \f$O(N^2)\f$ time for
    * boxlists with \f$N\f$ boxes.
    */
   void
   intersectBoxes(
      const BoxList& boxes);

   /**
    * Intersect the current boxlist against the boxes in the specified MappedBoxtree.
    *
    * MappedBoxTree has an efficient overlap search method so this
    * version of intersectBoxes is relatively fast.
    */
   void
   intersectBoxes(
      const MappedBoxTree& boxes);

   /**
    * Intersect the current boxlist against the boxes in the specified
    * MultiblockMappedBoxtree.  Use extra data to provide needed
    * information in a multiblock setting.
    *
    * MultiblockMappedBoxTree has an efficient overlap search method
    * so this version of intersectBoxes is relatively fast.
    *
    * @param[in] block_id Assume all boxes in this BoxList belong in
    * the index space specified this BlockId.
    *
    * @param[in] refinement_ratio Assume all boxes in this BoxList
    * belong in this refefinement ratio.
    *
    * @param[in] boxes The boxes to take intersect with this BoxList.
    *
    * @param[in] include_singularity_block_neighbors Whether to include
    * intersections with boxes in blocks that are neighbors of block
    * block_id across a multiblock singularity.
    */
   void
   intersectBoxes(
      const BlockId &block_id,
      const IntVector &refinement_ratio,
      const MultiblockMappedBoxTree& boxes,
      bool include_singularity_block_neighbors = false);

   /**
    * Combine any boxes in the list which may be coalesced.  Two boxes
    * may be coalesced if their union is a box (recall that boxes are not
    * closed under index set unions).  Empty boxes on the list are removed
    * during this process.  Note that this is potentially an expensive
    * calculation (e.g., it will require \f$(N-1)!\f$ box comparisons for a box
    * list with \f$N\f$ boxes in the worst possible case).  So this routine
    * should be used sparingly.  Also note that this routine is different
    * than simplifyBoxes() since it does not produce a canonical ordering.
    * In particular, this routine processes the boxes in the order in which
    * they appear on the list, rather than attempting to coalesce boxes
    * along specific coordinate directions before others.
    */
   void
   coalesceBoxes();

   /**
    * Count up the total number of indices in all the boxes in the list.
    */
   int
   getTotalSizeOfBoxes() const;

   /**
    * Check whether an index lies within the bounds of the collection
    * of boxes.
    */
   bool
   contains(
      const Index& p) const;

   /**
    * Grow all boxes in the box list by the specified ghost cell width.
    */
   void
   grow(
      const IntVector& ghosts);

   /**
    * Shift all boxes in the box list by the specified offset.
    */
   void
   shift(
      const IntVector& offset);

   /**
    * Rotate the individual boxes in the box array according to rotation_ident.
    * Currently works only in 2D.
    */
   void
   rotate(
      const Transformation::RotationIdentifier rotation_ident);

   /**
    * Refine the index space of each box in the box list by
    * the specified vector refinement ratio.
    */
   void
   refine(
      const IntVector& ratio);

   /**
    * Coarsen the index space of each box in the box list by
    * the specified vector coarsening ratio.
    */
   void
   coarsen(
      const IntVector& ratio);

   /**
    * Return true if there exists non-empty intersection among boxes in
    * list; otherwise, return false.
    */
   bool
   boxesIntersect() const;

   /**
    * Return the bounding box for all boxes in the box list.
    */
   Box
   getBoundingBox() const;

   /**
    * Print all class member data for this bounding box list object
    * to specified output stream.
    */
   void
   print(
      std::ostream& os = tbox::plog) const;

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

private:
   void
   burstBoxes(
      const Box& bursty,
      const Box& solid,
      const int dimension);

   void
   burstBoxes(
      const Box& bursty,
      const Box& solid,
      const int dimension,
      Iterator &itr);

   void
   removeIntersectionsFromSublist(
      const Box& takeaway,
      Iterator& sublist_start,
      Iterator& sublist_end,
      Iterator& insertion_pt);

   /**
    *
    * Sort boxes in list from largest to smallest in size with a heap sort.
    *
    */
   static void
   heapify(
      Box** heap,
      const int i,
      const int j);

   tbox::Dimension d_dim;

};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/BoxList.I"
#endif

#endif
