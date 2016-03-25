//
// File:	BoxList.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	A list of boxes with basic domain calculus operations
//

#ifndef included_hier_BoxList
#define included_hier_BoxList

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

namespace SAMRAI {
   namespace hier {

template<int DIM> class BoxArray;

/**
 * Class BoxList represents a linked list of boxes.  It defines
 * basic box calculus operations such as set intersection, union, and
 * difference.
 *
 * @see hier::Box
 * @see hier::BoxArray
 */

template<int DIM> class BoxList : public tbox::List< Box<DIM> >
{
public:
   /**
    * The iterator for class BoxList<DIM>.  This is a convenient
    * alias to the list iterator tbox::List< Box<DIM> >::Iterator.
    */
   typedef typename tbox::List< Box<DIM> >::Iterator Iterator;

   /**
    * Create an empty box list with no boxes.
    */
   BoxList();

   /**
    * Create a box list with one box in it.
    */
   BoxList(const Box<DIM>& box);

   /**
    * Create a box list and copy the boxes from the argument list.
    */
   BoxList(const BoxList<DIM>& list);

   /**
    * Create a box list and copy the boxes from the argument array.
    */
   BoxList(const BoxArray<DIM>& array);

   /**
    * Copy boxes from the argument list.
    */
   BoxList<DIM>& operator=(const BoxList<DIM>& list);

   /**
    * The destructor releases all list storage.
    */
   ~BoxList<DIM>();

   /**
    * Return integer number of boxes in the box list.  Note that this 
    * function merely calls the getNumberItems() function in the tbox::List
    * base class.  This function is provided for convenience and consistency
    * with the tbox::BoxArray class.
    */
   int getNumberOfBoxes() const;

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
   void simplifyBoxes();

   /**
    * Add the box to the list of boxes.  Note that this routine does not
    * simplify the box list.  Thus, the new box may overlap with boxes
    * that already reside on the list.
    */
   void unionBoxes(const Box<DIM>& box);

   /**
    * Add the boxes to the list of boxes.  Note that this routine does not
    * simplify the box list.  Thus, the new boxes may overlap with boxes
    * that already reside on the list.
    */
   void unionBoxes(const BoxList<DIM>& boxes);

   /**
    * Remove from the current boxlist the portions that intersect 
    * the box takeaway.  This operation can be thought of as a set 
    * difference defined over the abstract AMR box index space.  
    * Performing the set difference will require \f$O(N)\f$ time for a 
    * list with \f$N\f$ boxes.
    */
   void removeIntersections(const Box<DIM>& takeaway);

   /**
    * Remove from the current boxlist the portions that intersect the
    * boxes in the BoxList takeaway. 
    */
   void removeIntersections(const BoxList<DIM>& takeaway);

   /**
    * A special version for the case where the BoxList is empty initially,
    * this routine builds the list of boxes that get formed when intersecting
    * box with takeaway.  If the boxes do not intersect, box is added to 
    * the boxlist.  This routine is primarily suited for applications
    * which are looking only for the intersection of two boxes.
    */
   void removeIntersections(const Box<DIM>& box, 
                            const Box<DIM>& takeaway);

   /**
    * Intersect the current boxlist against the specified box.  Performing
    * the intersection will require \f$O(N)\f$ time for a list with \f$N\f$ boxes.
    */
   void intersectBoxes(const Box<DIM>& box);

   /**
    * Intersect the current boxlist against the specified boxlist.
    * The intersection calculation will require \f$O(N^2)\f$ time for
    * boxlists with \f$N\f$ boxes.
    */
   void intersectBoxes(const BoxList<DIM>& boxes);

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
   void coalesceBoxes();

   /**
    * Sort the boxes in the list from largest to smallest in size.  Recall 
    * that the size of a box is the number of cell indices it contains.  This
    * routine uses a heap sort algorithm which requires \f$O(N*logN)\f$ work
    * in both average and worst case scenarios.  Note that to obtain a 
    * box list from smallest to largest, the list elements can be reversed
    * after this sorting routine is applied. 
    */
   void sortDescendingBoxSizes();

   /**
    * Count up the total number of indices in all the boxes in the list.
    */
   int getTotalSizeOfBoxes() const;

   /**
    * Check whether an index lies within the bounds of the collection
    * of boxes.
    */
   bool contains(const Index<DIM>& p) const;

   /**
    * Grow all boxes in the box list by the specified ghost cell width.
    */
   void grow(const IntVector<DIM>& ghosts);

   /**
    * Shift all boxes in the box list by the specified offset.
    */
   void shift(const IntVector<DIM>& offset);

   /**
    * Refine the index space of each box in the box list by
    * the specified vector refinement ratio.
    */
   void refine(const IntVector<DIM>& ratio);

   /**
    * Coarsen the index space of each box in the box list by
    * the specified vector coarsening ratio.
    */
   void coarsen(const IntVector<DIM>& ratio);

   /**
    * Return true if there exists non-empty intersection among boxes in
    * list; otherwise, return false.
    */
   bool boxesIntersect() const; 

   /**
    * Return the bounding box for all boxes in the box list.
    */
   Box<DIM> getBoundingBox() const;

   /**
    * Print all class member data for this bounding box list object
    * to specified output stream.
    */
   void print(ostream& os = tbox::plog) const;

private:
   void burstBoxes(const Box<DIM>& bursty,
                   const Box<DIM>& solid,
                   const int dimension);
   void newBurstBoxes(const Box<DIM>& bursty,
                   const Box<DIM>& solid,
                   const int dimension);

   /**
    *
    * Sort boxes in list from largest to smallest in size with a heap sort. 
    *
    */
   static void heapify(Box<DIM>** heap, const int i, const int j);

};


}
}

#ifndef DEBUG_NO_INLINE
#include "BoxList.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxList.C"
#endif

#endif
