//
// File:        BoxTop.h
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 601 $
// Modified:    $Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: Utility class to reduce complexity of box calculus operations.
//

#ifndef included_hier_BoxTop
#define included_hier_BoxTop

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
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


/*!
 * Class BoxTop is a utility class that provides functionality
 * that can be used to reduce the runtime complexity of certain box
 * calculus operations.
 * 
 * A BoxTop object is constructed by passing the ctor a BoxArray;
 * this array is used internally to setup private data structures that
 * will be used in subsequent methods, described below. If the
 * BoxArray contains \b n boxes, then the setup phase (which is
 * invoked by the constructor) has a runtime complexity of O(n log(n)).
 * 
 * Following construction, two types of operations are supported.
 * First, the three related functions findOverlappingBoxes,
 * findOverlappingBoxIndices, and findOverlappingBoxesAndIndices
 * take a box as input, and compute the subset of the boxes, and/or
 * the indices of the boxes, that overlap with the box.
 * 
 * Second, the removeIntersections call takes a list of boxes as
 * input, and removes those portions that intersect with the boxes
 * in the BoxArray that was passed to the constructor.
 * 
 * @see hier::Box
 * @see hier::BoxArray
 */

template<int DIM> class BoxTop
{
public:

   //! @{
   /*!
    * @name BoxTop constructors
    */
   /*!
    * @brief Constructor for BoxTop.
    *
    * @param in_boxes input array of boxes.
    * @param in_shifts contains a shift vector associated with each boxes.
    */
   BoxTop(const BoxArray<DIM>& in_boxes,
                const tbox::Array< tbox::List< IntVector<DIM> > >& in_shifts);
   /*!
    * @brief Constructor for BoxTop.
    * 
    * @param in_boxes input array of boxes.
    */
   BoxTop(const BoxArray<DIM>& in_boxes);
   //@}

   /*!
    * The destructor releases privately held resources.
    */
   ~BoxTop<DIM>();

   /*!
    * @brief Compute the set of boxes that overlap with the specified \b box.
    * 
    * The return array, \b overlaps, contains the subset of the boxes
    * that were passed to the constructor that overlap the \b box.
    * 
    * @param overlaps the subset of boxes, from the array that was passed
    *                to the ctor, that overlap the specified box.
    * @param box the specified box whose overlaps are requested.
    */
   void findOverlappingBoxes(BoxArray<DIM> &overlaps,
                             const Box<DIM> & box);

   /*!
    * @brief Compute the indices of the boxes that overlap the specified \b box.
    * 
    * The return array, \b indices, contains the indices of the boxes
    * (in the BoxArray that was passed to the constructor) that overlap
    * the specified \b box.
    * 
    * @param indices the indices of the overlapping boxes.
    * @param box the specified box whose overlaps are requested.
    */
   void findOverlappingBoxIndices(tbox::Array<int> &indices,
                                  const Box<DIM> & box);

   /*!
    * @brief Compute the set of boxes, and their indices, that overlap with
    * the specified \b box.
    * 
    * The return array, \b overlaps, contains the subset of the boxes
    * that were passed to the constructor that overlap the \b box.
    * The returned array, \b indices, contains an index for each box
    * in \b overlaps; each index is the corresponding box's index
    * w.r.t. the BoxArray that was passed to the constructor.
    * 
    * @param overlaps the subset of boxes, from the array that was passed
    *                to the ctor, that overlap the specified box.
    * @param indices the indices of the overlapping boxes.
    * @param box the specified box whose overlaps are requested.
    */
   void findOverlappingBoxesAndIndices(BoxArray<DIM> &overlaps,
                                       tbox::Array<int> &indices,
                                       const Box<DIM> & box);

  /*!
   * @brief Remove from \b list the portions that intersect the boxes
   * in the BoxArray that was passed to the constructor.
   * 
   * CAUTION: the semantics of this call differ from that of
   * BoxList::removeIntersections(const BoxList<DIM> takeaway).
   * Here, the list that is being modified is the list that is
   * passed as an argument; the "takeaway" list is the list that
   * was passed when the BoxTop object was constructed.
   * 
   * @param list the list of boxes from which intersections are to be removed.
   */
  void removeIntersections(BoxList<DIM>& list);

  /*!
   * @brief Undocumented function, used during development and testing.
   */
  void print(ostream& os = tbox::plog);

private:

   //These refer to the manner in which an array of boxElts is sorted.
   enum{ ASCENDING, DESCENDING };

   //Private data structure; arrays of boxElts are used for
   //sorting the boxes by lhs coordinate, etc.
   struct boxElt {
      int        coord;        //coordinate of left-hand-side (etc) of box
      const Box<DIM>  *box;   //the corresponding box
      int        idx;          //box ordering in original list or array
   };

   /*!
    * @brief the copy constructor is not implemented.
    */
   BoxTop(const BoxTop<DIM>& bs);

   /*!
    * @brief the assignment operator is not implemented.
    */
   BoxTop<DIM>& operator=(const BoxTop<DIM>& boxtop);

   /*!
    * @brief This undefined constructor prevents the creation
    * of a "hidden" object.
    */
   BoxTop(const BoxList<DIM>& list);

   //called by constructors; sets up sorted arrays (private data)
   void setup();

   //called by findNabors() public functions.
   void findNaborsPrivate(const Box<DIM> & box);

   //called by setup()
   void sort(boxElt *&array, int len, int ordering);

   //called by findNabors()
   int binSearch(boxElt *tmp, int tmp_len, int coord, int ordering);

   //fills in entries in overlaps and indices by comparing the boxes
   //in the elements of d_shortest_list against the box.
   void buildShortestList(
      BoxArray<DIM> &overlaps,
      bool build_overlaps,
      tbox::Array<int> &indices,
      bool build_indices,
      const Box<DIM>& box);

   //needed for qsort().  Mneumonic: ``A''scending, ``D''escending.
   static int boxEltCompareA(const void *v, const void *w); 
   static int boxEltCompareD(const void *v, const void *w); 

   //called by print()
   void printEltArray(boxElt *&data, int len, ostream& os);

   //The array of boxes that will be tested against during findNabors();
   //same as the array of boxes that was passed to the constructor.
   BoxArray<DIM> d_boxes;

   //Arrays of boxes, sorted by upper or lower coordinate.
   boxElt  *d_sorted_lists[2*DIM];

   //The subset of d_long_list that is found by findNabors().
   //(todo: replace with tbox::Array< tbox::Array<boxElt> > list;)
   boxElt  *d_shortest_list;       
   int     d_shortest_list_length;
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxTop.C"
#endif

#endif 
