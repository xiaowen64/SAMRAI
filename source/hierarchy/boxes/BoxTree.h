//
// File:        BoxTree.h
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 601 $
// Modified:    $Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: Utility class to reduce complexity of box calculus operations.
//

#ifndef included_hier_BoxTree
#define included_hier_BoxTree

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_BoxTreeNode
#include "BoxTreeNode.h"
#endif
#ifndef included_hier_ProcessorMapping
#include "ProcessorMapping.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif


namespace SAMRAI {
   namespace hier {


/*!
 * Class BoxTree is a utility class that provides functionality
 * that can be used to reduce the runtime complexity of certain box
 * calculus operations.
 * 
 * A BoxTree object is constructed by passing the ctor a BoxArray;
 * this array is used internally to setup private data structures that
 * will be used in subsequent methods, described below. If the
 * BoxArray contains \b n boxes, then the setup phase (which is
 * invoked by the constructor) has an expected runtime complexity 
 * of O(n log(n)).
 * 
 * Following construction, two types of operations are supported by the 
 * member functions findOverlapIndices() and findLocalOverlapIndices().  
 * Each function takes an integer array and a box as input, and computes 
 * the indices of the subset of the boxes (in the array that was passed 
 * to the constructor) that overlap with the argument box.  On return from 
 * the function, the integer array holds the indices of the boxes in the 
 * array that overlap with the argument box.  The function findOverlapIndices()
 * returns indices of all boxes that overlap the argument box.  The function
 * findLocalOverlapIndices() returns only the indices of the boxes that are
 * assigned to this processor.
 * 
 * Second, the removeIntersections call takes a list of boxes as
 * input, and removes those portions that intersect with the boxes
 * in the BoxArray that was passed to the constructor.
 *
 * This class provides functionality that is similar to that provided
 * by the BoxTop and BoxGraph classes.  However, it employs different
 * algorithms to do so, and is expected to have different runtime
 * characteristics.  Questions such as "which is faster for
 * removeIntersections(): BoxTop, BoxGraph, or BoxTree?" depends on many
 * factors for which analysis is difficult.  The answer is largely a
 * function of a particular problem and domain topology, and is best 
 * answered experimentally.
 *
 * 
 * @see hier::BoxTop
 * @see hier::BoxGraph
 * @see hier::Box
 * @see hier::BoxArray
 */
template<int DIM> class BoxTree : public tbox::DescribedClass
{
public:

   //! @{
   /*!
    * @name BoxTree constructors
    */

   /*!
    * @brief Constructor for BoxTree.
    *
    * The primary difference between the constructors is that the first takes
    * processor mapping information while the others do not.  This is because
    * processor mapping information is required in the method
    * findLocalOverlapIndices(), but is not required for the other 
    * methods.
    *
    * If you do not pass processor mappying information, and subsequently call
    * findLocalOverlapIndices(), an unrecoverable assertion will be thrown.
    *
    * @param boxes input array of boxes.
    * @param shifts the amount by which each box is shifted when there
    *               are periodic boundary conditions.  If there are no
    *               periodic boundary conditions you can pass an array
    *               of length zero; otherwise the array must contain
    *               an entry for each box in \b boxes. 
    * @param box_mapping maps each box in \b boxes to a processor.
    * @param min_length controls the partitioning of \b boxes amongst 
    *                   child nodes in the tree.  
    *                   Setting to a larger value tends to
    *                   decrease the total number of nodes in the tree,
    *                   and hence reduces memory requirements; however, a 
    *                   larger value may also increase the cost of the
    *                   findOverlapBoxes and removeIntersections operations.
    *                   (and hence reduces memory requirements), but increase
    *                   the cost of findOverlapBoxes.  Unless you are
    *                   heavily invested in performance tweaking, please use
    *                   the default value.  Really, I mean it.
    */
   BoxTree(const BoxArray<DIM>& boxes,
                 const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
                 const ProcessorMapping& box_mapping,
                 int min_length = 10);

   BoxTree(const BoxArray<DIM>&  boxes,
                 const tbox::Array< tbox::List< IntVector<DIM> > >& box_shifts,
                 int min_length = 10);

   BoxTree(const BoxArray<DIM>& boxes,
                 int min_length = 10);
   //@}

   
   /*!
    * The dtor does nothing interesting.
    */
   ~BoxTree<DIM>();

   //! @{
   /*!
    * @name BoxTree findIndices
    */

   /*!
    * @brief Compute the indices of all boxes that overlap the specified \b box.
    *
    * The integer array will contain the indices of all boxes (in the BoxArray 
    * that was passed to the constructor) and that overlap the specified \b box.  
    *
    * @param indices           integer array to hold the array indices of 
    *                          the overlapping boxes in the box array.
    * @param box               box whose overlaps are requested.
    */
   void findOverlapIndices(tbox::Array<int>& indices,
                           const Box<DIM>& box) const;

   /*!
    * @brief Compute the indices of the boxes local to this processor that 
    * overlap the specified \b box.
    *
    * The integer array will contain the indices of the boxes (in the BoxArray
    * that was passed to the constructor) and that overlap the specified \b box
    * and which are assigned to this processor.  The processor assignment is
    * determined from the mapping information passed to the constructor.  This
    * routine will result in an unrecoverable exception if a constructor that
    * does not accept the mapping information was used to create this object.
    * See the comments for the constructors. 
    *
    * @param indices           integer array to hold the array indices of
    *                          the overlapping boxes in the box array.
    * @param box               box whose overlaps are requested.
    */
   void findLocalOverlapIndices(tbox::Array<int>& indices,
                                const Box<DIM>& box) const;

   //! @}

   /*!
    * @brief Compute list of boxes that overlap the specified \b box.
    *
    * @param overlap_boxes boxes in BoxArray that was passed to the constructor
    *                      that overlap with argument box.
    * @param box           box whose overlaps are requested.
    */
   void findOverlapBoxes(BoxList<DIM>& overlap_boxes, 
                         const Box<DIM>& box) const;

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
  void removeIntersections(BoxList<DIM>& list) const;

private:
   /*!
    * The copy ctor is not implemented.
    */
   BoxTree(const BoxTree<DIM> &box);

   /*!
    * The assignment operator is not implemented.
    */
   BoxTree<DIM>& operator=(const BoxTree<DIM>& box);

   /*!
    * Private routine to construct overlap indices.
    */
   void privateFindOverlapIndices(tbox::Array<int>& indices,
                                  const Box<DIM>& box,
                                  bool find_local_overlaps) const;

   tbox::Pointer< BoxTreeNode<DIM> > d_tree;

   BoxArray<DIM> d_boxes;

   bool d_have_mapping;
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxTree.C"
#endif

#endif

