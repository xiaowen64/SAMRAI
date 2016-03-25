//
// File:        BoxGraph.h
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Utility class to determines topographical box relationships
// 
  
#ifndef included_hier_BoxGraph
#define included_hier_BoxGraph

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_ProcessorMapping
#include "ProcessorMapping.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif


namespace SAMRAI {
   namespace hier {

/*!
 * Class BoxGraph is a utility class that provides functionality
 * that can be used to reduce the runtime complexity of certain box calculus
 * operations.  Two types of functionality are provided.  First,
 * there is a getSrcOverlapIndices() method; second, there is a
 * removeIntersections() method.  These are described below.
 * 
 * BoxGraph is designed around the premise that you have two sets
 * of boxes: an array or source boxes (src_box_array), and an array
 * of destination boxes (dst_box_array).  The dst boxes are grown by
 * some number of cells in each coordinate direction to determine overlaps 
 * with the source boxes.  Typically, this is related to the ghost width of 
 * that data that resides on AMR hierachy patches that correspond to the boxes.
 * The getSrcOverlapIndices() and getLocalSrcOverlapIndices()
 * methods provide efficient methods to determine which boxes from the 
 * src_box_array overlap with a specified box from the dst_box_array.  
 * getSrcOverlapIndices() returns the indices of all source boxes that
 * overlap with a specified destination box. getLocalSrcOverlapIndices()
 * returns the indices of all source boxes that overlap with a specified 
 * destination box and are mapped to the same processor as the specified 
 * destination box.
 * 
 * The removeIntersections() method is similar in concept to the
 * BoxList<DIM>::removeIntersections() method.  Here, however, we
 * operate on the two arrays of boxes that are passed to the constructor.
 * The BoxGraph::removeIntersections() method is functionally equivalent
 * to the following:
 * 
 * @verbatim
 *      BoxList<DIM> dst_box_list(dst_box_array);
 *      BoxList<DIM> src_box_list(src_box_array);
 *      dst_box_list.removeIntersections(src_box_list);
 * @endverbatim
 * 
 * Note that the dst_box_array may be the same as the src_box_array.
 * This occurs, for example, during regridding operations when both
 * both levels are the same.
 * 
 * The above sequence of calls, however, has O(N^2) runtime complexity,
 * assuming both arrays of boxes contain N items.  This class typically
 * performs removeIntersections in O(N^3/2), although worst case behavior
 * can be shown to be O(N^5/2).
 * 
 * Internally BoxGraph operates by constructing a bipartite graph;
 * (this is performed in the constructor).  One set of independent vertices 
 * contains a vertex for each box in the src_box_array; the other independent 
 * set contains vertices corresponding to the dst_box_array.
 * The graph contains an edge (i,j) if
 *    (1) vertices i and j belong to different independent sets;
 *    (2) the boxes corresponding to the vertices overlap each other.
 *
 * Methods in this class were originally designed to improve efficiency
 * in xfer_RefineSchedule<DIM>::generateCommunicationSchedule().
 */

template<int DIM> class BoxGraph : public tbox::DescribedClass
{
public:

   //! @{
   /*!
    * @name BoxGraph constructors
    */

   /*!
    * @brief Constructor for use with getSrcOverlapIndices().
    * 
    * This constructor must be used if you will later call
    * BoxGraph::getSrcOverlapIndices().  This constructor can also be used
    * if you will later call BoxGraph::removeIntersections().
    *
    * @param src_box_array the array of source boxes.
    * @param src_shifts the amount by which each box is shifted when there
    *                   are periodic boundary conditions.  If there are no
    *                   periodic boundary conditions you can pass an array
    *                   of length zero; otherwise the array must contain
    *                   an entry for each box in \b src_box_array.
    * @param src_mapping maps each box in \b src_box_array to a processor.
    * @param dst_box_array the array of destination boxes.
    * @param dst_grow the amount to grow each box in \b dst_box_array.
    * @param sort_dimension the direction (x, y, z) in which boxes will
    *                       be sorted when building the graph.  This parameter
    *                       is included for performance tuning and testing;
    *                       users should use the default setting.
    */
   BoxGraph(const BoxArray<DIM> & src_box_array,
                  const tbox::Array< tbox::List< IntVector<DIM> > >& src_shifts,
                  const ProcessorMapping & src_mapping,
                  const BoxArray<DIM> & dst_box_array,
                  const IntVector<DIM> & dst_grow = IntVector<DIM>(1),
                  int sort_dimension = -1);

   /*!
    * @brief Constructor for use with removeIntersections().
    * 
    * This constructor can be used when a BoxGraph object is employed
    * to perform BoxGraph::removeIntersections() operations.  It can
    * not be used if you intend to call BoxGraph::getSrcOverlapIndices().
    * 
    * @param src_box_array the array of source boxes.
    * @param src_shifts the amount by which each box is shifted when there
    *                   are periodic boundary conditions.  If there are no
    *                   periodic boundary conditions you can pass an array
    *                   of length zero; otherwise the array must contain
    *                   an entry for each box in \b src_box_array.
    * @param dst_box_array the array of destination boxes.
    * @param dst_grow the to grow each box in \b dst_box_array.
    * @param sort_dimension the direction (x, y, z) in which boxes will
    *                       be sorted when building the graph.  This parameter
    *                       is included for performance tuning and testing;
    *                       users should use the default setting.
    */
   BoxGraph(const BoxArray<DIM> & src_box_array,
                  const tbox::Array< tbox::List< IntVector<DIM> > >& src_shifts,
                  const BoxArray<DIM> & dst_box_array,
                  const IntVector<DIM>& dst_grow =  IntVector<DIM>(1),
                  int sort_dimension = -1);

   BoxGraph(const BoxArray<DIM> & src_box_array,
                  const BoxArray<DIM> & dst_box_array,
                  const IntVector<DIM>& dst_grow =  IntVector<DIM>(1),
                  int sort_dimension = -1);
   //@}

   /*!
    * The destructor releases privately held resources.
    */
   ~BoxGraph<DIM>();

   /*!
    * @brief Compute the indices of boxes from the src_box_array
    * that overlap with the box dst_box_array[index].
    * 
    * @param dst_index input; the index of the box from the \b dst_box_array
    *                  whose overlaps are requested.
    * @return the indices of the boxes from \b src_box_array
    *         that overlap with the box: \b dst_box_array[index].
    */
   const tbox::Array<int>& getSrcOverlapIndices(int dst_index);

   /*!
    * @brief Computes the indices of boxes from the src_box_array
    * that overlap with the box dst_box_array[index] and that are
    * mapped to the same processor.
    * 
    * @param dst_index input; the index of the box from the \b dst_box_array
    *                  whose overlaps are requested.
    * @return the indices of the boxes from \b src_box_array
    *         that overlap with the box: \b dst_box_array[index] and 
    *         that are mapped to the same processor.
    */
   const tbox::Array<int>& getLocalSrcOverlapIndices(int dst_index);

   /*!
    * @brief Remove from the dst_box_array the portions that intersect the
    * boxes in src_box_array.
    * 
    * Returns what you would get if you performed the following operations:
    * @verbatim
    *   BoxList<DIM> dst_box_list(dst_box_array);
    *   BoxList<DIM> src_box_list(src_box_array);
    *   dst_box_list.removeIntersections(src_box_list);
    * @endverbatim
    * 
    * @param list on return, contains the portions of the dst_box_array
    *                        that remain after removing the intersections
    *                        with the src_box_array.
    */
   void removeIntersections(BoxList<DIM> &list);

   /*!
    * @brief Undocumented function, used during development and testing.
    */
   void printGraph(ostream& os = tbox::plog);

private:

   BoxGraph(const BoxGraph<DIM>&);  //not implemented
   BoxGraph<DIM>& operator=(const BoxGraph<DIM>&);  //not implemented

   enum Locality { DST_IS_LOCAL, DST_IS_NONLOCAL };
 
   /*
    * The following structure stores information (coordiate, order, box, and
    * whether it is the source box) for each node of the graph.
    */
   struct GraphNode {
      int              coord;    //coordinate of left-hand-side (etc) of box
      int              idx;      //box ordering in original array
      const Box<DIM>  *box;     //the corresponding box
      bool             isSrcBox; //true, if box is on d_src; false is on d_dst
   };


   //private functions
   void addEdges(
      tbox::Array<typename BoxGraph<DIM>::GraphNode> & src, int src_ct,
      tbox::Array<typename BoxGraph<DIM>::GraphNode> & dst, int dst_ct,
      tbox::Array<tbox::List<int> > &adj_lists);

   void buildLocalOverlapArrays();

  //computes the index of the first and last+1 index (wrt nodesIN)
  //of the nodes in the current bin; computes rightEdgeOfBinOUT,
  //which is the cell, wrt the logical space in which the boxes
  //live, of the right edge of the bin.
   void findBinLimits(int binNumberIN, int dim, 
                      tbox::Array<typename BoxGraph<DIM>::GraphNode>& sortedNodesIN,
                      int & firstOUT, int & lastOUT, int & rightEdgeOfBinOUT);

   //returns dimension on which to sort the boxes
   void findSortDimension();

   //used by qsort
   static int qsortCompare(const void *v, const void *w);
   static int qsortIntCompare(const void *v, const void *w);
   void sortNodeList(tbox::Array<typename BoxGraph<DIM>::GraphNode> & list, int len);
   void sortIndices(tbox::Array< tbox::Array<int> > &graph);

   void buildGraph();

   void print(tbox::Array<typename BoxGraph<DIM>::GraphNode> & list, int len);

   //private data members
   BoxArray<DIM>         d_src_boxes;
   tbox::Array<int>        d_src_idx;
   tbox::Array<int>        d_src_map;
   BoxArray<DIM>         d_dst_boxes;

   IntVector<DIM> d_grow_dst_boxes;
   int d_sort_dimension;  //0,  1 or 2

   //Data structures for the computed graph;
   //d_adj[i][0..j] contain the indices, 0..j,
   //of boxes in d_src that overlap the (grown) box
   //d_dst_boxes[i].  The array d_adj_local[i] contains
   //the subset of indices in d_adj[i] that are mapped
   //to the same processor as the box d_dst_boxes[i].
   tbox::Array< tbox::Array<int> >  d_adj;
   tbox::Array< tbox::Array<int> >  d_adj_local;

   /*
    * Timers for performance measurement.
    */
   tbox::Pointer<tbox::Timer>  t_setup_total;
   tbox::Pointer<tbox::Timer>  t_setup_add_edges;
   tbox::Pointer<tbox::Timer>  t_setup_sort;

};


}
}

#ifndef DEBUG_NO_INLINE
#include "BoxGraph.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxGraph.C"
#endif

#endif //included_hier_BoxGraph

