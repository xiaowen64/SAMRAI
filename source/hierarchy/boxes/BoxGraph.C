//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/boxes/BoxGraph.C $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1846 $
// Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description: Utility class to determines topographical box relationships
//

#ifndef included_hier_BoxGraph_C
#define included_hier_BoxGraph_C

#include "BoxGraph.h"
#include "BoxGraphUtilities.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#ifdef DEBUG_NO_INLINE
#include "BoxGraph.I"
#endif

namespace SAMRAI {
   namespace hier {

#define  DEFAULT_ADJ_GUESS (DIM*4)

/*
 * ************************************************************************
 * 
 * Constructors and Destructor
 * 
 * ************************************************************************
 */

template<int DIM>  BoxGraph<DIM>::BoxGraph(
   const BoxArray<DIM> & src_box_array,
   const tbox::Array< tbox::List< IntVector<DIM> > >& src_shifts,
   const ProcessorMapping & src_mapping,
   const BoxArray<DIM> & dst_box_array,
   const IntVector<DIM> & dst_grow,
   int sort_dimension)
: 
d_dst_boxes(dst_box_array),
d_grow_dst_boxes(dst_grow),
d_sort_dimension(sort_dimension)
{
   t_setup_total =tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_total");
   t_setup_add_edges =tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_add_edges");
   t_setup_sort = tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_sort");

   d_src_map = src_mapping.getProcessorMapping();
   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
     d_src_boxes, d_src_idx, src_box_array, src_shifts);
   buildGraph();
}


template<int DIM>  BoxGraph<DIM>::BoxGraph(
   const BoxArray<DIM> & src_box_array,
   const tbox::Array< tbox::List< IntVector<DIM> > >& src_shifts,
   const BoxArray<DIM> & dst_box_array,
   const IntVector<DIM> & dst_grow,
   int sort_dimension)
:
d_dst_boxes(dst_box_array), 
d_grow_dst_boxes(dst_grow), 
d_sort_dimension(sort_dimension)
{
   t_setup_total =tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_total");
   t_setup_add_edges =tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_add_edges");
   t_setup_sort = tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_sort");

   BoxGraphUtilities<DIM>::makeBoxesPlusPeriodicBoxes(
      d_src_boxes, d_src_idx, src_box_array, src_shifts);
   buildGraph();
}

template<int DIM>  BoxGraph<DIM>::BoxGraph(
   const BoxArray<DIM> & src_box_array,
   const BoxArray<DIM> & dst_box_array,
   const IntVector<DIM> & dst_grow,
   int sort_dimension)
:
d_src_boxes(src_box_array),
d_dst_boxes(dst_box_array), 
d_grow_dst_boxes(dst_grow), 
d_sort_dimension(sort_dimension)
{
   t_setup_total =tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_total");
   t_setup_add_edges =tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_add_edges");
   t_setup_sort = tbox::TimerManager::getManager()->
      getTimer("hier::BoxGraph::setup_sort");

   buildGraph();
}


template<int DIM>  BoxGraph<DIM>::~BoxGraph()
{
}


/*
 * ******************************************************************
 * Constructs the graph that is used to determing
 * which boxes in d_dst are neighbors of which boxes
 * in d_src.  Two boxes are neighbors if their intersection
 * is non-null.
 * ******************************************************************
 */

template<int DIM> void BoxGraph<DIM>::buildGraph()
{

   t_setup_total->start();

   int src_len = d_src_boxes.getNumberOfBoxes();
   int dst_len = d_dst_boxes.getNumberOfBoxes();
   int total_len = src_len + dst_len;

   /*
    * grow all dst boxes by the specified value
    */
   for (int j=0; j<dst_len; ++j) {
      d_dst_boxes(j).grow(d_grow_dst_boxes);
   }

   /*
    * Build sorted list of nodes.
    */
   tbox::Array<typename BoxGraph<DIM>::GraphNode> sorted_nodes;
   t_setup_sort->start();
   findSortDimension();

   /*
    * Form an array of GraphNodes; there is one GraphNode for
    * each src box and one for every dst box.
    */
   sorted_nodes.resizeArray( total_len );
   int i;
   int idx = 0;
   for (i=0; i<dst_len; ++i) {
      sorted_nodes[idx].coord = d_dst_boxes(i).lower(d_sort_dimension);
      sorted_nodes[idx].idx = i;
      sorted_nodes[idx].box = &(d_dst_boxes(i));
      sorted_nodes[idx].isSrcBox = false;
      ++idx;
   }

   /*
    * If the BoxGraph is being instantiated only to perform
    * removeIntersections, then we don't have any indices;
    * to account for this case, the
    * next two loops are identical except for the statements
    * sorted_nodes[idx].idx = i;              (if)
    * sorted_nodes[idx].idx = d_src_idx[i];   (else)
    */
   if (d_src_idx.getSize() == 0) {
      for (i=0; i<src_len; ++i) {
         sorted_nodes[idx].coord = d_src_boxes(i).lower(d_sort_dimension);
         sorted_nodes[idx].idx = i;
         sorted_nodes[idx].box = &(d_src_boxes(i));
         sorted_nodes[idx].isSrcBox = true;
         ++idx;
      }
   } else {
      for (i=0; i<src_len; ++i) {
         sorted_nodes[idx].coord = d_src_boxes(i).lower(d_sort_dimension);
         sorted_nodes[idx].idx = d_src_idx[i];
         sorted_nodes[idx].box = &(d_src_boxes(i));
         sorted_nodes[idx].isSrcBox = true;
         ++idx;
      }
   }

   sortNodeList(sorted_nodes, total_len);
   t_setup_sort->stop();

   tbox::Array<typename BoxGraph<DIM>::GraphNode> overlap_list(total_len);
   tbox::Array<typename BoxGraph<DIM>::GraphNode> tmp_src_nodes(src_len);
   tbox::Array<typename BoxGraph<DIM>::GraphNode> tmp_dst_nodes(dst_len);

   //initialize storage for the graph's adjacency lists.
   tbox::Array<tbox::List<int> > tmp_graph(dst_len);

   //"processed_ct" is the number of nodes from the sorted
   //list that have been processed (i.e., removed from the list
   //and assigned to the dst or src list for the bin being processed).
   //During each iteration through this loop, we form the subgraph
   //imposed by the boxes that exist in the current bin.
   int first_node, last_node, right_edge_of_bin;
   int bin = -1;
   int overlap_ct = 0;
   int processed_ct = 0;
   while (processed_ct < total_len) {
      ++bin;  //the bin being processed.  

      findBinLimits(bin, d_sort_dimension, sorted_nodes,
                    first_node, last_node, right_edge_of_bin);

      int tmp_src_ct = 0;
      int tmp_dst_ct = 0;

      //Put any boxes that are in both this and the previous
      //bin in the appropriate list.  The boxes that straddle
      //two bins are stored in the "overlap_list."
      for (i=0; i<overlap_ct; ++i) {
         if (overlap_list[i].isSrcBox) {
            tmp_src_nodes[tmp_src_ct++] = overlap_list[i];
         } else {
            tmp_dst_nodes[tmp_dst_ct++] = overlap_list[i];
         }
      }

      //"Pack" the overlap list to remove boxes that do
      // not overlap with the next bin.
      int new_overlap_ct = 0;
      for (i=0; i<overlap_ct; ++i) {
         if ((overlap_list[i].box)->upper(d_sort_dimension) >=
            right_edge_of_bin) {
           overlap_list[new_overlap_ct++] = overlap_list[i];
         }
      }
      overlap_ct = new_overlap_ct;

      /*
       * Pick off the nodes whose left-hand-side is in this bin.
       * Each node is placed in either the dst or src list.
       * In graph terminology, these lists represent a pair of
       * independent sets.
       * Also, a node is appended to the overlap list if its
       * right-hand-side is in the next bin.
       */
      for (i=first_node; i<last_node; ++i) {
         if (sorted_nodes[i].isSrcBox) {
            tmp_src_nodes[tmp_src_ct++] = sorted_nodes[i];
         } else {
            tmp_dst_nodes[tmp_dst_ct++] = sorted_nodes[i];
         }

         if ((sorted_nodes[i].box)->upper(d_sort_dimension) >=
                                          right_edge_of_bin) {
            overlap_list[overlap_ct++] = sorted_nodes[i];
         }

         ++processed_ct;  //increments the outer loop (yuck!)
      }

      /*
       * Add edges that can be discovered in this bin to the graph.
       */
      t_setup_add_edges->start();
      addEdges(tmp_src_nodes, tmp_src_ct, 
               tmp_dst_nodes, tmp_dst_ct,
               tmp_graph);

      t_setup_add_edges->stop();

   } //while (processed_ct < total_len)

   /*
    * shrink all dst boxes by the specified value;
    * this must be done or removeIntersections will not work
    * properly if a "tryme" box has no nabors!
    */
   d_grow_dst_boxes *= -1;
   for (int j=0; j<dst_len; ++j) {
      d_dst_boxes(j).grow(d_grow_dst_boxes);
   }

   /* 
    * convert the lists to arrays; sort the indices
    * within each array; then build a set of adjacency
    * lists for the local boxes.
    */
   d_adj.resizeArray(dst_len);
   for (int x=0; x<dst_len; ++x) {
      int len = tmp_graph[x].getNumberItems();
      d_adj[x].resizeArray(len);

      int count = 0;
      for (tbox::List<int>::Iterator j(tmp_graph[x]); j; j++) {
        d_adj[x][count++] = j();
      }
   }

   sortIndices(d_adj);

   if (d_src_map.getSize()) {
      buildLocalOverlapArrays();
   }
   
   t_setup_total->stop();
}


/*
 * ***********************************************************************
 * When finding neighbors, the indices must be returned in
 * sorted order; this is because there are some undocumented
 * dependencies when constructing communication schedules.
 * The following functions are used for the sorting.
 * We also need to sort the node list during buildGraph.
 * ***********************************************************************
 */

template<int DIM> void BoxGraph<DIM>::sortIndices(tbox::Array< tbox::Array<int> > & graph)
{
   int size = graph.getSize();
   for (int j=0; j<size; ++j) {
      size_t len = graph[j].getSize();
      if (len) {
        int *tmp = graph[j].getPointer();
           qsort((void*)tmp, len, sizeof(int), BoxGraphUtilities<DIM>::qsortIntCompare);
      }
   }
}



template<int DIM> int BoxGraph<DIM>::qsortCompare(const void *v, const void *w)
{
      typename BoxGraph<DIM>::GraphNode *a = (typename BoxGraph<DIM>::GraphNode*) v;
      typename BoxGraph<DIM>::GraphNode *b = (typename BoxGraph<DIM>::GraphNode*) w;

      int i = a->coord;
      int j = b->coord;
      if (i > j) return 1;
      if (i < j) return -1;
      return 0;
}



template<int DIM> void BoxGraph<DIM>::sortNodeList(tbox::Array<BoxGraph<DIM>::GraphNode> & list, int len)
{
   typename BoxGraph<DIM>::GraphNode *tmp = new typename BoxGraph<DIM>::GraphNode[len];
   int h;
   for (h=0; h<len; ++h) tmp[h] = list[h];
   qsort((void*)tmp, len, sizeof(typename BoxGraph<DIM>::GraphNode), qsortCompare);
   for (h=0; h<len; ++h) list[h] = tmp[h];
   delete [] tmp;
}

/*
 * ************************************************************************
 * addEdges is is a private function called by buildGraph
 * ************************************************************************
 */

template<int DIM> void BoxGraph<DIM>::addEdges(
   tbox::Array<BoxGraph<DIM>::GraphNode> & src, int src_ct,
   tbox::Array<BoxGraph<DIM>::GraphNode> & dst, int dst_ct,
   tbox::Array<tbox::List<int> > &graph)
{
   for (int i=0; i<dst_ct; ++i) {
      typename BoxGraph<DIM>::GraphNode dst_node = dst[i];
      Box<DIM> dstbox = *(dst_node.box);
      int dst_idx = dst_node.idx;

      for (int j=0; j<src_ct; ++j) {
         typename BoxGraph<DIM>::GraphNode src_node = src[j];
         Box<DIM> srcbox = *(src_node.box);
         int src_idx = src_node.idx;

         if (dstbox.intersects(srcbox)) {

            //found edge (i,j); so need to insert j in adj(i),
            //if not previously inserted.  (The index could 
            //have been inserted if i and j both straddle 
            //the preceeding bin).

            //iterate over graph[dst_idx] to determine if
            //'src_idx' was previously inserted
            bool previously_discovered = false;
            for (tbox::List<int>::Iterator edge(graph[dst_idx]); edge; edge++) {
               if (edge() == src_idx) {
                  previously_discovered = true;
                  break;
               }
            }

            //insert the adjacency information if needed
            if (! previously_discovered) {
              graph[dst_idx].appendItem(src_idx);
            }
         }
      }
   }
}


template<int DIM> void BoxGraph<DIM>::buildLocalOverlapArrays()
{
   int myid = tbox::SAMRAI_MPI::getRank();
   int n = d_adj.getSize();
   d_adj_local.resizeArray(n);

   for (int j=0; j<n; ++j) {
      int count = 0;
      int len = d_adj[j].getSize();
      for (int k=0; k<len; ++k) {
         int idx = d_adj[j][k];
         if (d_src_map[idx] == myid) {
            ++count;
         }
      }

      d_adj_local[j].resizeArray(count);
      count = 0;
      for (int k=0; k<len; ++k) {
         int idx = d_adj[j][k];
         if (d_src_map[idx] == myid) {
            d_adj_local[j][count++] = idx;
         }
      }
   }
}


/*
 * ************************************************************************
 * findBinLimits is a private function called by buildGraph
 * ************************************************************************
 */

template<int DIM> void BoxGraph<DIM>::findBinLimits(
   int binNumberIN, int dim, tbox::Array<BoxGraph<DIM>::GraphNode> & nodesIN,
   int & firstOUT, int & lastOUT, int & rightEdgeOfBinOUT)
{
   int len = nodesIN.getSize();
   int bin_size = int(1 + sqrt((double)len));
   firstOUT = bin_size*binNumberIN;
   lastOUT = firstOUT + bin_size;
   lastOUT = (lastOUT > len) ? len : lastOUT;

   Box<DIM> box = *(nodesIN[lastOUT-1].box);
   rightEdgeOfBinOUT = box.lower(dim);
}

template<int DIM> void BoxGraph<DIM>::findSortDimension()
{
   //CASE 1: caller has explicitly specified a valid search direction
   if (d_sort_dimension >=0 && d_sort_dimension < DIM) {
      return;
   }

   //CASE 2: we will determine the search direction.  
   //        the following is a heuristic; the intent is to
   //        determine the direction in which the number of nodes
   //        that overlap bins is minimized; but optimally determining
   //        the direction would likely be costly---so we'll do
   //        something cheap.
   double rho = tbox::MathUtilities<double>::getMax();

   //Get bounding box for the lists
   int size_src = d_src_boxes.getNumberOfBoxes();
   int size_dst = d_dst_boxes.getNumberOfBoxes();
   int count = size_src + size_dst;
   Box<DIM> bbox;
   int i;

   for (i=0; i<size_src; ++i) {
      bbox += d_src_boxes(i);
   }
   for (i=0; i<size_dst; ++i) {
      bbox += d_dst_boxes(i);
   }

   Index<DIM> lower = bbox.lower();
   Index<DIM> upper = bbox.upper();
   for (int dir = 0; dir < DIM; ++dir) {
      int len = bbox.upper(dir) - bbox.lower(dir);
      double rho_tmp = (double)count/len;

      if (rho_tmp < rho) {
         d_sort_dimension = dir;
         rho = rho_tmp;
      }
   }
}

/*
 * ************************************************************************
 * Print functions; these are only needed for development/testing.
 * ************************************************************************
 */
template<int DIM> void BoxGraph<DIM>::print(tbox::Array<BoxGraph<DIM>::GraphNode> & list, int len)
{
    for (int i=0; i<len; ++i) {
      tbox::plog  << "   coord= " << list[i].coord
            << "; idx= " << list[i].idx
            << "; box= " << *(list[i].box)
            << "; isSrcBox= " << list[i].isSrcBox
            << std::endl;
    }
}

template<int DIM> void BoxGraph<DIM>::printGraph(std::ostream& os)
{
   int src_len = d_src_boxes.getNumberOfBoxes();
   int dst_len = d_dst_boxes.getNumberOfBoxes();
   int i;

   os << "\n===================================================\n";
   os << " START BoxGraph<DIM>::printGraph:\n";
   os << "===================================================\n";

   os << "\n---------------------------------------------------\n";
   os << "dst boxes (after growing):\n";
   for (i=0; i<dst_len; ++i) {
      os << "  dst box " << i << " = " << d_dst_boxes(i) << std::endl;
   }

   os << "\n---------------------------------------------------\n";
   os << "src boxes\n";
   for (i=0; i<src_len; ++i) {
      os << "  src box " << i << " = " << d_src_boxes(i) << std::endl;
   }

   os << "\n---------------------------------------------------\n";
   os << "Nabors of dst boxes:\n";
   for (i=0; i<dst_len; ++i) {
      os << "  dst box " << i << " = " << d_dst_boxes(i) << std::endl;
      os << "     src box nabors:\n";
      int adj_len = d_adj[i].getSize();
      for (int j=0; j<adj_len; ++j) {
         int idx = d_adj[i][j];

         if (idx >= src_len) {
             tbox::plog << "hysom's error!  idx= " << idx << ", should be < "
                  << src_len << std::endl;
         }

         os << "     idx= " << idx << std::endl;
      }
   }

   os << "---------------------------------------------------\n\n";
   os << "\n===================================================\n";
   os << " END BoxGraph<DIM>::printGraph:\n";
   os << "===================================================\n";
}


/*
 * ************************************************************************
 * removeIntersections
 * ************************************************************************
 */

template<int DIM> void BoxGraph<DIM>::removeIntersections(BoxList<DIM> &fragments)
{
   //This will hold boxes in "list" minus the intersection with "takeaway"
   fragments.clearItems();

   //Compare each box in "list (d_dst)" against all boxes in "takeaway (d_src)"
   //Note: this is opposite the approach of the current implementation

   int dst_len = d_dst_boxes.getNumberOfBoxes();

   for (int box_idx=0; box_idx<dst_len; ++box_idx) {
      Box<DIM> tryme = d_dst_boxes(box_idx);

      //Find a shorter list of boxes from "takeaway" that might intersect
      //with "tryme."  Recall that "takeaway" was the list with which
      //the BoxTop was instantiated; the next call returns a list of
      //boxes that in "takeaway," and are also nabors (i.e., overlap with)
      //the "tryme" box.
      //

      const tbox::Array<int> nabors = getSrcOverlapIndices(box_idx);
      int nabor_count = nabors.getSize();

      //if "tryme" doesn't intersect any boxes on takeaway
      //(i.e., if "nabors" is empty" then keep "tryme."
      if (nabor_count == 0) {
        fragments.appendItem(tryme);
      } 

      //else, burst up tryme against the nabors.
      else {

         BoxList<DIM> nabors_list;
         for (int j=0; j<nabor_count; ++j) {
            nabors_list.appendItem( d_src_boxes(nabors[j]) );
         }

         BoxList<DIM> tryme_burst(tryme);
         tryme_burst.removeIntersections(nabors_list);
         fragments.copyItems(tryme_burst);
      }
   }
}



}
}


#endif
