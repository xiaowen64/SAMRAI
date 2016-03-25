//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/clustering/BergerRigoutsos.C $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Class for Burger/Rigoutsos tagged cell clustering algorithm.
//

#ifndef included_mesh_BergerRigoutsos_C
#define included_mesh_BergerRigoutsos_C

#include "BergerRigoutsos.h"
#include "Index.h"
#include "BoxComm.h" 
#include "tbox/Database.h" 
#include "tbox/SAMRAI_MPI.h" 
#include "PatchHierarchy.h"
#include "BinaryTree.h"
#include "tbox/InputManager.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


namespace SAMRAI {
    namespace mesh {

/*
*************************************************************************
*                                                                       *
* Initialization for static data members.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer<tbox::Timer> BergerRigoutsos<DIM>::t_cluster;
template<int DIM>
tbox::Pointer<tbox::Timer> BergerRigoutsos<DIM>::t_globalize_boxes;
template<int DIM>
tbox::Pointer<tbox::Timer> BergerRigoutsos<DIM>::t_compute;
template<int DIM>
tbox::Pointer<tbox::Timer> BergerRigoutsos<DIM>::t_commwait;

template<int DIM> std::string
   BergerRigoutsos<DIM>::s_tag_reduce_method = "ORIGINAL";
template<int DIM> bool
   BergerRigoutsos<DIM>::s_tag_reduce_method_initialized = false;

/*
 * ************************************************************************
 *                                                                        *
 * Static function to set box intersection algorithm for schedules.       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void BergerRigoutsos<DIM>::setClusteringOption(
   const std::string& method)
{
   if (method.empty()) {

      if (tbox::SAMRAI_MPI::getNodes() == 1) {
         s_tag_reduce_method = "ORIGINAL";
      } else {
         s_tag_reduce_method = "BINARY_TREE";
      }

   } else {

      if ( !((method == "ORIGINAL") ||
             (method == "ORIG_FAST_REDUCE") ||
             (method == "COMMUNICATOR") ||
             (method == "BINARY_TREE")) ) {
         TBOX_ERROR("BergerRigoutsos<DIM>::setClusteringOption\n"
                    << "Given method string "
                    << method << " is invalid.\n Options are\n"
                    << "'ORIGINAL', 'ORIG_FAST_REDUCE', 'COMMUNICATOR', 'BINARY_TREE'."
                    << std::endl);
      }

      s_tag_reduce_method = method;
 
   }
}

/*
*******************************************************************************
*                                                                             *
* Default constructor and destructor for BergerRigoutsos<DIM>.               *
* The constructor queries the database to determine which version             *
* of findBoxesContainingTags will be used.                                    *
*                                                                             *
*******************************************************************************
*/

template<int DIM>  BergerRigoutsos<DIM>::BergerRigoutsos()
{
   if (!s_tag_reduce_method_initialized) {
      setClusteringOption();
      s_tag_reduce_method_initialized = true;
   }

   if ( t_cluster.isNull() ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( t_globalize_boxes.isNull() );
      TBOX_ASSERT( t_compute.isNull() );
      TBOX_ASSERT( t_commwait.isNull() );
#endif
      t_cluster = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::cluster");
      t_globalize_boxes = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::globalize_boxes");
      t_compute = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::compute");
      t_commwait = tbox::TimerManager::getManager()->
         getTimer("mesh::BergerRigoutsos::commwait");
      tbox::ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }
}

template<int DIM>  BergerRigoutsos<DIM>::~BergerRigoutsos()
{
}
 
/*
**************************************************************************
*                                                                        *
* Implementation of the Berger-Rigoutsos algorithm for creating a set    *
* of boxes boxes that cover all tags corresponding to the given          *
* patch level and descriptor index for the integer tag arrays.           *
* The tags must match the specified boolean value and reside within the  *
* bounding box. The boxes returned must be larger than the given         *
* minimum size, and the specified tolerances must be obeyed.             *
*                                                                        *
**************************************************************************
*/

template<int DIM> void BergerRigoutsos<DIM>::findBoxesContainingTags(
   hier::BoxList<DIM>& boxes, 
   const tbox::Pointer< hier::PatchLevel<DIM> > level, 
   const int index, 
   const int tag_val, 
   const hier::Box<DIM>& bound_box, 
   const hier::IntVector<DIM>& min_box,
   const double efficiency_tol, 
   const double combine_tol) const
{
   if (   s_tag_reduce_method == "ORIGINAL" 
       || s_tag_reduce_method == "ORIG_FAST_REDUCE") {
      findBoxesContainingTagsOriginal(boxes, level, index, tag_val,
                                      bound_box, min_box, 
                                      efficiency_tol, combine_tol);
   }

   else if (s_tag_reduce_method == "COMMUNICATOR") {
      findBoxesContainingTagsCommunicators(boxes, level, index, tag_val,
                                           bound_box, min_box, 
                                           efficiency_tol, combine_tol,
                                           0, tbox::SAMRAI_MPI::commWorld);
   }

   else if (s_tag_reduce_method == "BINARY_TREE") {
      findBoxesContainingTagsBinaryTree(boxes, level, index, tag_val,
                                        bound_box, min_box, 
                                        efficiency_tol, combine_tol, 0);
   }

   else {
      /*
       * should be impossible to be here!
       */
      TBOX_ERROR("unknown clustering algorithm: " 
                 << s_tag_reduce_method << std::endl);
   }
}


template<int DIM> void BergerRigoutsos<DIM>::findBoxesContainingTagsOriginal(
   hier::BoxList<DIM>& boxes, 
   const tbox::Pointer< hier::PatchLevel<DIM> > level, 
   const int index, 
   const int tag_val, 
   const hier::Box<DIM>& bound_box, 
   const hier::IntVector<DIM>& min_box,
   const double efficiency_tol, 
   const double combine_tol) const
{
   t_compute->start();
   HistogramBox<DIM> hist_box(bound_box);
   if (s_tag_reduce_method == "ORIG_FAST_REDUCE") {
      hist_box.setCommunicationMode(
         HistogramBox<DIM>::ORIG_FAST_REDUCE_MODE);
   }
   t_compute->stop();

   int num_tags = hist_box.computeTagHistogram(level, index, tag_val);

   if ( num_tags == 0 ) {
      boxes.clearItems();
   } 

   else {

      t_compute->start();
      hier::Box<DIM> tag_bound_box;
      hist_box.findBoundBoxForTags(tag_bound_box, min_box);

      int num_cells = tag_bound_box.size();
      double efficiency = ( num_cells == 0 ? 1.e0 :
                          ((double) num_tags)/((double) num_cells) );
      t_compute->stop();

      if (efficiency < efficiency_tol) {

         t_compute->start();
         hier::Box<DIM> box_lft;
         hier::Box<DIM> box_rgt;
         int is_split = splitTagBoundBox(box_lft, box_rgt, tag_bound_box, 
                                         hist_box, min_box);
         t_compute->stop();
         if ( is_split ) {
   
            /*
             * The bounding box "tag_bound_box" has been split into two
             * boxes, "box_lft" and "box_rgt".  Now, attempt to recursively
             * split these boxes further.
             */

            hier::BoxList<DIM> box_list_lft;
            hier::BoxList<DIM> box_list_rgt;
            findBoxesContainingTags(box_list_lft, 
                                    level, index, tag_val,
                                    box_lft, min_box,
                                    efficiency_tol, combine_tol);
            findBoxesContainingTags(box_list_rgt, 
                                    level, index, tag_val,
                                    box_rgt, min_box,
                                    efficiency_tol, combine_tol);
   
            /* 
             * If splitting the bounding box is good,
             * add the boxes to the list.
             */

            t_compute->start();
            if ( ((box_list_lft.getNumberItems() > 1) ||
                  (box_list_rgt.getNumberItems() > 1)) ||
                 ( (double) (box_list_lft.getFirstItem().size()
                           + box_list_rgt.getFirstItem().size())
                  < ((double) tag_bound_box.size())*combine_tol ) ) {

               boxes.catenateItems(box_list_lft);
               boxes.catenateItems(box_list_rgt);
            }
            t_compute->stop();

         }

      } 

      /*
       * If no good splitting is found, add bounding box to list.
       */

      if ( boxes.isEmpty() ) {
         t_compute->start();
         boxes.appendItem(tag_bound_box);
         t_compute->stop();
      }
   }
}


/*
**************************************************************************
*                                                                        *
* Reduces global communications by limiting recursive calls to those     *
* processors that have local patches that intersect with the bound box.  *
* This method explicitly manages communications using a specialized      *
* binary tree.                                                           *
*                                                                        *
**************************************************************************
*/

template<int DIM> void BergerRigoutsos<DIM>::findBoxesContainingTagsBinaryTree(
   hier::BoxList<DIM>& boxes, 
   const tbox::Pointer< hier::PatchLevel<DIM> > level, 
   const int index, 
   const int tag_val, 
   const hier::Box<DIM>& bound_box, 
   const hier::IntVector<DIM>& min_box,
   const double efficiency_tol, 
   const double combine_tol,
   int recurse_level) const
{
   if ( recurse_level == 0 ) t_cluster->start();

   t_compute->start();
   tbox::Pointer< hier::BinaryTree<DIM> > tree = level->getBinaryTree();

   /*
    * "participates" is true if this processor has at least one local
    * patch that intersects with bound_box.
    */
   const bool part = tree->participates(tbox::SAMRAI_MPI::getRank(), bound_box);
   t_compute->stop();
   if ( part ) {

      t_compute->start();
      HistogramBox<DIM> hist_box(bound_box);
      hist_box.setCommunicationMode(HistogramBox<DIM>::BINARY_TREE_MODE);
      int num_tags = hist_box.computeTagHistogram(level, index, tag_val);
      t_compute->stop();

      if ( num_tags == 0 ) {
         if (tbox::SAMRAI_MPI::getRank() == 0) {
            boxes.clearItems();
         }
      }

      else {

         ++recurse_level;

         int need_to_recursively_split = 0;
         hier::Box<DIM> box_lft;
         hier::Box<DIM> box_rgt;
         hier::Box<DIM> tag_bound_box;

         /*
          * only root can decide if recursive splitting is necessary,
          * since only root has the complete histogram
          */
         if (tbox::SAMRAI_MPI::getRank() == 0) {
            t_compute->start();
            hist_box.findBoundBoxForTags(tag_bound_box, min_box);
            int num_cells = tag_bound_box.size();
            double efficiency = ( num_cells == 0 ? 1.e0 :
                                ((double) num_tags)/((double) num_cells) );

            if (efficiency < efficiency_tol) {
               need_to_recursively_split = splitTagBoundBox(
                  box_lft, box_rgt, tag_bound_box, hist_box, min_box);
            }
            t_compute->stop();
         }

         /*
          * root must tell the other processors if recursive splitting
          * is needed.
          */
         t_commwait->start();
         tree->partialBcast(bound_box, need_to_recursively_split);
         t_commwait->stop();

         if (need_to_recursively_split) {

            /*
             * only P_0 knows the boxes, so they must be broadcast
             * to others.
             */
            t_commwait->start();
            tree->partialBcast(bound_box, box_lft);
            tree->partialBcast(bound_box, box_rgt);
            t_commwait->stop();

            t_compute->start();
            hier::BoxList<DIM> box_list_rgt;
            hier::BoxList<DIM> box_list_lft;
            const bool lpart = tree->participates(tbox::SAMRAI_MPI::getRank(), box_lft);
            const bool rpart = tree->participates(tbox::SAMRAI_MPI::getRank(), box_rgt);
            t_compute->stop();

            /*
             * recursive call to split box_lft
             */
            if (lpart) {
               findBoxesContainingTagsBinaryTree(
                                       box_list_lft, 
                                       level, index, tag_val,
                                       box_lft, min_box,
                                       efficiency_tol, combine_tol,
                                       recurse_level);
            }

            /*
             * recursive call to split box_rgt
             */
            if (rpart) {
               findBoxesContainingTagsBinaryTree(
                                       box_list_rgt, 
                                       level, index, tag_val,
                                       box_rgt, min_box,
                                       efficiency_tol, combine_tol,
                                       recurse_level);
            }
   
            /* 
             * If splitting the bounding box is good,
             * P_0 adds the boxes to the list.
             */

            t_compute->start();
            if (tbox::SAMRAI_MPI::getRank() == 0) {
               if ( ((box_list_lft.getNumberItems() > 1) ||
                     (box_list_rgt.getNumberItems() > 1)) ||
                    ( (double) (box_list_lft.getFirstItem().size()
                              + box_list_rgt.getFirstItem().size())
                     < ((double) tag_bound_box.size())*combine_tol ) ) {

                  boxes.catenateItems(box_list_lft);
                  boxes.catenateItems(box_list_rgt);
               }
            } //if (tbox::SAMRAI_MPI::getRank() == 0) 
            t_commwait->stop();
         } //if (need_to_recursively_split)

         /*
          * If no good splitting is found, P_0 adds the bounding box to the list.
          */
         if (tbox::SAMRAI_MPI::getRank() == 0) {
            if ( boxes.isEmpty() ) {
               boxes.appendItem(tag_bound_box);
            }
         }

         --recurse_level;
      }

   } //if (participates)


   /*
    * If we're exiting from the outermost recursion level, then the
    * list of boxes is complete, but only P_0 has a copy; therefore,
    * we must broadcast the list to all other processors.
    */
   if (recurse_level == 0) {
      t_cluster->stop();
      t_globalize_boxes->start();
      hier::BoxComm<DIM>::bcastBoxList(boxes);
      t_globalize_boxes->stop();
   }
}


/*
**************************************************************************
*                                                                        *
* Reduces global communications by limiting recursive calls to those     *
* processors that have local patches that intersect with the bound box.  *
* This method forms new tbox::MPI communicators to handle the communications   *
* in the subsets of processors.                                          *
*                                                                        *
**************************************************************************
*/

template<int DIM> void BergerRigoutsos<DIM>::findBoxesContainingTagsCommunicators(
   hier::BoxList<DIM>& boxes, 
   const tbox::Pointer< hier::PatchLevel<DIM> > level, 
   const int index, 
   const int tag_val, 
   const hier::Box<DIM>& bound_box, 
   const hier::IntVector<DIM>& min_box,
   const double efficiency_tol, 
   const double combine_tol,
   int recurse_level,
   tbox::SAMRAI_MPI::comm comm) const
{
   tbox::Pointer< hier::BinaryTree<DIM> > tree = level->getBinaryTree();

   /*
    * Form a new communicator, "comm_bound_box," that only contains
    * processors that have at least one patch that intersects bound_box.
    * Note that we only need to do this at the outermost recursion
    * layer; this is because we also form new communicators before
    * the left-box and right-box recursive calls; doing so again
    * here would duplicate that effort.
    * 
    */
   tbox::SAMRAI_MPI::comm comm_bound_box;
   tbox::SAMRAI_MPI::group group_bound_box;

   tree->buildParticipatingCommunicator(bound_box, comm,
                                        group_bound_box,
                                        comm_bound_box);

   /*
    * "tree->participates" is true if this processor has at least
    * one local patch that intersects with bound_box.
    */
   if ( tree->participates(tbox::SAMRAI_MPI::getRank(), bound_box) ) {

      HistogramBox<DIM> hist_box(bound_box);
      hist_box.setCommunicationMode(HistogramBox<DIM>::COMMUNICATOR_MODE);
      hist_box.setCommunicator(comm_bound_box);
      int num_tags = hist_box.computeTagHistogram(level, index, tag_val);


      if ( num_tags == 0) {
         if (tbox::SAMRAI_MPI::getRank() == 0) {
            boxes.clearItems();
         }
      }
   
      else {
   
         ++recurse_level;
   
         int need_to_recursively_split = 0;
         hier::Box<DIM> box_lft;
         hier::Box<DIM> box_rgt;
         hier::Box<DIM> tag_bound_box;
   
         /*
          * only root can decide if recursive splitting is necessary,
          * since only root has the complete histogram
          */
         if (tbox::SAMRAI_MPI::getRank() == 0) {
            hist_box.findBoundBoxForTags(tag_bound_box, min_box);
            int num_cells = tag_bound_box.size();
            double efficiency = ( num_cells == 0 ? 1.e0 :
                                ((double) num_tags)/((double) num_cells) );

            if (efficiency < efficiency_tol) {
               need_to_recursively_split = splitTagBoundBox(
                  box_lft, box_rgt, tag_bound_box, hist_box, min_box);
            }
         }

         /*
          * root must tell the other processors if recursive splitting
          * is needed.
          */

         tbox::SAMRAI_MPI::comm keep_me = tbox::SAMRAI_MPI::getCommunicator();
         tbox::SAMRAI_MPI::setCommunicator(comm_bound_box);
         need_to_recursively_split = 
            tbox::SAMRAI_MPI::bcast(need_to_recursively_split, 0);
         tbox::SAMRAI_MPI::setCommunicator(keep_me);

         if (need_to_recursively_split) {
   
            /*
             * only P_0 knows the boxes, so they must be broadcast
             * to others.
             */
            hier::BoxComm<DIM>::bcastBox(box_lft, 0, comm_bound_box);
            hier::BoxComm<DIM>::bcastBox(box_rgt, 0, comm_bound_box);
   
            hier::BoxList<DIM> box_list_rgt;
            hier::BoxList<DIM> box_list_lft;
   
            tbox::SAMRAI_MPI::comm comm_lft_box, comm_rgt_box;
            tbox::SAMRAI_MPI::group group_lft_box, group_rgt_box;

            /*
             * recursive call to split box_lft
             */
            tree->buildParticipatingCommunicator(box_lft, comm_bound_box,
                                                 group_lft_box, comm_lft_box);

            if ( tree->participates(tbox::SAMRAI_MPI::getRank(), box_lft) ) {
               findBoxesContainingTagsCommunicators(
                                       box_list_lft, 
                                       level, index, tag_val,
                                       box_lft, min_box,
                                       efficiency_tol, combine_tol, 
                                       recurse_level, comm_lft_box);
            }

            if (comm_lft_box != tbox::SAMRAI_MPI::commNull) {
#ifdef HAVE_MPI	      
               MPI_Comm_free(&comm_lft_box);
               MPI_Group_free(&group_lft_box);
#endif
            }

            /*
             * recursive call to split box_rgt
             */
            tree->buildParticipatingCommunicator(box_rgt, comm_bound_box,
                                                 group_rgt_box, comm_rgt_box);
   
            if ( tree->participates(tbox::SAMRAI_MPI::getRank(), box_rgt) ) {
               findBoxesContainingTagsCommunicators(
                                       box_list_rgt, 
                                       level, index, tag_val,
                                       box_rgt, min_box,
                                       efficiency_tol, combine_tol, 
                                       recurse_level, comm_rgt_box);
            }

            if (comm_rgt_box != tbox::SAMRAI_MPI::commNull) {
#ifdef HAVE_MPI	      
               MPI_Comm_free(&comm_rgt_box);
               MPI_Group_free(&group_rgt_box);
#endif
            }

            /* 
             * If splitting the bounding box is good,
             * P_0 adds the boxes to the list.
             */
            if (tbox::SAMRAI_MPI::getRank() == 0) {
               if ( ((box_list_lft.getNumberItems() > 1) ||
                     (box_list_rgt.getNumberItems() > 1)) ||
                    ( (double) (box_list_lft.getFirstItem().size()
                              + box_list_rgt.getFirstItem().size())
                     < ((double) tag_bound_box.size())*combine_tol ) ) {
   
                  boxes.catenateItems(box_list_lft);
                  boxes.catenateItems(box_list_rgt);
               }
            } //if (tbox::SAMRAI_MPI::getRank() == 0) 
   
         } //if (need_to_recursively_split)
   
         /*
          * If no good splitting is found, P_0 adds the bounding box
          * to the list.
          */
   
         if (tbox::SAMRAI_MPI::getRank() == 0) {
            if ( boxes.isEmpty() ) {
               boxes.appendItem(tag_bound_box);
            }
         }
   
         --recurse_level;
      }
   } //if (participates)

//   if (recurse_level > 0) {
      if (comm_bound_box != tbox::SAMRAI_MPI::commNull) {
#ifdef HAVE_MPI	      
         MPI_Comm_free(&comm_bound_box);
         MPI_Group_free(&group_bound_box);
#endif
      }
//   }

   /*
    * If we're exiting from the outermost recursion level, then the
    * list of boxes is complete, but only P_0 has a copy; therefore,
    * we must broadcast the list to all other processors.
    */
   if (recurse_level == 0) {
      hier::BoxComm<DIM>::bcastBoxList(boxes);
   }

}



/*
*************************************************************************
*                                                                       *
* Attempt to split the box bounding a collection of tagged cells into   *
* two boxes. If an appropriate splitting is found, the two smaller      *
* boxes are returned and the return value of the function is true.      *
* Otherwise, false is returned. Note that the bounding box must be      *
* contained in the histogram box and that the two splitting boxes must  *
* be larger than some smallest box size.                                *
*                                                                       *
* Note that it is assumed that box indices are cell indices.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BergerRigoutsos<DIM>::splitTagBoundBox(
   hier::Box<DIM>& box_lft, 
   hier::Box<DIM>& box_rgt,
   const hier::Box<DIM>& bound_box, 
   const HistogramBox<DIM>& hist_box, 
   const hier::IntVector<DIM>& min_box) const 
{

   int cut_pt = -(tbox::MathUtilities<int>::getMax());
   int tmp_dim = -1; 
   int id = -1;
   hier::IntVector<DIM> num_cells(bound_box.numberCells()); 
   hier::Index<DIM> box_lo(bound_box.lower()); 
   hier::Index<DIM> box_hi(bound_box.upper());

   /*
    * Sort the bound box dimensions from largest to smallest.
    */
   hier::IntVector<DIM> sorted;
   for (id = 0; id < DIM; id++) {
      sorted(id) = id;
   }
   for (int id0 = 0; id0 < DIM-1; id0++) {
      for (int id1 = id0+1; id1 < DIM; id1++) {
         if ( num_cells(sorted(id0)) < num_cells(sorted(id1)) ) {
            tmp_dim = sorted(id0);
            sorted(id0) = sorted(id1);
            sorted(id1) = tmp_dim; 
         }
      }
   }

   /*
    * Determine number of coordinate directions in bounding box
    * that are splittable according to the minimum box size restriction.
    */

   int nsplit;
   for (nsplit = 0; nsplit < DIM; nsplit++) {
      tmp_dim = sorted(nsplit);
// Which test??
      if ( num_cells(tmp_dim) < 2*min_box(tmp_dim) ) {
         break;
      }
//    if ( (num_cells(tmp_dim) < 2*min_box(tmp_dim)) ||
//         (num_cells(tmp_dim) < num_cells(sorted(0))/2) ) {
//       break;
//    }
   }
   if (nsplit == 0) return(false);

   /*
    * Attempt to split box at a zero interior point in the histogram.
    * Check each splittable direction, from largest to smallest, until
    * zero point found.
    */

   for (id = 0; id < nsplit; id++) {
      tmp_dim = sorted(id);
      if ( findZeroCutPoint(cut_pt, 
                            tmp_dim, box_lo(tmp_dim), box_hi(tmp_dim), 
                            hist_box, min_box(tmp_dim)) ) {
         break;
      }
   }

   /*
    * If no zero point found, try Laplacian on longest side of bound box.
    */

   if (id == nsplit) {
      tmp_dim = sorted(0);
      cutAtLaplacian(cut_pt, 
                     tmp_dim, box_lo(tmp_dim), box_hi(tmp_dim), 
                     hist_box, min_box(tmp_dim));
   }

   /*
    * Split bound box at cut_pt; tmp_dim is splitting dimension.
    */
   hier::Index<DIM> lft_hi(box_hi); 
   hier::Index<DIM> rgt_lo(box_lo); 
   lft_hi(tmp_dim) = cut_pt;
   rgt_lo(tmp_dim) = cut_pt + 1;

   box_lft = hier::Box<DIM>(box_lo, lft_hi);  
   box_rgt = hier::Box<DIM>(rgt_lo, box_hi);

   return(true);
}

/*
*************************************************************************
*                                                                       *
* Attempt to find a zero histogram value near the middle of the index   *
* interval (lo, hi) in the given coordinate direction. Note that the    *
* cut_pt is kept more than a minimium distance from the endpoints of    *
* of the index interval. Since box indices are cell-centered, the cut   *
* point value corresponds to the right edge of the cell whose index     *
* is equal to the cut point.                                            *
*                                                                       *
* Note that it is assumed that box indices are cell indices.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BergerRigoutsos<DIM>::findZeroCutPoint(
   int& cut_pt, 
   const int id, 
   const int lo, 
   const int hi, 
   const HistogramBox<DIM>& hist_box, 
   const int min_size) const
{
   const int cut_lo = lo + min_size - 1;
   const int cut_hi = hi - min_size;
   const int cut_mid = (lo + hi)/2;

   for (int ic = 0; ((cut_mid-ic>=cut_lo) && (cut_mid+ic<=cut_hi)); ic++) {
      if (hist_box.histogramElt(id, cut_mid-ic) == 0) {
         cut_pt = cut_mid-ic;
         return(true);
      }
      if (hist_box.histogramElt(id, cut_mid+ic) == 0) {
         cut_pt = cut_mid+ic;
         return(true);
      }
   }
   
   return(false);
}

/*
***************************************************************************
*                                                                         *
* Attempt to find a point in the given coordinate direction near an       *
* inflection point in the histogram for that direction. Note that the     *
* cut point is kept more than a minimium distance from the endpoints      *
* of the index interval (lo, hi).  Also, the box must have at least       *
* three cells along a side to apply the Laplacian test.  If no            *
* inflection point is found, the mid-point of the interval is returned    *
* as the cut point.                                                       *
*                                                                         *
* Note that it is assumed that box indices are cell indices.              *
*                                                                         *
***************************************************************************
*/

template<int DIM> void BergerRigoutsos<DIM>::cutAtLaplacian(
   int& cut_pt, 
   const int id,
   const int lo, 
   const int hi, 
   const HistogramBox<DIM>& hist_box,
   const int min_size) const
{
   int loc_zero = (lo + hi)/2;

   if ( (hi - lo + 1) >= 3 ) {

      int max_zero = 0;
      const int infpt_lo = tbox::MathUtilities<int>::Max(lo + min_size - 1, 
                                                         lo + 1);  
      const int infpt_hi = tbox::MathUtilities<int>::Min(hi - min_size, 
                                                         hi - 2);  

      int last_lap = hist_box.histogramElt(id, infpt_lo - 1) 
                     - 2 * hist_box.histogramElt(id, infpt_lo)
                     + hist_box.histogramElt(id, infpt_lo + 1);
  
      for (int ic = infpt_lo+1; ic <= infpt_hi+1; ic++) {
         int new_lap = hist_box.histogramElt(id, ic - 1)
                       - 2 * hist_box.histogramElt(id, ic)
                       + hist_box.histogramElt(id, ic + 1);

         if ( ((new_lap < 0) && (last_lap >= 0)) ||
              ((new_lap >= 0) && (last_lap < 0)) ) {

            int delta = new_lap - last_lap;

            if ( delta < 0 ) {
               delta = -delta; 
            }

            if ( delta > max_zero ) {
               loc_zero = ic - 1;
               max_zero = delta;
            }
         }
     
         last_lap = new_lap;
      }

   }

   cut_pt = loc_zero;
}




/*
***************************************************************************
*                                                                         *
* Release static timers.  To be called by shutdown registry to make sure  *
* memory for timers does not leak.                                        *
*                                                                         *
***************************************************************************
*/
template<int DIM>
void BergerRigoutsos<DIM>::freeTimers()
{
   t_cluster.setNull();
   t_globalize_boxes.setNull();
   t_compute.setNull();
   t_commwait.setNull();
   return;
}

}
}
#endif
