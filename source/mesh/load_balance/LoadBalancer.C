//
// File:        LoadBalancer.C
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 454 $
// Modified:    $Date: 2005-06-16 12:44:09 -0700 (Thu, 16 Jun 2005) $
// Description: Load balance routines for uniform and non-uniform workloads.
//

#ifndef included_mesh_LoadBalancer_C
#define included_mesh_LoadBalancer_C

#include "LoadBalancer.h"

#include <stdlib.h>
#include <fstream>
using namespace std;
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif
#include "BoxUtilities.h"
#include "BoxComm.h"
#include "PatchDescriptor.h"
#include "VariableDatabase.h"
#include "BalanceUtilities.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellDoubleConstantRefine.h"
#include "tbox/Array.h"
#include "tbox/List.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"

namespace SAMRAI {
    namespace mesh {

/*
*************************************************************************
*                                                                       *
* Constructors and destructor for LoadBalancer<DIM>.                   *
*                                                                       *
*************************************************************************
*/

template<int DIM>  LoadBalancer<DIM>::LoadBalancer(
   const string& name,
   tbox::Pointer<tbox::Database> input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif

   d_object_name = name;

   d_processor_layout_specified = false;
   d_processor_layout = hier::IntVector<DIM>(0);

   d_master_workload_data_id = -1;

   d_master_max_workload_factor = 1.0;

   d_master_bin_pack_method = "SPATIAL";

   d_workload_data_id.resizeArray(0);
   d_max_workload_factor.resizeArray(0);
   d_bin_pack_method.resizeArray(0);

   d_ignore_level_box_union_is_single_box = false;

   getFromInput(input_db);

   t_load_balance_boxes = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::loadBalanceBoxes()");
   t_load_balance_boxes_remove_intersection = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::loadBalanceBoxes()_remove_intersection");
   t_bin_pack_boxes = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::binPackBoxes()");
   t_bin_pack_boxes_sort = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::binPackBoxes()_sort");
   t_bin_pack_boxes_pack = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::binPackBoxes()_pack");
   t_chop_boxes = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::chop_boxes");
   
}

template<int DIM>  LoadBalancer<DIM>::LoadBalancer(
   tbox::Pointer<tbox::Database> input_db)
{
   d_object_name = "LoadBalancer";

   d_processor_layout_specified = false;
   d_processor_layout = hier::IntVector<DIM>(0);

   d_master_workload_data_id = -1;

   d_master_max_workload_factor = 1.0;

   d_master_bin_pack_method = "SPATIAL";

   d_workload_data_id.resizeArray(0);
   d_max_workload_factor.resizeArray(0);
   d_bin_pack_method.resizeArray(0);

   getFromInput(input_db);

   t_load_balance_boxes = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::loadBalanceBoxes()");
   t_load_balance_boxes_remove_intersection = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::loadBalanceBoxes()_remove_intersection");
   t_bin_pack_boxes = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::binPackBoxes()");
   t_bin_pack_boxes_sort = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::binPackBoxes()_sort");
   t_bin_pack_boxes_pack = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::binPackBoxes()_pack");
   t_chop_boxes = tbox::TimerManager::getManager() ->
      getTimer("mesh::LoadBalancer::chop_boxes");

}

template<int DIM>  LoadBalancer<DIM>::~LoadBalancer()
{
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to get/set load balancing parameters.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool LoadBalancer<DIM>::getLoadBalanceDependsOnPatchData(int level_number) const
{
   return( getWorkloadDataId(level_number) < 0 ? false : true );
}

template<int DIM> void LoadBalancer<DIM>::setMaxWorkloadFactor(
   double factor,
   int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(factor > 0.0);
#endif
   if (level_number >= 0) {
      int asize = d_max_workload_factor.getSize();
      if (asize < level_number+1) {
         d_max_workload_factor.resizeArray(level_number+1);
         for (int i = asize; i < level_number-1; i++) {
            d_max_workload_factor[i] =
               d_master_max_workload_factor;
         }
         d_max_workload_factor[level_number] = factor;
      }
   } else {
      d_master_max_workload_factor = factor;
      for (int ln = 0; ln < d_max_workload_factor.getSize(); ln++) {
         d_max_workload_factor[ln] = d_master_max_workload_factor;
      } 
   }
}

template<int DIM> void LoadBalancer<DIM>::setWorkloadPatchDataIndex(
   int data_id,
   int level_number)
{
   tbox::Pointer< pdat::CellDataFactory<DIM,double> > datafact = 
      hier::VariableDatabase<DIM>::getDatabase()->getPatchDescriptor()->
         getPatchDataFactory(data_id);
   if (datafact.isNull()) {
      TBOX_ERROR(d_object_name << " error: "
                 << "\n   data_id " << data_id << " passed to "
                 << "setWorkloadPatchDataIndex()"
                 << " does not refer to cell-centered double patch data. " << endl);
   }

   if (level_number >= 0) {
      int asize = d_workload_data_id.getSize();
      if (asize < level_number+1) {
         d_workload_data_id.resizeArray(level_number+1);
         for (int i = asize; i < level_number-1; i++) {
            d_workload_data_id[i] =
               d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = data_id;
      }
   } else {
      d_master_workload_data_id = data_id;
      for (int ln = 0; ln < d_workload_data_id.getSize(); ln++) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}


template<int DIM> void LoadBalancer<DIM>::setUniformWorkload(
   int level_number)
{
   if (level_number >= 0) {
      int asize = d_workload_data_id.getSize();
      if (asize < level_number+1) {
         d_workload_data_id.resizeArray(level_number+1);
         for (int i = asize; i < level_number-1; i++) {
            d_workload_data_id[i] =
               d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = -1;
      }
   } else {
      d_master_workload_data_id = -1;
      for (int ln = 0; ln < d_workload_data_id.getSize(); ln++) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}

template<int DIM> void LoadBalancer<DIM>::setBinPackMethod(
   const string& method,
   int level_number)
{

   if ( !(method == "GREEDY" ||
          method == "SPATIAL") ) {
      TBOX_ERROR(d_object_name << " error: "
                 << "\n   String " << method << " passed to setBinPackMethod()"
                 << " is not a valid method string identifier." << endl);

   }
   
   if (level_number >= 0) {
      int asize = d_bin_pack_method.getSize(); 
      if (asize < level_number+1) {
         d_bin_pack_method.resizeArray(level_number+1);
         for (int i = asize; i < level_number-1; i++) {
            d_bin_pack_method[i] = d_master_bin_pack_method;
         }  
         d_bin_pack_method[level_number] = method;
      }
   } else {
      d_master_bin_pack_method = method;
      for (int ln = 0; ln < d_bin_pack_method.getSize(); ln++) {
         d_bin_pack_method[ln] = d_master_bin_pack_method;
      }
   }
}

template<int DIM> void 
LoadBalancer<DIM>::setIgnoreLevelDomainIsSingleBox(bool flag)
{
   d_ignore_level_box_union_is_single_box = flag;
}

/*
*************************************************************************
*                                                                       *
* This main load balance routine performs either uniform or             *
* non-uniform load balancing operations on the given level depending    *
* on user specifications.   In either case, the goal is to produce      *
* a set of boxes and a mapping of those boxes to processors so that     *
* the workload on each processor is close to the average workload.      *
* The average workload is the total computational workload divided      *
* by the number of processors.  In the uniform load balance case        *
* (default), the workload is the number of cells in each box.  In the   *
* non-uniform case, the workload is computed using weight data on the   *
* grid hierarchy (i.e., a cell-centered double array on each patch.     *
*                                                                       *
* Typically, any box whose load is larger than the average is chopped.  *
* A user can prescribe a parameter (the 'max workload factor') to alter *
* average load used in this computation. Chopping is done using the     *
* BalanceUtilities<DIM>::recursiveBisection()) method which is similar *
* to the Berger-Rigoutsos algorithm.                                    *
*                                                                       *
* Once the boxes are chopped into a collection os smaller boxes, they   *
* are assigned to processors by a bin packing algorithm.                *
*                                                                       *
* The algorithm is summarized as follows:
*                                                                       *
* 1) Compute the estimated workload associated with each box.  In the   *
*    uniform workload case, the load in each box region is the number   *
*    of cells in the region.  Otherwise, the workload is computed using *
*    patch data defined by the d_workload_data_id array set by the user.*
*                                                                       *
* 2) Compute the maximum workload allowed on any box.  This quantity is *
*    by default the total workload divided by the number of processors. *
*    The user may provide a maximum workload factor, either through the *
*    input file or through a member function, which can alter the       *
*    average workload used in this computation.                         *
*                                                                       *
* 3) Chop each box whose workload is more than the max allowed into a   *
*    smaller set of boxes.                                              *
*                                                                       *
* 4) Check constraints placed on the boxes by the problem - i.e.        *
*    verify boxes are within the maximum and minimum box size           *
*    constraints and maintain a correct cut factor.                     *
*                                                                       *
* 5) Sort boxes largest to smallest and form an array.  Also form an    *
*    array of the workloads associated with each box.                   *
*                                                                       *
* 6) Use a bin packing procedure to construct a processor mapping for   *
*    the set of boxes.                                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::loadBalanceBoxes(
   hier::BoxArray<DIM>& out_boxes,
   hier::ProcessorMapping& mapping,
   const hier::BoxList<DIM>& in_boxes,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int level_number,
   const hier::BoxArray<DIM>& physical_domain,
   const hier::IntVector<DIM>& ratio_to_hierarchy_level_zero, 
   const hier::IntVector<DIM>& min_size,
   const hier::IntVector<DIM>& max_size,
   const hier::IntVector<DIM>& cut_factor,
   const hier::IntVector<DIM>& bad_interval) const
{
   t_load_balance_boxes->start();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert(level_number >= 0);
   assert(physical_domain.getNumberOfBoxes() > 0);
   assert(min_size > hier::IntVector<DIM>(0));   
   assert(max_size >= min_size);   
   assert(cut_factor > hier::IntVector<DIM>(0));   
   assert(bad_interval >= hier::IntVector<DIM>(0));   
#if 0
   /*
    * Additional assertion check for debugging - make sure given boxes
    * satisfy min size, cut_factor, bad_interval, and physical domain
    * constraints.  If the supplied boxes don't satisfy these constraints,
    * the load balancer cannot generate boxes that will.
    */
   for (hier::BoxList<DIM>::Iterator bi(in_boxes); bi; bi++) {
      hier::BoxUtilities<DIM>::checkBoxConstraints(bi(),
                                              min_size,
                                              cut_factor,
                                              bad_interval,
                                              physical_domain);

   }
#endif
#endif

   /*
    * If uniform load balancing is used and the level domain can be
    * expressed as a single box, we can construct an optimal box 
    * layout across processors without more involved chopping operations.
    *
    * Otherwise, we chop each box individually to construct the new array
    * of boxes and associated array of workloads based on either uniform
    * or nonuniform workload estimates.
    */

   int wrk_indx = getWorkloadDataId(level_number);

   tbox::Array<double> workloads;

   if ( (wrk_indx < 0) || (hierarchy->getNumberLevels() == 0) ) {

      if ( !d_ignore_level_box_union_is_single_box ) {
         hier::Box<DIM> bbox = in_boxes.getBoundingBox();
         hier::BoxList<DIM> difference(bbox);
         t_load_balance_boxes_remove_intersection->start();
         difference.removeIntersections(in_boxes);
         t_load_balance_boxes_remove_intersection->stop();

         if (difference.isEmpty()) {

            t_chop_boxes->start();
            chopUniformSingleBox(out_boxes,
                                 workloads,
                                 bbox,
                                 min_size,
                                 max_size,
                                 cut_factor,
                                 bad_interval,
                                 physical_domain);
            t_chop_boxes->stop();

         } else {

            t_chop_boxes->start();
            chopBoxesWithUniformWorkload(out_boxes,
                                         workloads,
                                         in_boxes,
                                         hierarchy,
                                         level_number,
                                         min_size,
                                         max_size,
                                         cut_factor,
                                         bad_interval,
                                         physical_domain); 
            t_chop_boxes->stop();

         }
      }
      else {
         t_chop_boxes->start();
         chopBoxesWithUniformWorkload(out_boxes,
                                      workloads,
                                      in_boxes,
                                      hierarchy,
                                      level_number,
                                      min_size,
                                      max_size,
                                      cut_factor,
                                      bad_interval,
                                      physical_domain); 
         t_chop_boxes->stop();
      }

   } else {

      t_chop_boxes->start();
      chopBoxesWithNonuniformWorkload(out_boxes,
                                      workloads, 
                                      in_boxes,
                                      hierarchy,
                                      level_number,
                                      ratio_to_hierarchy_level_zero,
                                      wrk_indx,
                                      min_size, 
                                      max_size, 
                                      cut_factor,
                                      bad_interval,
                                      physical_domain);
      t_chop_boxes->stop();

   }

#if 0
   /*
    * Assertion check for additional debugging - make sure new boxes 
    * satisfy min size, cut_factor, bad_interval, and physical domain 
    * constraints.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   const int nboxes = out_boxes.getNumberOfBoxes();
   for (int ib = 0; ib < nboxes; ib++) {
      hier::BoxUtilities<DIM>::checkBoxConstraints(out_boxes.getBox(ib),
                                              min_size,
                                              cut_factor,
                                              bad_interval,
                                              physical_domain);

   }
#endif
#endif

   /*
    * Generate mapping of boxes to processors using workload estimates.
    * Note boxes in array may be reordered during this operation.
    */

   binPackBoxes(out_boxes,
                mapping,
                workloads,
                getBinPackMethod(level_number)); 

   t_load_balance_boxes->stop();
}

/*
*************************************************************************
*                                                                       *
* Private function that chops a single box into a set of boxes          *
* that will map to the array of processors as uniformly as possible.    *
* The routine assumes a spatially-uniform workload distribution.        *
* Note that no error checking is done.                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::chopUniformSingleBox(
   hier::BoxArray<DIM>& out_boxes,
   tbox::Array<double>& out_workloads,
   const hier::Box<DIM>& in_box,
   const hier::IntVector<DIM>& min_size,
   const hier::IntVector<DIM>& max_size,
   const hier::IntVector<DIM>& cut_factor,
   const hier::IntVector<DIM>& bad_interval,
   const hier::BoxArray<DIM>& physical_domain) const
{

   /*
    * Determine processor layout that corresponds to box dimensions.
    */

   hier::IntVector<DIM> processor_distribution;
   if (d_processor_layout_specified) {
      processor_distribution = d_processor_layout;
   } else {
      BalanceUtilities<DIM>::computeDomainDependentProcessorLayout(
         processor_distribution,
         tbox::MPI::getNodes(),
         in_box);
   }

   /*
    * The ideal box size will be the dimensions of the input box divided
    * by the number of processors in each direction.  Compute this
    * ideal size and then adjust as necessary to fit within min/max size
    * constraints.
    */

   hier::IntVector<DIM> ideal_box_size;
   for (int i = 0; i < DIM; i++) {
      ideal_box_size(i) = (in_box.numberCells(i)/processor_distribution(i));
      ideal_box_size(i) = (ideal_box_size(i) > max_size(i) ?
                           max_size(i) : ideal_box_size(i));
      ideal_box_size(i) = (ideal_box_size(i) < min_size(i) ?
                           min_size(i) : ideal_box_size(i));
   }

   /*
    * Chop the single input box into a set of smaller boxes using the
    * ideal_box_size as the maximum size of each of the smaller boxes.
    */

   hier::BoxList<DIM> tmp_box_list(in_box);   

   hier::BoxUtilities<DIM>::chopBoxes(tmp_box_list,
                                 ideal_box_size,
                                 min_size,
                                 cut_factor,
                                 bad_interval,
                                 physical_domain);

   /*
    * Set output box array to list of chopped boxes and set workload array.
    */

   out_boxes = tmp_box_list;
   
   const int nboxes = out_boxes.getNumberOfBoxes();
   out_workloads.resizeArray(nboxes);
   for (int ib = 0; ib < nboxes; ib++) {
      out_workloads[ib] = (double)(out_boxes.getBox(ib).size());
   }

}

/*
*************************************************************************
*                                                                       *
* Private function that chops boxes on a list into another list of      *
* boxes all of which have approximately or less than an average         *
* workload estimate.  The routine assumes a spatially-uniform           *
* workload distribution.   Note that no error checking is done.         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::chopBoxesWithUniformWorkload(
   hier::BoxArray<DIM>& out_boxes,
   tbox::Array<double>& out_workloads,
   const hier::BoxList<DIM>& in_boxes,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int level_number,
   const hier::IntVector<DIM>& min_size,
   const hier::IntVector<DIM>& max_size,
   const hier::IntVector<DIM>& cut_factor,
   const hier::IntVector<DIM>& bad_interval,
   const hier::BoxArray<DIM>& physical_domain) const
{
   NULL_USE(hierarchy); 

   /*
    * Create copy of input box list to prevent changing it.
    */

   hier::BoxList<DIM> tmp_in_boxes_list(in_boxes);

   /*
    * Chop any boxes in input box list that are larger than max box size
    * if possible.
    */

   hier::BoxUtilities<DIM>::chopBoxes(tmp_in_boxes_list,
                                 max_size,
                                 min_size,
                                 cut_factor,
                                 bad_interval,
                                 physical_domain);

   double total_work = 0.0;
   for (typename hier::BoxList<DIM>::Iterator ib0(tmp_in_boxes_list); ib0; ib0++) {
      total_work += ib0().size();
   }

   double work_factor = getMaxWorkloadFactor(level_number);
   double average_work = work_factor * total_work / tbox::MPI::getNodes();

   hier::BoxList<DIM> tmp_box_list;
   tbox::List<double> tmp_work_list;
   BalanceUtilities<DIM>::recursiveBisectionUniform(tmp_box_list,
                                                     tmp_work_list,
                                                     tmp_in_boxes_list,
                                                     average_work,
                                                     min_size,
                                                     cut_factor,
                                                     bad_interval,
                                                     physical_domain);

   if (tmp_box_list.getNumberItems() != tmp_work_list.getNumberItems()) {
      TBOX_ERROR(d_object_name << ": "
        << "Number of boxes generated != number of workload values generated."
        << endl);
   }

   /*
    * Set output box array to list of chopped boxes and set workload array.
    */

   out_boxes = tmp_box_list;

   out_workloads.resizeArray(out_boxes.getNumberOfBoxes());
   int i = 0;
   for (tbox::List<double>::Iterator il(tmp_work_list); il; il++) {
      out_workloads[i] = il();
      i++;
   }

}

/*
*************************************************************************
*                                                                       *
* Private function that chops boxes on a list into another list of      *
* boxes all of which have approximately or less than an average         *
* workload estimate.  The routine assumes a spatially-nonuniform        *
* workload distribution.  Note that no error checking is done.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::chopBoxesWithNonuniformWorkload(
   hier::BoxArray<DIM>& out_boxes,
   tbox::Array<double>& out_workloads,
   const hier::BoxList<DIM>& in_boxes,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   int level_number,
   const hier::IntVector<DIM>& ratio_to_hierarchy_level_zero,
   int wrk_indx,
   const hier::IntVector<DIM>& min_size,
   const hier::IntVector<DIM>& max_size,
   const hier::IntVector<DIM>& cut_factor,
   const hier::IntVector<DIM>& bad_interval,
   const hier::BoxArray<DIM>& physical_domain) const
{

   /*
    * Create copy of input box list to prevent changing it.
    */

   hier::BoxList<DIM> tmp_in_boxes_list(in_boxes);

   hier::BoxUtilities<DIM>::chopBoxes(tmp_in_boxes_list,
                                 max_size,
                                 min_size,
                                 cut_factor,
                                 bad_interval,
                                 physical_domain);

   /*
    * Create temporary patch level from in_boxes, distributed using simple 
    * uniform workload estimate.  Then, fill the patch data on the level 
    * corresponding to the non-uniform workload estimate.  Next, accumulate
    * the total work for the set of boxes.
    */

   hier::BoxArray<DIM> tmp_level_boxes(tmp_in_boxes_list);

   const int num_tmp_patches = tmp_level_boxes.getNumberOfBoxes();
   tbox::Array<double> tmp_level_workloads(num_tmp_patches);
   for (int i = 0; i < num_tmp_patches; i++) {
      tmp_level_workloads[i] = tmp_level_boxes.getBox(i).size();
   }

   hier::ProcessorMapping tmp_level_mapping;

   binPackBoxes(tmp_level_boxes,
                tmp_level_mapping,
                tmp_level_workloads,
                "GREEDY"); 

   tbox::Pointer< hier::PatchLevel<DIM> > tmp_level = 
      new hier::PatchLevel<DIM>(tmp_level_boxes,
                           tmp_level_mapping,
                           ratio_to_hierarchy_level_zero,
                           hierarchy->getGridGeometry(),
                           hierarchy->getPatchDescriptor());

   tmp_level->allocatePatchData(wrk_indx);

   xfer::RefineAlgorithm<DIM> fill_work_algorithm;

   tbox::Pointer< xfer::RefineOperator<DIM> > work_refine_op =
      new pdat::CellDoubleConstantRefine<DIM>();

   fill_work_algorithm.registerRefine(wrk_indx,
                                      wrk_indx,
                                      wrk_indx,
                                      work_refine_op);  

   tbox::Pointer< hier::PatchLevel<DIM> > current_level(NULL);
   if (level_number <= hierarchy->getFinestLevelNumber()) {
      current_level = hierarchy->getPatchLevel(level_number);
   }
 
   fill_work_algorithm.createSchedule(tmp_level,
                                      current_level,
                                      level_number-1,
                                      hierarchy,
                                      (xfer::RefinePatchStrategy<DIM>*)NULL)->fillData(0.0);

   double local_work = 0;
   for (typename hier::PatchLevel<DIM>::Iterator ip(tmp_level); ip; ip++) {
      tbox::Pointer< hier::Patch<DIM> > patch = tmp_level->getPatch(ip());

      double patch_work = 
         BalanceUtilities<DIM>::computeNonUniformWorkload(patch,
                                                           wrk_indx,
                                                           patch->getBox());

      local_work += patch_work;
   }
            
   double total_work = tbox::MPI::sumReduction(local_work);

   double work_factor = getMaxWorkloadFactor(level_number);
   double average_work = work_factor * total_work / tbox::MPI::getNodes();

   hier::BoxList<DIM> tmp_box_list;
   tbox::List<double> tmp_work_list;
   BalanceUtilities<DIM>::recursiveBisectionNonuniform(tmp_box_list,
                                                        tmp_work_list,
                                                        tmp_level,
                                                        wrk_indx,
                                                        average_work,
                                                        min_size,
                                                        cut_factor,
                                                        bad_interval,
                                                        physical_domain);

   tmp_level->deallocatePatchData(wrk_indx);
   tmp_level.setNull();

   if (tmp_box_list.getNumberItems() != tmp_work_list.getNumberItems()) {
      TBOX_ERROR(d_object_name << ": "
        << "Number of boxes generated != number of workload values generated."
        << endl);
   }

   /*
    * Set local box array to list of chopped boxes and set local workload array.
    */

   hier::BoxArray<DIM> local_out_boxes(tmp_box_list);

   tbox::Array<double> local_out_workloads(local_out_boxes.getNumberOfBoxes());

   int i = 0;
   for (tbox::List<double>::Iterator il(tmp_work_list); il; il++) {
      local_out_workloads[i] = il();
      i++;
   }

   /*
    * Gather local box and work arrays so that each processor has a copy.
    */

   hier::BoxComm<DIM>::exchangeBoxArraysAndWeightArrays(local_out_boxes,
                                                   local_out_workloads,
                                                   out_boxes,
                                                   out_workloads);

}

/*
*************************************************************************
*                                                                       *
* Print out all attributes of class instance for debugging.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::printClassData(ostream& os) const
{
   os << "\nLoadBalancer<DIM>::printClassData..." << endl;
   os << "\nLoadBalancer<DIM>: this = "
      << (LoadBalancer<DIM>*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_processor_layout_specified = "
      << d_processor_layout_specified << endl;
   os << "d_processor_layout = "
      << d_processor_layout << endl;
   os << "d_master_workload_data_id = "
      << d_master_workload_data_id << endl;
   os << "d_master_max_workload_factor = "
      << d_master_max_workload_factor << endl;
   os << "d_master_bin_pack_method = "
      << d_master_bin_pack_method << endl;

   int ln;

   os << "d_workload_data_id..." << endl;
   for (ln = 0; ln < d_workload_data_id.getSize(); ln++) {
      os << "    d_workload_data_id[" << ln << "] = "
         << d_workload_data_id[ln] << endl;
   }
   os << "d_max_workload_factor..." << endl;
   for (ln = 0; ln < d_max_workload_factor.getSize(); ln++) {
      os << "    d_max_workload_factor[" << ln << "] = "
         << d_max_workload_factor[ln] << endl;
   }
   os << "d_bin_pack_method..." << endl;
   for (ln = 0; ln < d_bin_pack_method.getSize(); ln++) {
      os << "    d_bin_pack_method[" << ln << "] = "
         << d_bin_pack_method[ln] << endl;
   }
   os << "d_ignore_level_box_union_is_single_box = " << d_ignore_level_box_union_is_single_box << endl;

}

/*
*************************************************************************
*                                                                       *
* Read values (described in the class header) from input database.      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{

   if (!db.isNull()) {

      d_master_bin_pack_method = 
         db->getStringWithDefault("bin_pack_method", 
                                  d_master_bin_pack_method);
      if ( !(d_master_bin_pack_method == "GREEDY" ||
             d_master_bin_pack_method == "SPATIAL") ) {
         TBOX_WARNING(d_object_name << ": "
             << "Unknown 'bin_pack_method' " << d_master_bin_pack_method 
             << " found in input. \nDefault 'GREEDY' will be used." << endl);
         d_master_bin_pack_method = "GREEDY";
      }

      d_master_max_workload_factor = 
         db->getDoubleWithDefault("max_workload_factor", 
                                  d_master_max_workload_factor);
      if (d_master_max_workload_factor <= 0.0) {
         TBOX_WARNING(d_object_name << ": "
             << "Invalid 'max_workload_factor' value " 
             << d_master_max_workload_factor 
             << " found in input. \nDefault 1.0 will be used." << endl);
         d_master_max_workload_factor = 1.0;
      }

      d_ignore_level_box_union_is_single_box =
         db->getBoolWithDefault("ignore_level_box_union_is_single_box",
                                d_ignore_level_box_union_is_single_box);

      d_processor_layout_specified = false;
      int temp_processor_layout[DIM];
      if (db->keyExists("processor_layout")) {
         db->getIntegerArray("processor_layout", temp_processor_layout, DIM);

         /* consistency check */
         int totprocs = 1;
         for (int n = 0; n < DIM; n++) {
            totprocs *= temp_processor_layout[n];
         }
  
         if (totprocs != tbox::MPI::getNodes()) {
            TBOX_WARNING(d_object_name << ": "
               << "Input values for 'processor_layout' are inconsistent with"
               << "\nnumber of processors.  Processor layout information will"
               << "\nbe generated when needed." << endl);
         } else {
            for (int n = 0; n < DIM; n++) {
               d_processor_layout(n) = temp_processor_layout[n];
            }
            d_processor_layout_specified = true;
         }
      }

     d_ignore_level_box_union_is_single_box = 
         db->getBoolWithDefault("ignore_level_box_union_is_single_box", 
                                d_ignore_level_box_union_is_single_box);

   } 

}

/*
*************************************************************************
*                                                                       *
* Private utility function to map boxes to processors.                  *
* Note that no error checking is done.                                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void LoadBalancer<DIM>::binPackBoxes(
   hier::BoxArray<DIM>& boxes,
   hier::ProcessorMapping& mapping,
   tbox::Array<double>& workloads,
   const string& bin_pack_method) const
{
   t_bin_pack_boxes->start();
   /*
    * Sort boxes in order of highest to lowest workload and assign
    * boxes to processors by bin packing.
    */
   t_bin_pack_boxes_sort->start();
   BalanceUtilities<DIM>::sortDescendingBoxWorkloads(boxes,
                                                     workloads);
   t_bin_pack_boxes_sort->stop();

   /*
    * Finally, assign boxes to processors by bin packing.
    */

   int num_procs = tbox::MPI::getNodes();

   t_bin_pack_boxes_pack->start();
   if (bin_pack_method == "SPATIAL") {

     (void) BalanceUtilities<DIM>::spatialBinPack(mapping,
                                                   workloads,
                                                   boxes,
                                                   num_procs);

   } else if (bin_pack_method == "GREEDY") {

     (void) BalanceUtilities<DIM>::binPack(mapping,
                                            workloads,
                                            num_procs);

   } else {

     TBOX_ERROR(d_object_name << ": "
        << "Unknown bin pack method " << bin_pack_method << " found." << endl);

   }
   t_bin_pack_boxes_pack->stop();

   t_bin_pack_boxes->stop();
}

/*
*************************************************************************
*                                                                       *
* Private utility functions to determine parameter values for level.    *
*                                                                       *
*************************************************************************
*/

template<int DIM> int LoadBalancer<DIM>::getWorkloadDataId(int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number >= 0);
#endif
   int wrk_id = (level_number < d_workload_data_id.getSize() ?
                 d_workload_data_id[level_number] :
                 d_master_workload_data_id);

   return(wrk_id);
}

template<int DIM> double LoadBalancer<DIM>::getMaxWorkloadFactor(int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number >= 0);
#endif
   double factor = (level_number < d_max_workload_factor.getSize() ?
                    d_max_workload_factor[level_number] :
                    d_master_max_workload_factor);

   return(factor);
}

template<int DIM> string LoadBalancer<DIM>::getBinPackMethod(int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number >= 0);
#endif
   string factor = (level_number < d_bin_pack_method.getSize() ?
                    d_bin_pack_method[level_number] :
                    d_master_bin_pack_method);

   return(factor);
}

}
}
#endif
