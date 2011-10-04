/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer_C
#define included_mesh_TreeLoadBalancer_C

#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxContainerOrderedIterator.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxSet.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cmath>

#ifndef SAMRAI_INLINE
#include "SAMRAI/mesh/TreeLoadBalancer.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int TreeLoadBalancer::TreeLoadBalancer_LOADTAG0 = 12120001;
const int TreeLoadBalancer::TreeLoadBalancer_LOADTAG1 = 12120002;
const int TreeLoadBalancer::TreeLoadBalancer_EDGETAG0 = 12120003;
const int TreeLoadBalancer::TreeLoadBalancer_EDGETAG1 = 12120004;
const int TreeLoadBalancer::TreeLoadBalancer_PREBALANCE0 = 12120005;
const int TreeLoadBalancer::TreeLoadBalancer_PREBALANCE1 = 12120006;
const int TreeLoadBalancer::TreeLoadBalancer_FIRSTDATALEN = 1000;

const int TreeLoadBalancer::d_default_data_id = -1;

/*
 *************************************************************************
 * TreeLoadBalancer constructor.
 *************************************************************************
 */

TreeLoadBalancer::TreeLoadBalancer(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::Pointer<tbox::Database> input_db):
   d_dim(dim),
   d_mpi_dup(tbox::SAMRAI_MPI::commNull),
   d_n_root_cycles(-1),
   d_degree(2),
   d_master_workload_data_id(d_default_data_id),
   d_balance_penalty_wt(1.0),
   d_surface_penalty_wt(1.0),
   d_slender_penalty_wt(1.0),
   d_slender_penalty_threshold(3.0),
   d_precut_penalty_wt(1.0),
   // Data shared during balancing.
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_min_size(d_dim),
   d_max_size(d_dim),
   d_domain_boxes(d_dim),
   d_bad_interval(d_dim),
   d_cut_factor(d_dim),
   // Output control.
   d_report_load_balance(false),
   // Performance evaluation.
   d_barrier_before(false),
   d_barrier_after(false),
   d_print_steps(false),
   d_print_break_steps(false),
   d_print_swap_steps(false),
   d_print_edge_steps(false),
   d_check_connectivity(false),
   d_check_map(false)
{
   TBOX_ASSERT(!name.empty());

   d_object_name = name;
   d_using_all_procs = true;

   getFromInput(input_db);

   setTimers();

}



/*
 *************************************************************************
 * TreeLoadBalancer constructor.
 *************************************************************************
 */

TreeLoadBalancer::~TreeLoadBalancer()
{
   freeMPICommunicator();
}



/*
 *************************************************************************
 * Accessory functions to get/set load balancing parameters.
 *************************************************************************
 */

bool TreeLoadBalancer::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}



/*
**************************************************************************
**************************************************************************
*/
void TreeLoadBalancer::setWorkloadPatchDataIndex(
   int data_id,
   int level_number)
{
   tbox::Pointer<pdat::CellDataFactory<double> > datafact =
      hier::VariableDatabase::getDatabase()->getPatchDescriptor()->
      getPatchDataFactory(data_id);
   if (datafact.isNull()) {
      TBOX_ERROR(
         d_object_name << " error: "
                       << "\n   data_id " << data_id << " passed to "
                       << "setWorkloadPatchDataIndex()"
                       << " does not refer to cell-centered double patch data. " << std::endl);
   }

   if (level_number >= 0) {
      int asize = d_workload_data_id.getSize();
      if (asize < level_number + 1) {
         d_workload_data_id.resizeArray(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
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



/*
**************************************************************************
**************************************************************************
*/
void TreeLoadBalancer::setUniformWorkload(
   int level_number)
{
   d_workload_data_id[level_number] = -1;
}



/*
 *************************************************************************
 *************************************************************************
 */
void TreeLoadBalancer::loadBalanceBoxLevel(
   hier::BoxLevel& balance_mapped_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const hier::Connector& balance_to_attractor,
   const hier::Connector& attractor_to_balance,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_mapped_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   NULL_USE(balance_to_attractor);
   NULL_USE(attractor_to_balance);
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   TBOX_ASSERT(anchor_to_balance.isInitialized() ==
      balance_to_anchor.isInitialized());
   if (anchor_to_balance.isInitialized()) {
      TBOX_ASSERT(anchor_to_balance.isTransposeOf(balance_to_anchor));
   }
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS6(d_dim,
      balance_mapped_box_level,
      min_size,
      max_size,
      domain_mapped_box_level,
      bad_interval,
      cut_factor);
   if (!hierarchy.isNull()) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);
   }

   d_mpi =
      d_mpi_dup.getCommunicator() == tbox::SAMRAI_MPI::commNull ?
      balance_mapped_box_level.getMPI() : d_mpi_dup;

   if (d_print_steps ||
       d_print_break_steps) {
      tbox::plog << "TreeLoadBalancer::loadBalanceBoxLevel called with:"
                 << "\n  min_size = " << min_size
                 << "\n  max_size = " << max_size
                 << "\n  bad_interval = " << bad_interval
                 << "\n  cut_factor = " << cut_factor
                 << std::endl << balance_mapped_box_level.format("", 2);
   }

   {
      /*
       * Periodic image Box should be ignored during load
       * balancing.  The load-balanced results should contain no
       * periodic images.
       *
       * To avoid need for special logic to skip periodic
       * images while load balancing, we just remove
       * periodic images in the balance_mapped_box_level and
       * all periodic edges in anchor<==>balance.
       */

      balance_mapped_box_level.removePeriodicImageBoxes();
      if (balance_to_anchor.isInitialized()) {

         anchor_to_balance.removePeriodicRelationships();
         anchor_to_balance.initialize(
            anchor_to_balance.getBase(),
            balance_mapped_box_level,
            anchor_to_balance.getConnectorWidth(),
            anchor_to_balance.getNeighborhoodSets());

         balance_to_anchor.removePeriodicRelationships();
         balance_to_anchor.initialize(
            balance_mapped_box_level,
            balance_to_anchor.getHead(),
            balance_to_anchor.getConnectorWidth(),
            balance_to_anchor.getNeighborhoodSets());

      }
   }

   const hier::OverlapConnectorAlgorithm oca;

   if (d_barrier_before) {
      t_barrier_before->start();
      d_mpi.Barrier();
      t_barrier_before->stop();
   }

   if (!rank_group.containsAllRanks()) {
      prebalanceBoxLevel(balance_mapped_box_level,
         balance_to_anchor,
         anchor_to_balance,
         rank_group);
   }

   t_load_balance_mapped_box_level->start();

   /*
    * Shadow data is internal duplicates of some parameters
    * so they do not have to appear in all the interfaces.
    */
   setShadowData(min_size,
      max_size,
      domain_mapped_box_level,
      bad_interval,
      cut_factor,
      balance_mapped_box_level.getRefinementRatio());

   if (d_print_steps) {
      tbox::plog << "Pre balanced:\n" << balance_mapped_box_level.format("", 2);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (balance_to_attractor.isInitialized()) {
      /*
       * If balance_to_attractor is given, sanity-check it.
       */
      if (&balance_mapped_box_level != &balance_to_attractor.getBase() &&
          !(balance_mapped_box_level == balance_to_attractor.getBase())) {
         TBOX_ERROR(
            "TreeLoadBalancer::loadBalanceBoxLevel: balance_mapped_box_level\n"
            << "does not match the base of balance_to_attractor.");
      }
      if (!balance_to_attractor.isTransposeOf(attractor_to_balance)) {
         TBOX_ERROR("TreeLoadBalancer::loadBalanceBoxLevel:\n"
            << "attractor_to_balance and balance_to_attractor\n"
            << "are not transposes of each other.");
      }
   }
#endif

   double local_load, global_sum_load;

   t_compute_local_load->start();
   local_load = computeLocalLoads(balance_mapped_box_level);
   t_compute_local_load->stop();

   size_t nproc_with_initial_load =
      balance_mapped_box_level.getLocalNumberOfBoxes() == 0 ? 0 : 1;

   {
      /*
       * Global reduction for total load and number of procs that has
       * any initial load.
       *
       * TODO: If there's more than one rank group, shouldn't this a
       * global reduction for each rank group instead of a single one
       * for all?
       */
      t_compute_global_load->start();
      if (d_mpi.getSize() > 1) {
         double dtmp[2], dtmp_sum[2];
         dtmp[0] = local_load;
         dtmp[1] = (double)nproc_with_initial_load;
         d_mpi.Allreduce(dtmp,
            dtmp_sum,
            2,
            MPI_DOUBLE,
            MPI_SUM);
         global_sum_load = dtmp_sum[0];
         nproc_with_initial_load = (size_t)dtmp_sum[1];
      } else {
         global_sum_load = local_load;
      }
      t_compute_global_load->stop();
      if (d_print_steps) {
         tbox::plog << "TreeLoadBalancer::loadBalanceBoxLevel balancing "
                    << global_sum_load << " (initially born on "
                    << nproc_with_initial_load << " procs) across all "
                    << d_mpi.getSize()
                    << " procs, averaging " << global_sum_load / d_mpi.getSize()
                    << " or " << pow(global_sum_load / d_mpi.getSize(), 1.0 / d_dim.getValue())
                    << "^" << d_dim << " per proc." << std::endl;
      }
   }

   int output_nproc = rank_group.size();
   d_global_avg_load = global_sum_load / output_nproc;

   int number_of_cycles = d_n_root_cycles;
   if (number_of_cycles < 0) {
      /*
       * User requested automatic number of cycles.
       *
       * Heuristic algorithm: Use 1 cycle unless too few (less than
       * sqrt(nproc)) processors own the initial load.  Exception: use
       * 1 cycle if nproc is small.
       */
      if (balance_mapped_box_level.getMPI().getSize() < 64 ||
          int(nproc_with_initial_load * nproc_with_initial_load) >
          balance_mapped_box_level.getMPI().getSize()) {
         number_of_cycles = 1;
      } else {
         number_of_cycles = 2;
      }
   }

   // Temporary object.
   hier::BoxLevel tmp_mapped_box_level(d_dim);

   for (int n = 0; n < number_of_cycles; ++n) {

      if (n > 0) {
         t_compute_local_load->start();
         local_load = computeLocalLoads(balance_mapped_box_level);
         t_compute_local_load->stop();
      }

      if (d_report_load_balance) {
         // Debugging: check overall load balance at intermediate cycles.
         tbox::plog
         << "TreeLoadBalancer::loadBalanceBoxLevel results before cycle "
         << n << ":" << std::endl;
         TreeLoadBalancer::gatherAndReportLoadBalance(
            local_load,
            balance_mapped_box_level.getMPI());
      }

      bool using_this_rank = rank_group.isMember(d_mpi.getRank());

      /*
       * 1.  Rebalance balance_mapped_box_level and store in temporary tmp_mapped_box_level.
       * 2.  Use the balancing map to modify connectors to anchor.
       * 3.  Reset balance_mapped_box_level to tmp_mapped_box_level.
       */
      hier::Connector balance_to_tmp, tmp_to_balance;
      t_get_map->start();
      if (using_this_rank) {
         loadBalanceBoxLevel_rootCycle(n,
            number_of_cycles,
            local_load,
            global_sum_load,
            balance_mapped_box_level,
            rank_group,
            tmp_mapped_box_level,
            balance_to_tmp,
            tmp_to_balance);
      } else {
         //initialize empty last 3 args
         tmp_mapped_box_level.initialize(
            balance_mapped_box_level.getRefinementRatio(),
            balance_mapped_box_level.getGridGeometry(),
            balance_mapped_box_level.getMPI());
         balance_to_tmp.initialize(
            balance_mapped_box_level,
            tmp_mapped_box_level,
            hier::IntVector::getZero(d_dim));
         tmp_to_balance.initialize(
            tmp_mapped_box_level,
            balance_mapped_box_level,
            hier::IntVector::getZero(d_dim));
      }
      if (d_barrier_after) {
         t_barrier_after->start();
         d_mpi.Barrier();
         t_barrier_after->stop();
      }
      t_get_map->stop();

      t_use_map->start();
      if (anchor_to_balance.isInitialized()) {
         if (0) {
            tbox::plog << "loadBalanceBoxLevel: unbalanced:\n"
                       << balance_mapped_box_level.format("--> ", 3)
                       << "loadBalanceBoxLevel: balanced:\n" << tmp_mapped_box_level.format(
               "--> ",
               3)
            << "loadBalanceBoxLevel: unbalanced_to_balanced:\n"
            << balance_to_tmp.format("--> ", 3)
            << "loadBalanceBoxLevel: balanced_to_unbalanced:\n"
            << tmp_to_balance.format("--> ", 3);
         }
         const hier::MappingConnectorAlgorithm mca;
         mca.modify(anchor_to_balance,
            balance_to_anchor,
            balance_to_tmp,
            tmp_to_balance,
            &balance_mapped_box_level,
            &tmp_mapped_box_level);
         if (0) {
            tbox::plog << "loadBalanceBoxLevel: anchor:\n"
                       << anchor_to_balance.getBase().format("--> ", 3)
                       << "loadBalanceBoxLevel: balanced:\n" << tmp_mapped_box_level.format(
               "--> ",
               3)
            << "loadBalanceBoxLevel: balance_to_anchor:\n"
            << balance_to_anchor.format("--> ", 3)
            << "loadBalanceBoxLevel: anchor_to_balance:\n"
            << anchor_to_balance.format("--> ", 3);
         }
      } else {
         BoxLevel::swap(balance_mapped_box_level, tmp_mapped_box_level);
      }
      if (d_barrier_after) {
         t_barrier_after->start();
         d_mpi.Barrier();
         t_barrier_after->stop();
      }
      t_use_map->stop();

   }

   /*
    * If max_size is given (positive), constrain boxes to the given
    * max_size.  If not given, skip the enforcement step to save some
    * communications.
    */
   if (max_size > hier::IntVector::getZero(d_dim)) {
      t_constrain_size->start();
      hier::Connector unconstrained_to_constrained;
      mapOversizedBoxes(
         balance_mapped_box_level,
         tmp_mapped_box_level,
         unconstrained_to_constrained);
      if (anchor_to_balance.isInitialized()) {
         const hier::MappingConnectorAlgorithm mca;
         mca.modify(anchor_to_balance,
            balance_to_anchor,
            unconstrained_to_constrained,
            &balance_mapped_box_level,
            &tmp_mapped_box_level);
      } else {
         BoxLevel::swap(balance_mapped_box_level, tmp_mapped_box_level);
      }
      t_constrain_size->stop();
      if (d_print_steps) {
         tbox::plog << " TreeLoadBalancer completed applying unconstrained_to_constrained map"
                    << "\n";
      }
   }

   unsetShadowData();

   t_load_balance_mapped_box_level->stop();

   local_load = computeLocalLoads(balance_mapped_box_level);
   d_load_stat.push_back(local_load);
   d_box_count_stat.push_back(
      static_cast<int>(balance_mapped_box_level.getBoxes().size()));

   if (d_report_load_balance) {
      t_report_loads->start();
      tbox::plog
      << "TreeLoadBalancer::loadBalanceBoxLevel results after "
      << number_of_cycles << " cycles:" << std::endl;
      TreeLoadBalancer::gatherAndReportLoadBalance(local_load,
         balance_mapped_box_level.getMPI());
      t_report_loads->stop();
   }

   if (d_check_connectivity && anchor_to_balance.isInitialized()) {
      tbox::plog << "TreeLoadBalancer checking balance-anchor connectivity."
                 << std::endl;
      int errs = 0;
      if (oca.checkOverlapCorrectness(anchor_to_balance, false, true, true)) {
         ++errs;
         tbox::perr << "Error found in anchor_to_balance!\n";
      }
      if (oca.checkOverlapCorrectness(balance_to_anchor, false, true, true)) {
         ++errs;
         tbox::perr << "Error found in balance_to_anchor!\n";
      }
      if (anchor_to_balance.checkTransposeCorrectness(
             balance_to_anchor)) {
         ++errs;
         tbox::perr << "Error found in balance-anchor transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "anchor_mapped_box_level:\n" << anchor_to_balance.getBase().format("", 2)
            << "balance_mapped_box_level:\n" << balance_mapped_box_level.format("", 2)
            << "anchor_to_balance:\n" << anchor_to_balance.format("", 2)
            << "balance_to_anchor:\n" << balance_to_anchor.format("", 2));
      }
      tbox::plog << "TreeLoadBalancer checked balance-anchor connectivity."
                 << std::endl;
   }

   if (d_barrier_after) {
      t_barrier_after->start();
      d_mpi.Barrier();
      t_barrier_after->stop();
   }

   d_mpi.setCommunicator(tbox::SAMRAI_MPI::commNull);
}



/*
 *************************************************************************
 *************************************************************************
 */
void TreeLoadBalancer::mapOversizedBoxes(
   const hier::BoxLevel& unconstrained,
   hier::BoxLevel& constrained,
   hier::Connector& unconstrained_to_constrained) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, unconstrained, constrained);

   t_map_big_boxes->start();

   if (d_print_break_steps) {
      tbox::plog << "Mapping oversized boxes starting with "
                 << unconstrained.getBoxes().size() << " boxes."
                 << std::endl;
   }

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   constrained.initialize(unconstrained.getRefinementRatio(),
      unconstrained.getGridGeometry(),
      unconstrained.getMPI());
   hier::NeighborhoodSet unconstrained_eto_constrained(d_dim);

   const hier::BoxSet& unconstrained_mapped_boxes = unconstrained.getBoxes();

   hier::LocalId next_available_index = unconstrained.getLastLocalId() + 1;

   for (hier::BoxSet::OrderedConstIterator ni = unconstrained_mapped_boxes.orderedBegin();
        ni != unconstrained_mapped_boxes.orderedEnd(); ++ni) {

      const hier::Box& mapped_box = *ni;

      const hier::IntVector mapped_box_size = mapped_box.numberCells();

      if (mapped_box_size <= d_max_size) {
         if (d_print_break_steps) {
            tbox::plog << "    Not oversized: " << mapped_box
                       << mapped_box.numberCells() << "\n";
         }
         constrained.addBox(mapped_box);
      } else {

         if (d_print_break_steps) {
            tbox::plog << "    Breaking oversized " << mapped_box
                       << mapped_box.numberCells() << " ->";
         }
         hier::BoxList chopped(mapped_box);
         hier::BoxUtilities::chopBoxes(chopped,
            d_max_size,
            d_min_size,
            d_cut_factor,
            d_bad_interval,
            d_domain_boxes);

         if (true || chopped.size() != 1) {

            hier::BoxSet& constrained_nabrs =
               unconstrained_eto_constrained.getNeighborSet(mapped_box.getId(),
                                                            d_dim);

            for (hier::BoxList::Iterator li(chopped);
                 li != chopped.end(); ++li) {

               const hier::Box new_box = *li;

               const hier::Box new_mapped_box(new_box,
                                              next_available_index++,
                                              d_mpi.getRank(),
                                              (*ni).getBlockId());

               if (d_print_break_steps) {
                  tbox::plog << "  " << new_mapped_box
                             << new_mapped_box.numberCells();
               }

               constrained.addBox(new_mapped_box);

               constrained_nabrs.insert(constrained_nabrs.orderedEnd(),
                  new_mapped_box);

            }

            if (d_print_break_steps) {
               tbox::plog << "\n";
            }

         } else {
            if (d_print_break_steps) {
               tbox::plog << " Unbreakable!" << "\n";
            }
            constrained.addBox(mapped_box);
         }

      }

   }

   unconstrained_to_constrained.swapInitialize(
      unconstrained,
      constrained,
      zero_vector,
      unconstrained_eto_constrained);
   unconstrained_to_constrained.setConnectorType(hier::Connector::MAPPING);

   if (d_print_steps) {
      tbox::plog
      << " TreeLoadBalancer::mapOversizedBoxes completed building unconstrained_to_constrained"
      << "\n";
   }
   t_map_big_boxes->stop();
}



/*
 *************************************************************************
 * Balance a BoxLevel by using a single cycle of the tree load
 * balancer algorithm.
 *************************************************************************
 */
void TreeLoadBalancer::loadBalanceBoxLevel_rootCycle(
   const int cycle_number,
   const int number_of_cycles,
   const double local_load,
   const double global_sum_load,
   const hier::BoxLevel& unbalanced_mapped_box_level,
   const tbox::RankGroup& rank_group,
   hier::BoxLevel& balanced_mapped_box_level,
   hier::Connector& unbalanced_to_balanced,
   hier::Connector& balanced_to_unbalanced) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim,
      unbalanced_mapped_box_level,
      balanced_mapped_box_level);

   const tbox::SAMRAI_MPI &mpi(unbalanced_mapped_box_level.getMPI());

   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* child_comms = NULL;
   tbox::AsyncCommPeer<int>* parent_send = NULL;
   tbox::AsyncCommPeer<int>* parent_recv = NULL;
   int num_children = 0;

   int g_count, g_id, g_rank, g_nproc;
   createBalanceSubgroups(cycle_number,
      number_of_cycles,
      rank_group,
      g_count,
      g_id,
      g_rank,
      g_nproc,
      num_children,
      child_comms,
      parent_send,
      parent_recv,
      comm_stage);

   double group_sum_load, group_avg_load;

   comm_stage.setCommunicationWaitTimer(t_MPI_wait);


   /*
    * Compute the load for the group.  If this is the last cycle,
    * the group must include all processes, and the group's load
    * is the global sum load.  Else, use all-reduce to get the
    * group load.
    */
   t_compute_tree_load->start();
   if (cycle_number == number_of_cycles - 1) {
      group_sum_load = global_sum_load;
   } else {
      switch (cycle_number) {
         case 0: t_compute_tree_load0->start();
            break;
         case 1: t_compute_tree_load1->start();
            break;
         case 2: t_compute_tree_load2->start();
            break;
      }
      /*
       * Use MPI's vector all-reduce to get individual group loads.
       * This gives more info than the process needs, but because the
       * number of groups << number of procs, it is still faster
       * (probably) than hand coded conmunication.
       */
      std::vector<double> group_loads(g_count, 0.0);
      group_loads[g_id] = local_load;
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&group_loads[0],
            static_cast<int>(group_loads.size()),
            MPI_SUM);
      }
      group_sum_load = group_loads[g_id];
      switch (cycle_number) {
         case 0: t_compute_tree_load0->stop();
            break;
         case 1: t_compute_tree_load1->stop();
            break;
         case 2: t_compute_tree_load2->stop();
            break;
      }
   }
   t_compute_tree_load->stop();

   group_avg_load = group_sum_load / g_nproc;

   if (d_print_steps) {
      tbox::plog << "Initial local,group_sum,group_avg loads:"
                 << "  " << local_load
                 << "  " << group_sum_load
                 << "  " << group_avg_load
                 << std::endl;
   }


   /*
    * Before the last cycle, it is possible for the group average load
    * to be below the global average, if the group just happens to have
    * underloaded processors.  However, there is no point in driving the
    * processor loads down below the global average just to have it
    * brought back up by the last cycle.  It just unnecessarily fragments
    * the boxes and costs more to do.  To prevent this, reset the group
    * average to the global average if it is below.
    */
   group_avg_load =
      tbox::MathUtilities<double>::Max(group_avg_load, d_global_avg_load);

   // List of all ranks we communicate with.
   std::vector<int> peer_ranks(num_children + 1);
   for (int i = 0; i < num_children; ++i) {
      peer_ranks[i] = child_comms[i].getPeerRank();
   }
   peer_ranks[num_children] = parent_send ==
      NULL ? -1 : parent_send->getPeerRank();

   const hier::BoxSet& unbalanced_mapped_boxes =
      unbalanced_mapped_box_level.getBoxes();

   /*
    * Outline for load balancing:
    *
    * 1. For each child:
    * Receive data from subtree rooted at child (number in
    * subtree, excess work, remaining work in subtree, etc.).
    *
    * 2. Compute data for subtree rooted at self by combining
    * local data with children subtree data.
    *
    * 3. If parent exists:
    * Compute subtree info (number in subtree, excess work,
    * remaining work in subtree, etc.) and send to parent.
    *
    * 4. If not root:
    * Receive additional work (if any) from parent.
    *
    * 5. Partition additional work among children and self.
    *
    * 6. For each child:
    * Send additional work (if any).
    *
    * Outline for connecting pre-balance and post-balance Boxes:
    */

   tbox::AsyncCommStage::MemberVec completed;

   /*
    * Post receive for data from subtree rooted at children.
    * We have to do a few local setups, but post the receive
    * now to overlap communication.
    */
   t_get_load_from_children->start();
   for (int c = 0; c < num_children; ++c) {
      child_comms[c].beginRecv();
      if (child_comms[c].isDone()) {
         completed.insert(completed.end(), &child_comms[c]);
      }
   }
   t_get_load_from_children->stop();

   t_misc1->start();

   /*
    * Data for storing and transfering subtree info.
    */
   SubtreeLoadData* child_load_data = new SubtreeLoadData[num_children];
   SubtreeLoadData my_load_data;

   t_misc1->stop();
   t_misc2->start();

   /*
    * For tracking available indices that can be assigned to newly
    * created Boxes: The indices for Boxes from this proc
    * and those from children procs are interlaced to remove results'
    * dependence on which message arrives first.
    *
    * Nodes from the first child start with the first unused index
    * divisible by (2+d_degree).  Nodes from child c start with the
    * same number plus c.  Locally created Boxes start with the
    * same number plus d_degree.  Nodes from parent start with the
    * same number plus d_degree+1.
    *
    * Each time an index is used, next_available_index is increased by 2+degree.
    */
   std::vector<hier::LocalId> next_available_index(2 + d_degree);
   next_available_index[0] = unbalanced_mapped_box_level.getLastLocalId() + 1;
   next_available_index[0] += (next_available_index[0] % (1 + d_degree));
   for (int c = 1; c < d_degree + 2; ++c) {
      next_available_index[c] = next_available_index[0] + c;
   }

   /*
    * Initialize output objects.
    */
   balanced_mapped_box_level.initialize(
      unbalanced_mapped_box_level.getRefinementRatio(),
      unbalanced_mapped_box_level.getGridGeometry(),
      unbalanced_mapped_box_level.getMPI());
   hier::NeighborhoodSet balanced_eto_unbalanced(d_dim);
   hier::NeighborhoodSet unbalanced_eto_balanced(d_dim);
   balanced_to_unbalanced.initialize(
      balanced_mapped_box_level,
      unbalanced_mapped_box_level,
      hier::IntVector::getZero(d_dim));
   unbalanced_to_balanced.initialize(
      unbalanced_mapped_box_level,
      balanced_mapped_box_level,
      hier::IntVector::getZero(d_dim));

   t_misc2->stop();
   t_misc3->start();

   /*
    * Unassigned NodeInTransit.  First set to local excess Box
    * (if any).  If children send up excess NodeInTransit, add them.  Result
    * should be taken by local proc and sent down the tree to children
    * (if needed).  The rest (if any) are sent up to parent.
    */
   TransitSet unassigned;

   /*
    * Compute local proc's Boxes and loads and store in
    * my_load_data.  This will eventually include data for the subtree.
    * We will add the rest of the subtree's work when we receive that
    * data from the children.
    */
   my_load_data.num_procs = 1;
   my_load_data.total_work = (int)computeLocalLoads(unbalanced_mapped_box_level);

   t_misc3->stop();

   /*
    * Decide whether local parts of unbalanced_mapped_box_level should
    * be reassigned.  Populate balanced_mapped_box_level from local loads.
    */

   t_local_balancing->start();

   if (my_load_data.total_work /* currently excluding children */ <=
       group_avg_load) {
      /*
       * Local process is underloaded, so put all of unbalanced_mapped_box_level into
       * the balanced_mapped_box_level (and add more later).
       */
      for (hier::BoxSet::OrderedConstIterator ni = unbalanced_mapped_boxes.orderedBegin();
           ni != unbalanced_mapped_boxes.orderedEnd(); ++ni) {
         balanced_mapped_box_level.addBox(*ni);
      }
   } else {
      /*
       * Local process is overloaded, so reassign some loads:
       * - sort NodeInTransit by load
       * - reassignLoads (put excess loads in unassigned container) and
       * - put remainder in balanced_mapped_box_level.
       *
       * Note: This algorithm would also work if we put all local
       * Boxes into unassigned (instead of just the excess load).
       * In the end,
       * all remaining unassigned Boxes get assigned to the local
       * process anyway.  In fact, having more Boxes in unassigned
       * may let reassignLoads do a better job in reassigning loads
       * to children and parents, because it would have more choices.
       * The reason we place only the excess load into unassigned
       * is to help preserve locality, assuming that current local
       * Boxes may have more local neighbors.  However, practical
       * evidence so far suggest that locality is not that
       * important to performance, at least for explicit solvers.
       * Transfering all local Boxes
       * into unassigned at the start may affect the performance
       * of this algorithms though, because it bypasses one
       * reassignLoads call but makes the unassigned Box set
       * bigger, which may make other calls to reassignLoads
       * work harder.  Which one is a bigger effect and whether
       * there are any significant effects at all remains to
       * be seen.
       */

      TransitSet
      sorted_loads(unbalanced_mapped_boxes.orderedBegin(), unbalanced_mapped_boxes.orderedEnd());
      int ideal_transfer = int(0.5 + my_load_data.total_work - group_avg_load);
      if (d_print_steps) {
         tbox::plog << "  Reassigning initial overload of " << ideal_transfer
                    << " to unassigned.\n";
      }
      int actual_transfer;
      reassignLoads(sorted_loads,
         unassigned,
         ideal_transfer,
         actual_transfer,
         next_available_index[d_degree]);

      for (TransitSet::const_iterator
           ni = sorted_loads.begin(); ni != sorted_loads.end(); ++ni) {
         const BoxInTransit& mapped_box_in_transit = *ni;
         balanced_mapped_box_level.addBox(mapped_box_in_transit.mapped_box);
         /*
          * Create local edges only for mapped_box_in_transit that changed.
          * Semilocal edges are created by final owner
          * and sent back.
          */
         if (mapped_box_in_transit.mapped_box.getLocalId() !=
             mapped_box_in_transit.orig_mapped_box.getLocalId()) {
            balanced_eto_unbalanced.insertNeighbor(
               mapped_box_in_transit.mapped_box.getId(),
               mapped_box_in_transit.orig_mapped_box);
            unbalanced_eto_balanced.insertNeighbor(
               mapped_box_in_transit.orig_mapped_box.getId(),
               mapped_box_in_transit.mapped_box);
         }
      }

      if (d_print_steps) {
         tbox::plog << "    Unassigning " << unassigned.size()
                    << " boxes (" << actual_transfer << " / "
                    << ideal_transfer << " units) "
                    << "for " << (g_nproc - 1) << " procs:";
         for (TransitSet::const_iterator
              ni = unassigned.begin(); ni != unassigned.end(); ++ni) {
            tbox::plog << "  " << *ni;
         }
         tbox::plog << std::endl;
      }
   }

   t_local_balancing->stop();



   /*
    * Complete load communications with children.
    * Add imported BoxInTransit to unassigned bin.
    */
   t_get_load_from_children->start();
   do {
      for (unsigned int i = 0; i < completed.size(); ++i) {
         tbox::AsyncCommPeer<int>* peer_comm =
            dynamic_cast<tbox::AsyncCommPeer<int> *>(completed[i]);
         TBOX_ASSERT(peer_comm >= child_comms);
         TBOX_ASSERT(peer_comm <= child_comms + num_children);
         const int* received_data = peer_comm->getRecvData();
         int received_data_length = peer_comm->getRecvSize();
         const int cindex = static_cast<int>(peer_comm - child_comms);

         TBOX_ASSERT(cindex >= 0 && cindex < num_children);

         /*
          * Extract data from the children.  Initialize data for tree
          * rooted at children.  If child sent up any excess
          * Box, put them in unassigned.
          */
         int old_size = static_cast<int>(unassigned.size());
         unpackSubtreeLoadData(received_data,
            received_data_length,
            child_load_data[cindex],
            unassigned,
            next_available_index[cindex]);
         child_load_data[cindex].ideal_work =
            int(group_avg_load * child_load_data[cindex].num_procs + 0.5);
         if (d_print_steps) {
            TransitSet::const_iterator
               ni = unassigned.begin();
            for (int ii = 0; ii < old_size; ++ii) {
               ++ni;
            }
            if (d_print_steps) {
               tbox::plog << "    Got " << unassigned.size() - old_size
                          << " boxes (" << child_load_data[cindex].load_imported
                          << " units) from child "
                          << peer_comm->getPeerRank() << ":";
               for ( ; ni != unassigned.end(); ++ni) {
                  const BoxInTransit& mapped_box_in_transit = *ni;
                  tbox::plog << "  " << mapped_box_in_transit;
               }
               tbox::plog << std::endl;
            }
         }
         // Sum children load into my_load_data.
         my_load_data.num_procs += child_load_data[cindex].num_procs;
         my_load_data.total_work +=
            child_load_data[cindex].total_work
            + child_load_data[cindex].load_imported;
      }
      completed.clear();
      t_children_load_comm->start();
      comm_stage.advanceSome(completed);
      t_children_load_comm->stop();
   } while (!completed.empty());



   // We should have received everything at this point.
   TBOX_ASSERT(!comm_stage.hasPendingRequests());

   my_load_data.ideal_work = int(group_avg_load * my_load_data.num_procs + 0.5);

   if (d_print_steps) {
      tbox::plog << "Received children subtree data." << std::endl;
      for (int c = 0; c < num_children; ++c) {
         tbox::plog << "Child " << child_comms[c].getPeerRank()
                    << " subtree data: " << child_load_data[c].total_work
                    << "/" << child_load_data[c].ideal_work
                    << " for " << child_load_data[c].num_procs << " procs averaging "
                    << child_load_data[c].total_work / child_load_data[c].num_procs
                    << " after sending up " << child_load_data[c].load_imported
                    << std::endl;
      }
      tbox::plog << "Initial subtree data: " << my_load_data.total_work
                 << " / " << my_load_data.ideal_work
                 << " for " << my_load_data.num_procs << " procs averaging "
                 << my_load_data.total_work / my_load_data.num_procs
                 << " before sending up anything."
                 << std::endl;
   }

   t_get_load_from_children->stop();



   /*
    * Send subtree info and excess work (if any) up to parent.
    */
   t_send_load_to_parent->start();
   if (parent_send != NULL) {
      /*
       * Compute the excess work we want to send to parent.
       */
      int ideal_transfer = my_load_data.total_work > my_load_data.ideal_work ?
         my_load_data.total_work - my_load_data.ideal_work : 0;
      if (ideal_transfer > 0) {
         if (d_print_steps) {
            tbox::plog << "  Attempting to reassign " << ideal_transfer
                       << " of unassigned load to parent.\n";
         }
         int actual_transfer;
         reassignLoads(unassigned,
            my_load_data.for_export /* to parent */,
            ideal_transfer,
            actual_transfer,
            next_available_index[d_degree]);
         my_load_data.load_exported = actual_transfer;
         my_load_data.total_work -= actual_transfer;

         if (d_print_steps) {
            tbox::plog << "  Giving " << my_load_data.for_export.size()
                       << " boxes (" << actual_transfer << " / "
                       << ideal_transfer << " units) to parent "
                       << parent_send->getPeerRank() << " for "
                       << (g_nproc - my_load_data.num_procs) << " procs:";
            for (TransitSet::const_iterator
                 ni = my_load_data.for_export.begin();
                 ni != my_load_data.for_export.end(); ++ni) {
               tbox::plog << "  " << *ni;
            }
            tbox::plog << std::endl;
         }
      }
      /*
       * Send local process's load info up to parent.
       */
      std::vector<int> msg;
      packSubtreeLoadData(msg, my_load_data);
      parent_send->beginSend(&msg[0], static_cast<int>(msg.size()));
   }
   t_send_load_to_parent->stop();



   /*
    * Finish the send-up.
    * To preclude sending work in both directions, the parent
    * will *not* send a work message down if we sent work up.
    */
   if (parent_send != NULL && my_load_data.load_exported == 0) {
      t_parent_load_comm->start();
      t_get_load_from_parent->start();

      parent_recv->beginRecv();

      t_get_load_from_parent->stop();
      t_parent_load_comm->stop();
   }


   /*
    * May do some things here that do not depend on message from parents.
    */


   if (parent_recv != NULL && my_load_data.load_exported == 0) {
      t_get_load_from_parent->start();

      t_parent_load_comm->start();
      parent_recv->completeCurrentOperation();
      t_parent_load_comm->stop();

      int old_size = static_cast<int>(unassigned.size());
      SubtreeLoadData parent_load_data;
      unpackSubtreeLoadData(parent_recv->getRecvData(),
         parent_recv->getRecvSize(),
         parent_load_data,
         unassigned,
         next_available_index[1 + d_degree]);
      my_load_data.load_imported = parent_load_data.load_imported;
      my_load_data.total_work += parent_load_data.load_imported;

      if (d_print_steps) {
         TransitSet::const_iterator
            ni = unassigned.begin();
         for (int i = 0; i < old_size; ++i) {
            ++ni;
         }
         if (d_print_steps) {
            tbox::plog << "    Got " << unassigned.size() - old_size
                       << " boxes (" << parent_load_data.load_imported
                       << " units) from parent "
                       << parent_recv->getPeerRank() << ":";
            for ( ; ni != unassigned.end(); ++ni) {
               const BoxInTransit& mapped_box_in_transit = *ni;
               tbox::plog << "  " << mapped_box_in_transit;
            }
            tbox::plog << std::endl;
         }
      }

      t_get_load_from_parent->stop();
   }

   if (d_print_steps) {
      tbox::plog << "  After parent/children, my total work is "
                 << my_load_data.total_work << " / "
                 << my_load_data.ideal_work << std::endl;
   }

   /*
    * Divide unassigned between self and children subtrees.
    *
    * Give as much as possible of the excess work to children trees
    * so they can start their partitioning process.  The remaining
    * excess is to be kept locally.
    */

   t_send_load_to_children->start();

   /*
    * Compute number of processors that would have to share
    * the unassigned loads.  (Processor subtrees that have
    * detached themselves do not count because they will
    * not receive any unassigned loads.)
    */
   int remaining_num_procs = 1;
   for (int c = 0; c < num_children; ++c) {
      if (child_load_data[c].load_imported == 0) {
         remaining_num_procs += child_load_data[c].num_procs;
      }
   }

   for (int c = 0; c < num_children; ++c) {

      SubtreeLoadData& recip_data = child_load_data[c];

      /*
       * To preclude unneeded messages,
       * we do *not* send a work message down if the child
       * sent work up to us.
       */
      if (recip_data.load_imported == 0) {
         int ideal_transfer = recip_data.ideal_work - recip_data.total_work;
         int actual_transfer = 0;
         if (d_print_steps) {
            tbox::plog << "  Attempting to reassign " << ideal_transfer
                       << " of unassigned load to child "
                       << child_comms[c].getPeerRank() << "\n";
         }
         if (ideal_transfer > 0) {
            remaining_num_procs -= recip_data.num_procs;
            reassignLoads(unassigned,
               recip_data.for_export,
               ideal_transfer,
               actual_transfer,
               next_available_index[d_degree]);
            recip_data.load_exported += actual_transfer;
            recip_data.total_work += actual_transfer;
         }

         if (d_print_steps) {
            tbox::plog << "  Giving " << recip_data.for_export.size()
                       << " boxes (" << actual_transfer << " / " << ideal_transfer
                       << " units) to child " << c << ':'
                       << child_comms[c].getPeerRank() << " for "
                       << recip_data.num_procs
                       << " procs:";
            for (TransitSet::const_iterator
                 ni = recip_data.for_export.begin();
                 ni != recip_data.for_export.end(); ++ni) {
               tbox::plog << "  " << *ni;
            }
            tbox::plog << std::endl;
         }
         std::vector<int> msg;
         packSubtreeLoadData(msg, recip_data);
         child_comms[c].beginSend(&msg[0], static_cast<int>(msg.size()));

      }

   }

   t_send_load_to_children->stop();


   /*
    * Populate balanced_mapped_box_level with unassigned Boxes
    * and generate relationships in balanced<==>unbalanced.
    *
    * Unassigned BoxInTransit that differ from their orig_mapped_box
    * are imported or left-overs from breaking.  These Boxes
    * should have relationships in balanced<==>unbalanced.
    *
    * Unassigned BoxInTransit that are identical to their
    * orig_mapped_box are Boxes that were intended for export but
    * never exported.  These Boxes should NOT have relationships in
    * balanced<==>unbalanced.
    */
   t_local_balancing->start();

   for (TransitSet::iterator
        ni = unassigned.begin();
        ni != unassigned.end(); /* incremented in loop */) {
      const BoxInTransit& mapped_box_in_transit = *ni;
      balanced_mapped_box_level.addBox(mapped_box_in_transit.mapped_box);
      if (!mapped_box_in_transit.mapped_box.isIdEqual(mapped_box_in_transit.orig_mapped_box)) {
         // Create relationship for changed Box.
         balanced_eto_unbalanced.insertNeighbor(
            mapped_box_in_transit.mapped_box.getId(),
            mapped_box_in_transit.orig_mapped_box);
         ++ni;
      } else {
         // Remove unchanged mapped_box_in_transit so packedgeMessages
         // will not create edge.
         TransitSet::iterator sni = ni;
         ++ni;
         unassigned.erase(sni);
      }
   }

   balanced_to_unbalanced.swapInitialize(
      balanced_mapped_box_level,
      unbalanced_mapped_box_level,
      hier::IntVector::getZero(d_dim),
      balanced_eto_unbalanced);
   balanced_to_unbalanced.setConnectorType(hier::Connector::MAPPING);
   unbalanced_to_balanced.swapInitialize(
      unbalanced_mapped_box_level,
      balanced_mapped_box_level,
      hier::IntVector::getZero(d_dim),
      unbalanced_eto_balanced);
   unbalanced_to_balanced.setConnectorType(hier::Connector::MAPPING);

   t_local_balancing->stop();

   /*
    * Now that unchanged BoxInTransit have been removed from
    * unassigned, call it "changed", for clarity.  Changed
    * BoxInTransit need to have edges built.
    */
   TransitSet& changed = unassigned;

   /*
    * Finish messages before starting edge info exchange.
    * We have only sends to complete, so it should not take
    * long to advance them all to completion.
    */
   t_finish_comms->start();
   comm_stage.advanceAll(completed);
   t_finish_comms->stop();
   completed.clear();
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < num_children; ++i) {
      TBOX_ASSERT(child_comms[i].isDone());
   }
   if (parent_send != NULL) {
      TBOX_ASSERT(parent_send->isDone());
      TBOX_ASSERT(parent_recv->isDone());
   }
#endif

#if 1
   /*
    * This is not needed if the communications are bug-free.
    * As a debugging measure, switch the MPI tags for edge communication.
    */
   for (int i = 0; i < num_children; ++i) {
      child_comms[i].setMPITag(TreeLoadBalancer_EDGETAG0,
         TreeLoadBalancer_EDGETAG1);
   }
   if (parent_send != NULL) {
      parent_send->setMPITag(TreeLoadBalancer_EDGETAG0,
         TreeLoadBalancer_EDGETAG1);
      parent_recv->setMPITag(TreeLoadBalancer_EDGETAG0,
         TreeLoadBalancer_EDGETAG1);
   }
#endif

   /*
    * Determine edges from unbalanced to balanced BoxLevels
    * by sending balanced Boxes back to the owners of
    * the unbalanced Boxes that originated them.  Use
    * the proc_hist data from the balanced BoxInTransit
    * to send along the tree.
    *
    * Post receives for edge data sent from children.
    * Pack relationship data for locally owned Boxes.
    * Then complete receives  from children and send
    * data up to parents.
    *
    * Four steps for passing relationship data to the owners of
    * originating Boxes:
    * 1. Receive relationships from all children to whom we exported work.
    * 1a. Note relationships to local Boxes on balanced_mapped_box_level.
    * 2. Send relationship data to parent if we imported work from parent.
    * 3. Receive relationship data from parent if we exported work to parent.
    * 4. Send relationships to children from whom we imported work.
    *
    * TODO: Implement alternate edge generation method where terminal
    * owners send edge info directly back to initial owners.  No need
    * for trails.  Initial owners receive from MPI_ANY and stops
    * receiving when it can account for the number of cells it started
    * with.  We should retain both methods for edge generation because
    * we don't know which is faster on any given machine.  BTNG.
    */

   std::vector<int>* messages = new std::vector<int>[num_children + 1];

   t_get_edge_from_children->start();
   for (int c = 0; c < num_children; ++c) {
      if (child_load_data[c].load_exported > 0) {
         child_comms[c].beginRecv();
         if (child_comms[c].isDone()) {
            completed.insert(completed.end(), &child_comms[c]);
         }
      }
   }

   t_local_edges->start();
   packNeighborhoodSetMessages(messages,
      unbalanced_to_balanced,
      changed,
      peer_ranks);
   t_local_edges->stop();

   /*
    * Advance edge messages from children.
    */

   if (completed.empty()) {
      t_children_edge_comm->start();
      comm_stage.advanceSome(completed);
      t_children_edge_comm->stop();
   }
   while (!completed.empty()) {
      for (unsigned int i = 0; i < completed.size(); ++i) {
         tbox::AsyncCommPeer<int>* peer_comm =
            dynamic_cast<tbox::AsyncCommPeer<int> *>(completed[i]);
         const int* received_data = peer_comm->getRecvData();
         int received_data_length = peer_comm->getRecvSize();
         unpackAndRouteNeighborhoodSets(received_data,
            received_data_length,
            messages,
            unbalanced_to_balanced,
            peer_ranks);
      }
      completed.clear();
      t_children_edge_comm->start();
      comm_stage.advanceSome(completed);
      t_children_edge_comm->stop();
   }
   t_get_edge_from_children->stop();

   if (my_load_data.load_imported > 0) {
      /*
       * We have received BoxInTransit from the parent, so the parent
       * is expecting us to send relationship data.
       */
      t_send_edge_to_parent->start();
      t_parent_edge_comm->start();
      parent_send->completeCurrentOperation();
      parent_send->beginSend(&messages[num_children][0],
         static_cast<int>(messages[num_children].size()));
      t_parent_edge_comm->stop();
      t_send_edge_to_parent->stop();
   } else if (my_load_data.load_exported > 0) {
      /*
       * We have sent BoxInTransit to the parent, so the parent would
       * send back relationship data.
       */
      t_get_edge_from_parent->start();
      t_parent_edge_comm->start();
      parent_recv->beginRecv();
      parent_recv->completeCurrentOperation();
      t_parent_edge_comm->stop();
      unpackAndRouteNeighborhoodSets(parent_recv->getRecvData(),
         parent_recv->getRecvSize(),
         messages,
         unbalanced_to_balanced,
         peer_ranks);
      t_get_edge_from_parent->stop();
   }

   t_send_edge_to_children->start();
   for (int c = 0; c < num_children; ++c) {
      if (child_load_data[c].load_imported > 0) {
         child_comms[c].beginSend(&messages[c][0],
            static_cast<int>(messages[c].size()));
      }
   }
   t_send_edge_to_children->stop();

   if (d_check_connectivity) {
      const hier::OverlapConnectorAlgorithm oca;
      tbox::plog
      << "TreeLoadBalancer checking unbalanced-balanced connectivity."
      << std::endl;
      int errs = 0;
      if (oca.checkOverlapCorrectness(unbalanced_to_balanced, true, true)) {
         ++errs;
         tbox::perr << "Error found in unbalanced_to_balanced!\n";
      }
      if (oca.checkOverlapCorrectness(balanced_to_unbalanced, true, true)) {
         ++errs;
         tbox::perr << "Error found in balanced_to_unbalanced!\n";
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         ++errs;
         tbox::perr << "Error found in balanced-unbalanced transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "unbalanced_mapped_box_level:\n" << unbalanced_mapped_box_level.format("", 2)
            << "balanced_mapped_box_level:\n" << balanced_mapped_box_level.format("", 2)
            << "unbalanced_to_balanced:\n" << unbalanced_to_balanced.format("", 2)
            << "balanced_to_unbalanced:\n" << balanced_to_unbalanced.format("", 2));
      }
   }

   if (d_check_map) {
      const hier::MappingConnectorAlgorithm mca;
      if (mca.findMappingErrors(unbalanced_to_balanced) != 0) {
         TBOX_ERROR(
            "TreeLoadBalancer::loadBalanceBoxLevel_rootCycle Mapping errors found in unbalanced_to_balanced!");
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         TBOX_ERROR(
            "TreeLoadBalancer::loadBalanceBoxLevel_rootCycle Transpose errors found!");
      }
   }

   delete[] child_load_data;
   delete[] messages;

   /*
    * Wrap up communications.
    * No follow-up needed because all the receives have been
    * processed--only the sends may still be out.
    *
    * TODO: This step should take no significant time, because the sends
    * should all be well on their way.  However, I'm seeing some
    * aparent scaling issues.  I may want to use MPI_Request_free
    * to just let the sends finish quietly without furthe actions.
    * Not sure if that would help.
    */
   t_finish_comms->start();
   comm_stage.advanceAll(completed);
   t_finish_comms->stop();

   destroyBalanceSubgroups(child_comms, parent_send, parent_recv);
}



/*
 *************************************************************************
 *
 * The tree load balancer's reassignment method, essentially a two-bin
 * load balancer.  Given two sets of BoxInTransit (the bins) and
 * an amount of work to move from one set to the other, this method
 * makes a best effort to effect the work transfer between the two
 * bins.  This can move BoxInTransit between given sets and, if
 * needed, break some NodeInTransit up to move part of the work.
 *
 * Reassign some BoxInTransit from one set into another by moving
 * them between src and dst containers of BoxInTransit objects.
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 * Effectiveness of reassignLoads is key to getting good balance.
 * Any deviation from the ideal can accumulate each generation in
 * the tree.
 *
 *************************************************************************
 */
void TreeLoadBalancer::reassignLoads(
   TransitSet& src,
   TransitSet& dst,
   int ideal_transfer,
   int& actual_transfer,
   hier::LocalId& next_available_index) const
{
   if (d_print_steps) {
      tbox::plog << "  reassignLoads attempting to reassign "
                 << ideal_transfer << " from src to dst."
                 << std::endl;
   }

   actual_transfer = 0;

   if ((ideal_transfer >= 0 && src.empty()) ||
       (ideal_transfer <= 0 && dst.empty())) {
      return;
   }

   t_reassign_loads->start();

   /*
    * The algorithm cycles through a do-loop.  Each time around, we try
    * to swap some BoxInTransit between src and dst until we cannot improve the
    * actual_transfer any further.  Then, we try breaking up some BoxInTransit to
    * improve the results.  If we break some BoxInTransit, we generate some more
    * swapping options that were not there before, so we loop back to
    * try swapping again.
    *
    * If a break phase does not break any Box (and does not generate more
    * swap options), the loop will stop making changes.  We exit the loop
    * at that point (and whenever we reached the ideal transfer).
    */
   do {

      /*
       * Try to balance load through swapping.
       */
      int swap_transfer;
      shiftLoadsBySwapping(
         src,
         dst,
         ideal_transfer - actual_transfer,
         swap_transfer);
      actual_transfer += swap_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            src,
            dst,
            actual_transfer - ideal_transfer);
         tbox::plog << "  Balance penalty after shiftLoadsBySwapping = "
                    << balance_penalty
                    << ", needs " << (ideal_transfer - actual_transfer)
                    << " more with " << src.size() << " source and "
                    << dst.size() << " dst Boxes remaining."
                    << std::endl;
      }

      if (actual_transfer == ideal_transfer) break;

      /*
       * Assuming that we did the best we could, swapping
       * some BoxInTransit without breaking any, we now break up a Box
       * in the overloaded side for partial transfer to the
       * underloaded side.
       */
      int brk_transfer;
      shiftLoadsByBreaking(src,
         dst,
         ideal_transfer - actual_transfer,
         brk_transfer,
         next_available_index);
      actual_transfer += brk_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            src,
            dst,
            actual_transfer - ideal_transfer);
         tbox::plog << "    Balance penalty after shiftLoadsByBreaking = "
                    << balance_penalty
                    << ", needs " << (ideal_transfer - actual_transfer)
                    << " more with " << src.size() << " source and "
                    << dst.size() << " dst Boxes remaining."
                    << std::endl;
      }
      if (brk_transfer == 0) {
         /*
          * If no Box can be broken to improve the actual_transfer,
          * there is nothing further we can do.  (The swap phase, tried
          * before the break phase, also generated no transfer.)
          */
         break;
      }

      /*
       * Now that we have broken up a Box, redo this loop to
       * see if swapping can produce a better result.
       */
   } while (ideal_transfer != actual_transfer);

   t_reassign_loads->stop();
}



/*
 *************************************************************************
 *************************************************************************
 */
void TreeLoadBalancer::packSubtreeLoadData(
   std::vector<int>& msg,
   const SubtreeLoadData& load_data) const
{
   t_pack_load->start();
   msg.insert(msg.end(), load_data.num_procs);
   msg.insert(msg.end(), load_data.total_work);
   msg.insert(msg.end(), load_data.load_exported);
   const TransitSet& for_export = load_data.for_export;
   msg.insert(msg.end(), static_cast<int>(for_export.size()));
   for (TransitSet::const_iterator
        ni = for_export.begin(); ni != for_export.end(); ++ni) {
      const BoxInTransit& mapped_box_in_transit = *ni;
      int i0 = static_cast<int>(msg.size());
      msg.insert(msg.end(), mapped_box_in_transit.commBufferSize(), 0);
      mapped_box_in_transit.putToIntBuffer(&msg[i0]);
   }
   t_pack_load->stop();
}



/*
 *************************************************************************
 *************************************************************************
 */
void TreeLoadBalancer::unpackSubtreeLoadData(
   const int* received_data,
   int received_data_length,
   SubtreeLoadData& load_data,
   TransitSet& receiving_bin,
   hier::LocalId& next_available_index) const
{
   (void)received_data_length;
   t_unpack_load->start();
   load_data.num_procs = *(received_data++);
   load_data.total_work = *(received_data++);
   load_data.load_imported = *(received_data++);
   const int num_transits = *(received_data++);
   for (int i = 0; i < num_transits; ++i) {
      BoxInTransit received_mapped_box_in_transit(d_dim);
      received_mapped_box_in_transit.getFromIntBuffer(received_data);
      BoxInTransit renamed_mapped_box_in_transit(received_mapped_box_in_transit,
                                                 received_mapped_box_in_transit.getBox(),
                                                 d_mpi.getRank(),
                                                 next_available_index);
      next_available_index += 2 + d_degree;
      receiving_bin.insert(renamed_mapped_box_in_transit);
      received_data += received_mapped_box_in_transit.commBufferSize();
   }
   t_unpack_load->stop();
}



/*
 *************************************************************************
 * This method sets up local relationships in unbalanced_to_balanced.
 * For remote relationships, it packs information into messages to be
 * sent out.
 *
 * The input mapped_boxes_in_transit contains local BoxInTransit that
 * will remain on the local process.  This method builds the
 * relationships from the BoxInTransit to its origin (in
 * unbalanced_to_balanced) if the origin is local.
 * If the origin is remote, pack the BoxInTransit in a message to send
 * in the direction of the originating rank.  The process history
 * tells us whether to send the BoxInTransit to the parent or one of
 * the children.
 *************************************************************************
 */
void TreeLoadBalancer::packNeighborhoodSetMessages(
   std::vector<int> messages[],
   hier::Connector& unbalanced_to_balanced,
   const TransitSet& mapped_boxes_in_transit,
   const std::vector<int>& peer_ranks) const
{
   t_pack_edge->start();

   const int num_children = static_cast<int>(peer_ranks.size()) - 1;
   const int parent_rank = peer_ranks[num_children];
   const std::vector<int>& child_ranks = peer_ranks;

   hier::NeighborhoodSet unbalanced_eto_balanced(d_dim);
   unbalanced_to_balanced.swapInitialize(
      unbalanced_to_balanced.getBase(),
      unbalanced_to_balanced.getHead(),
      unbalanced_to_balanced.getConnectorWidth(),
      unbalanced_eto_balanced);

   for (TransitSet::const_iterator
        ni = mapped_boxes_in_transit.begin(); ni != mapped_boxes_in_transit.end(); ++ni) {

      BoxInTransit mapped_box_in_transit = *ni;

      // If there is no process history, then local process must be originator.
      if (mapped_box_in_transit.proc_hist.empty()) {

         if (d_print_edge_steps) {
            tbox::plog << "Keeping edge " << mapped_box_in_transit << std::endl;
         }
         TBOX_ASSERT(mapped_box_in_transit.orig_mapped_box.getOwnerRank() == d_mpi.getRank());
         unbalanced_eto_balanced.insertNeighbor(
            mapped_box_in_transit.orig_mapped_box.getId(),
            mapped_box_in_transit.mapped_box);

      } else {

         TBOX_ASSERT(mapped_box_in_transit.orig_mapped_box.getOwnerRank() != d_mpi.getRank());

         // Remember the previous owner and remove it from the proc_hist.
         const int prev_owner = mapped_box_in_transit.proc_hist.back();
         mapped_box_in_transit.proc_hist.pop_back();

         // Find message corresponding to prev_owner.
         int source_i = -1;
         if (parent_rank != -1 && prev_owner == parent_rank) {
            source_i = num_children; // Means source is parent.
         } else {
            for (source_i = 0; source_i != num_children; ++source_i) {
               if (prev_owner == child_ranks[source_i]) {
                  break;
               }
            }
         }

         TBOX_ASSERT(source_i < num_children + 1);

         std::vector<int>& msg = messages[source_i];

         // Pack mapped_box_in_transit into msg.
         if (d_print_edge_steps) {
            tbox::plog << "Pack edge to "
                       << (source_i < num_children ? "child " : "parent ")
                       << (source_i < num_children ? child_ranks[source_i] : parent_rank)
                       << ": " << mapped_box_in_transit
                       << std::endl;
         }
         msg.insert(msg.end(), mapped_box_in_transit.commBufferSize(), 0);
         mapped_box_in_transit.putToIntBuffer(&msg[msg.size()
                                                   - mapped_box_in_transit.commBufferSize()]);
      }

   }
   unbalanced_to_balanced.swapInitialize(
      unbalanced_to_balanced.getBase(),
      unbalanced_to_balanced.getHead(),
      unbalanced_to_balanced.getConnectorWidth(),
      unbalanced_eto_balanced);
   t_pack_edge->stop();
}



/*
 *************************************************************************
 * Unpack BoxInTransit objects from relationship message and either
 * use data to set up local edges in unbalanced--->balanced or route
 * data in the direction of the origin of the BoxInTransit.
 *************************************************************************
 */
void TreeLoadBalancer::unpackAndRouteNeighborhoodSets(
   const int* received_data,
   int received_data_length,
   std::vector<int> outgoing_messages[],
   hier::Connector& unbalanced_to_balanced,
   const std::vector<int>& peer_ranks) const
{
   t_unpack_edge->start();

   const int num_children = static_cast<int>(peer_ranks.size()) - 1;
   const int parent_rank = peer_ranks[num_children];
   const std::vector<int>& child_ranks = peer_ranks;

   hier::NeighborhoodSet unbalanced_eto_balanced(d_dim);
   unbalanced_to_balanced.swapInitialize(
      unbalanced_to_balanced.getBase(),
      unbalanced_to_balanced.getHead(),
      unbalanced_to_balanced.getConnectorWidth(),
      unbalanced_eto_balanced);

   const int* beg = received_data;
   const int* end = received_data + received_data_length;

   while (beg < end) {

      BoxInTransit mapped_box_in_transit(d_dim);
      mapped_box_in_transit.getFromIntBuffer(beg);
      beg += mapped_box_in_transit.commBufferSize();
      if (d_print_edge_steps) {
         tbox::plog << "Unpacked edge " << mapped_box_in_transit << std::endl;
      }

      TBOX_ASSERT(beg <= end);


      /*
       * If mapped_box_in_transit originated on this process (its
       * proc_hist is empty), the relationship should live here.
       * Else, send mapped_box_in_transit toward its originating
       * process by packing it into the message to route it to its
       * previous owner.
       */
      if (mapped_box_in_transit.proc_hist.empty()) {

         if (d_print_edge_steps) {
            tbox::plog << "Keeping edge " << mapped_box_in_transit << std::endl;
         }

         TBOX_ASSERT(mapped_box_in_transit.orig_mapped_box.getOwnerRank() == d_mpi.getRank());

         unbalanced_eto_balanced.insertNeighbor(
            mapped_box_in_transit.orig_mapped_box.getId(),
            mapped_box_in_transit.mapped_box);

      } else {

         TBOX_ASSERT(mapped_box_in_transit.orig_mapped_box.getOwnerRank() != d_mpi.getRank());

         // Remember the previous owner and remove it from the proc_hist.
         const int prev_owner = mapped_box_in_transit.proc_hist.back();
         mapped_box_in_transit.proc_hist.pop_back();

         // Find outgoing message corresponding to prev_owner.
         int source_i = -1;
         if (parent_rank != -1 && prev_owner == parent_rank) {
            source_i = num_children; // Means source is parent.
         } else {
            for (source_i = 0; source_i != num_children; ++source_i) {
               if (prev_owner == child_ranks[source_i]) {
                  break;
               }
            }
         }

         TBOX_ASSERT(source_i < num_children + 1);

         std::vector<int>& msg = outgoing_messages[source_i];

         // Pack mapped_box_in_transit into msg.
         if (d_print_edge_steps) {
            tbox::plog << "Pack edge to "
                       << (source_i < num_children ? "child " : "parent ")
                       << (source_i < num_children ? child_ranks[source_i] : parent_rank)
                       << ": " << mapped_box_in_transit
                       << std::endl;
         }
         msg.insert(msg.end(), mapped_box_in_transit.commBufferSize(), 0);
         mapped_box_in_transit.putToIntBuffer(&msg[msg.size()
                                                   - mapped_box_in_transit.commBufferSize()]);
      }
   }
   unbalanced_to_balanced.swapInitialize(
      unbalanced_to_balanced.getBase(),
      unbalanced_to_balanced.getHead(),
      unbalanced_to_balanced.getConnectorWidth(),
      unbalanced_eto_balanced);
   t_unpack_edge->stop();
}



/*
 *************************************************************************
 * Set the MPI commuicator.  If there's a private communicator, free
 * it first.  It's safe to free the private communicator because no
 * other code have access to it.
 *************************************************************************
 */
void TreeLoadBalancer::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& samrai_mpi)
{
   if (samrai_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      TBOX_ERROR("TreeLoadBalancer::setMPICommunicator error: Given\n"
         << "communicator is invalid.");
   }

   d_mpi_dup.freeCommunicator();

   // Enable private communicator.
   d_mpi_dup.dupCommunicator(samrai_mpi);
}



/*
 *************************************************************************
 * Set the MPI commuicator.
 *************************************************************************
 */
void TreeLoadBalancer::freeMPICommunicator()
{
   // Free the private communicator (if MPI has not been finalized).
   int flag;
   tbox::SAMRAI_MPI::Finalized(&flag);
   if (!flag) {
      d_mpi_dup.freeCommunicator();
   }
}



/*
 *************************************************************************
 * Divide all ranks into groups.  With each cycle,
 * the group size grows geometrically while the number of groups
 * shrink geometrically.
 * The last cycle_number has a single group of d_mpi.getSize() processors.
 *
 * Let m = number of cycles
 * i = cycle number, [0,m)
 * p = communicator size
 * r = rank in communicator
 * q = group size
 * c = number of groups
 * g = group id
 * s = rank in group
 * d = tree degree (d==2 means binary tree)
 *
 * c = p^( 1-(i+1)/m )
 *
 * With interlacing
 * g = r%q
 * r = q*s + g
 * s = (r+g)/q
 *
 * Without interlacing
 * q = ceil(p/c)
 * g=r/q
 * s = r%c
 *
 * Note that for the last cycle, p == q and r == s.
 *
 * The tree is built from the processor ranks in its group, as shown:
 *         0
 *        / \
 *       /   \
 *      /     \
 *     1       2
 *    / \     / \
 *   3   4   5   6
 *  /|  /|   |
 * 7 8 9 10 11 ...
 *
 * Communication with children are separated by communication
 * with parents, so we can reuse child_comms.  Sends and
 * receives with parents can overlap, so we have two objects
 * (parent_send and parent_recv) whose activities may
 * overlap.
 *
 *************************************************************************
 */
void TreeLoadBalancer::createBalanceSubgroups(
   const int cycle_number,
   const int number_of_cycles,
   const tbox::RankGroup& rank_group,
   int& g_count,
   int& g_id,
   int& g_rank,
   int& g_nproc,
   int& num_children,
   tbox::AsyncCommPeer<int> *& child_comms,
   tbox::AsyncCommPeer<int> *& parent_send,
   tbox::AsyncCommPeer<int> *& parent_recv,
   tbox::AsyncCommStage& comm_stage) const
{
   if (d_mpi.getSize() == 1) {
      child_comms = parent_send = parent_recv = NULL;
      g_count = 1;
      g_nproc = d_mpi.getSize();
      g_id = 0;
      g_rank = 0;
      return;
   }

   TBOX_ASSERT(d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull);
   TBOX_ASSERT(child_comms == NULL);
   TBOX_ASSERT(parent_send == NULL);
   TBOX_ASSERT(parent_recv == NULL);

   /*
    * Groups can be contiguous intervals in [0..nproc) or they can be
    * interlaced.
    *
    * We assume that contiguous numbers of processors tend to be closer
    * on the network.  We usually want contiguous cycles to help reduce
    * the communication hops made by transfered loads.  However,
    * multiple cycles, used for severe initial imbalances, are meant to
    * spread the load across the machine as fast as possible so we want
    * contiguous groups only in the final cycle.
    */
   const bool interlace_groups = cycle_number < number_of_cycles - 1;

   TBOX_ASSERT(rank_group.isMember(d_mpi.getRank()));
   int my_output_id = rank_group.getMapIndex(d_mpi.getRank());
   int output_nproc = rank_group.size();

   /*
    * Compute group sizes, membership, etc.
    *
    * Tiny groups tend to leave the members with
    * possibly large overloads.
    * In order to make all groups similar in size
    * we tend to round down the number of groups
    * (and round up the group size).
    */
   // Number of groups:
   g_count = (int)pow((double)output_nproc,
         1.0 - double(cycle_number + 1) / number_of_cycles);
   if (interlace_groups) {
      // Group id (group number):
      g_id = my_output_id % g_count;
      // Rank in group g_id:
      g_rank = my_output_id / g_count;
      // Population of group g_id:
      g_nproc = output_nproc / g_count + (g_id < (output_nproc % g_count));
   } else {
      /*
       * All groups will have a nominal population count.
       * Remainder is distributed to a subset of groups, starting from group 0,
       * so these groups will have one more than nominal.
       */
      int nominal_g_nproc = output_nproc / g_count;
      int first_nominal_group = output_nproc % g_count;
      int first_nominal_rank_id = first_nominal_group * (1 + nominal_g_nproc);
      if (my_output_id < first_nominal_rank_id) {
         // Group id (group number):
         g_id = my_output_id / (1 + nominal_g_nproc);
         // Rank in group g_id:
         g_rank = my_output_id % (1 + nominal_g_nproc);
         // Population of group g_id:
         g_nproc = nominal_g_nproc + 1;
      } else {
         // Group id (group number):
         g_id = first_nominal_group
            + (my_output_id - first_nominal_rank_id) / nominal_g_nproc;
         // Rank in group g_id:
         g_rank = (my_output_id - first_nominal_rank_id) % nominal_g_nproc;
         // Population of group g_id:
         g_nproc = nominal_g_nproc;
      }
   }

#if 1
   /*
    * Arrange the group ranks in a BalancedDepthFirstTree in order to get
    * the parent/children in the group.
    *
    * By the way, the BalancedDepthFirstTree currently assumes a binary
    * tree, d_degree = 2.
    */
   TBOX_ASSERT(d_degree == 2);
   tbox::BalancedDepthFirstTree bdfs(0, g_nproc - 1, g_rank, true);
   num_children = bdfs.getNumberOfChildren();
   child_comms = new tbox::AsyncCommPeer<int>[num_children];
   for (int i = 0; i < num_children; ++i) {
      const int child_g_rank = bdfs.getChildRank(i);
      int child_true_rank;
      if (interlace_groups) {
         child_true_rank = rank_group.getMappedRank(child_g_rank * g_count + g_id);
      } else {
         int nominal_g_nproc = output_nproc / g_count;
         int first_nominal_group = output_nproc % g_count;
         int group_first_rank = g_id < first_nominal_group ?
            g_id * (1 + nominal_g_nproc) :
            first_nominal_group
            * (1
               + nominal_g_nproc)
            + (g_id - first_nominal_group) * nominal_g_nproc;
         child_true_rank = rank_group.getMappedRank(child_g_rank + group_first_rank);
      }
      child_comms[i].initialize(&comm_stage);
      child_comms[i].setPeerRank(child_true_rank);
      child_comms[i].setMPI(d_mpi);
      child_comms[i].setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      child_comms[i].limitFirstDataLength(TreeLoadBalancer_FIRSTDATALEN);
   }
   if (bdfs.getParentRank() != bdfs.getInvalidRank()) {
      const int parent_g_rank = bdfs.getParentRank();
      int parent_true_rank;
      if (interlace_groups) {
         parent_true_rank = rank_group.getMappedRank(parent_g_rank * g_count + g_id);
      } else {
         int nominal_g_nproc = output_nproc / g_count;
         int first_nominal_group = output_nproc % g_count;
         int group_first_rank = g_id < first_nominal_group ?
            g_id * (1 + nominal_g_nproc) :
            first_nominal_group
            * (1
               + nominal_g_nproc)
            + (g_id - first_nominal_group) * nominal_g_nproc;
         parent_true_rank = rank_group.getMappedRank(parent_g_rank + group_first_rank);
      }

      parent_send = new tbox::AsyncCommPeer<int>;
      parent_send->initialize(&comm_stage);
      parent_send->setPeerRank(parent_true_rank);
      parent_send->setMPI(d_mpi);
      parent_send->setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      parent_send->limitFirstDataLength(TreeLoadBalancer_FIRSTDATALEN);

      parent_recv = new tbox::AsyncCommPeer<int>;
      parent_recv->initialize(&comm_stage);
      parent_recv->setPeerRank(parent_true_rank);
      parent_recv->setMPI(d_mpi);
      parent_recv->setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      parent_recv->limitFirstDataLength(TreeLoadBalancer_FIRSTDATALEN);
   }
#else
   /*
    * Arrange the group ranks in a breadth-first tree in order to get
    * the parent/children in the group.
    */
   for (num_children = 0; num_children < d_degree; ++num_children) {
      const int child_g_rank = d_degree * (g_rank + 1) + num_children - 1;
      const int child_true_rank = child_g_rank * g_count + g_id;
      if (child_true_rank >= d_mpi.getSize()) break;
   }
   child_comms = new tbox::AsyncCommPeer<int>[num_children];
   for (int i = 0; i < num_children; ++i) {
      const int child_g_rank = d_degree * (g_rank + 1) + i - 1;
      int child_true_rank;
      if (interlace_groups) {
         child_true_rank = child_g_rank * g_count + g_id;
      } else {
         int nominal_g_nproc = d_mpi.getSize() / g_count;
         int first_nominal_group = d_mpi.getSize() % g_count;
         // int first_nominal_rank = first_nominal_group*(1+nominal_g_nproc);
         int group_first_rank = g_id < first_nominal_group ?
            g_id * (1 + nominal_g_nproc) :
            first_nominal_group
            * (1
               + nominal_g_nproc)
            + (g_id - first_nominal_group) * nominal_g_nproc;
         child_true_rank = child_g_rank + group_first_rank;
      }
      child_comms[i].initialize(&comm_stage);
      child_comms[i].setPeerRank(child_true_rank);
      child_comms[i].setSAMRAI_MPI(d_mpi);
      child_comms[i].setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      child_comms[i].limitFirstDataLength(TreeLoadBalancer_FIRSTDATALEN);
   }
   if (g_rank > 0) {
      const int parent_g_rank = (g_rank + 1) / d_degree - 1;
      int parent_true_rank;
      if (interlace_groups) {
         parent_true_rank = parent_g_rank * g_count + g_id;
      } else {
         int nominal_g_nproc = d_mpi.getSize() / g_count;
         int first_nominal_group = d_mpi.getSize() % g_count;
         // int first_nominal_rank = first_nominal_group*(1+nominal_g_nproc);
         int group_first_rank = g_id < first_nominal_group ?
            g_id * (1 + nominal_g_nproc) :
            first_nominal_group
            * (1
               + nominal_g_nproc)
            + (g_id - first_nominal_group) * nominal_g_nproc;
         parent_true_rank = parent_g_rank + group_first_rank;
      }
      parent_send = new tbox::AsyncCommPeer<int>;
      parent_send->initialize(&comm_stage);
      parent_send->setPeerRank(parent_true_rank);
      parent_send->setSAMRAI_MPI(d_mpi);
      parent_send->setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      parent_send->limitFirstDataLength(TreeLoadBalancer_FIRSTDATALEN);

      parent_recv = new tbox::AsyncCommPeer<int>;
      parent_recv->initialize(&comm_stage);
      parent_recv->setPeerRank(parent_true_rank);
      parent_recv->setSAMRAI_MPI(d_mpi);
      parent_recv->setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      parent_recv->limitFirstDataLength(TreeLoadBalancer_FIRSTDATALEN);
   }
#endif
   if (d_print_steps) {
      tbox::plog << "Groups for cycle " << cycle_number
                 << '/' << number_of_cycles
                 << ": g_count=" << g_count
                 << "  g_nproc=" << g_nproc
                 << "  g_id=" << g_id
                 << "  g_rank=" << g_rank
                 << "  parent="
                 << (parent_send ? parent_send->getPeerRank() : -1)
                 << "  children (" << num_children << ") =";
      for (int i = 0; i < num_children; ++i) {
         tbox::plog << ' ' << child_comms[i].getPeerRank();
      }
      tbox::plog << "\n";
   }

}



/*
 *************************************************************************
 *************************************************************************
 */
void TreeLoadBalancer::destroyBalanceSubgroups(
   tbox::AsyncCommPeer<int> *& child_comms,
   tbox::AsyncCommPeer<int> *& parent_send,
   tbox::AsyncCommPeer<int> *& parent_recv) const
{
   if (d_mpi.getSize() == 1) {
      TBOX_ASSERT(child_comms == NULL);
      TBOX_ASSERT(parent_send == NULL);
      TBOX_ASSERT(parent_recv == NULL);
   } else {
      delete[] child_comms;
      delete parent_send;
      delete parent_recv;
      child_comms = parent_send = parent_recv = NULL;
   }
}



/*
 *************************************************************************
 * Sort the sizes of an IntVector from smallest to largest value.
 *************************************************************************
 */
void TreeLoadBalancer::sortIntVector(
   hier::IntVector& ordered_dims,
   const hier::IntVector& vector) const
{
   const hier::IntVector num_cells = vector;

   for (int d = 0; d < d_dim.getValue(); d++) {
      ordered_dims(d) = d;
   }
   for (int d0 = 0; d0 < d_dim.getValue() - 1; d0++) {
      for (int d1 = d0 + 1; d1 < d_dim.getValue(); d1++) {
         if (vector(ordered_dims(d0)) > vector(ordered_dims(d1))) {
            int tmp_d = ordered_dims(d0);
            ordered_dims(d0) = ordered_dims(d1);
            ordered_dims(d1) = tmp_d;
         }
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int d = 0; d < d_dim.getValue() - 1; d++) {
      TBOX_ASSERT(vector(ordered_dims(d)) <= vector(ordered_dims(d + 1)));
   }
#endif
}



/*
 *************************************************************************
 * Attempt to shift a specified ammount of load from one TransitSet to
 * another by breaking a single box from the overloaded set.  Examine
 * multiple NodeInTransit and breakages to
 *
 * - filter out the ones that worsens the balance (see breakOffLoad) for
 * reasons why this can happen.
 *
 * - find the one that results in the smallest combinedBreakingPenalty.
 *
 * Return whether any changes were made.
 *************************************************************************
 */
bool TreeLoadBalancer::shiftLoadsByBreaking(
   TransitSet& src,
   TransitSet& dst,
   const int ideal_transfer,
   int& actual_transfer,
   hier::LocalId& next_available_index) const
{
   if (ideal_transfer < 0) {
      bool rval = shiftLoadsByBreaking(dst,
            src,
            -ideal_transfer,
            actual_transfer,
            next_available_index);
      actual_transfer = -actual_transfer;
      return rval;
   }

   TBOX_ASSERT(src.size() + dst.size() > 0);

   t_shift_loads_by_breaking->start();

   if (d_print_steps) {
      tbox::plog << "    shiftLoadsByBreaking asked to break off "
                 << ideal_transfer
                 << " from one of " << src.size()
                 << " source Boxes to add to set of " << dst.size()
                 << " Boxes."
                 << std::endl;
   }

   actual_transfer = 0;

   /*
    * The best results so far (from not transfering anything).
    */
   TransitSet best_src = src;
   TransitSet best_dst = dst;
   int best_actual_transfer = 0;

   double best_combined_penalty = combinedBreakingPenalty(
         computeBalancePenalty(best_src, best_dst, (double)ideal_transfer
            - best_actual_transfer),
         computeSurfacePenalty(best_src, best_dst),
         computeSlenderPenalty(best_src, best_dst));

   /*
    * Scale the pre-cut penalty.  Scaling makes this method more
    * agressive about producing a cut.  Sometimes the uncut penalty can
    * be low enough to prevent a cut, but the cut may be required for
    * reasonable balancing.  The down side of agressive cutting is that
    * it tends to produce more slivers (but not terribly many).
    * This exchange of load balance and slivers may not be easily
    * eliminated by the weights in the penalty functions.
    */
   if (d_print_steps) {
      tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
      tbox::plog.precision(6);
      tbox::plog << "    Uncut penalty: " << best_combined_penalty
                 << " scaled by " << d_precut_penalty_wt << " to "
                 << best_combined_penalty * d_precut_penalty_wt
                 << std::endl;
   }
   best_combined_penalty *= d_precut_penalty_wt;

   bool found_breakage = false;

   std::vector<hier::Box> breakoff;
   std::vector<hier::Box> leftover;
   double breakoff_amt;
   for (TransitSet::iterator si = src.begin();
        si != src.end() &&
        (!found_breakage || si->load >= ideal_transfer); ++si) {

      const BoxInTransit& candidate = *si;

      breakOffLoad(
         candidate.mapped_box,
         ideal_transfer,
         breakoff,
         leftover,
         breakoff_amt);

      if (!breakoff.empty()) {

         const BoxInTransit& brk_mapped_box_in_transit = *si;

         const bool improves_balance =
            tbox::MathUtilities<double>::Abs(
               (double)ideal_transfer - breakoff_amt) <
            (ideal_transfer - tbox::MathUtilities<double>::getEpsilon());

         if (d_print_steps) {
            tbox::plog << "    Potential to replace " << brk_mapped_box_in_transit << " with "
                       << breakoff.size() << " breakoff Boxes and "
                       << leftover.size() << " leftover Boxes."
                       << "  improves_balance=" << improves_balance
                       << std::endl;
         }

         if (!improves_balance) {
            /*
             * Consider only options that changes the balance for the
             * better.  Options that degrades or leave alone the
             * balance put us in an infinite loop bouncing between
             * better balance and better combined penalty.
             */
            continue;
         }

         /*
          * Trial modifications of src and dst, for evaluating penalties.
          */
         TransitSet trial_src = src;
         TransitSet trial_dst = dst;
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         int trial_actual_transfer = actual_transfer;

         /*
          * Put breakoff in trial_dst and put leftover back into trial_src.
          */
         for (std::vector<hier::Box>::const_iterator bi = breakoff.begin();
              bi != breakoff.end();
              ++bi) {
            BoxInTransit give_mapped_box_in_transit(
               brk_mapped_box_in_transit,
               *bi,
               d_mpi.getRank(),
               next_available_index);
            give_mapped_box_in_transit.load = (int)computeLoad(
                  give_mapped_box_in_transit.orig_mapped_box,
                  give_mapped_box_in_transit.getBox());
            next_available_index += 2 + d_degree;
            trial_dst.insert(give_mapped_box_in_transit);
            trial_actual_transfer += give_mapped_box_in_transit.load;
            if (d_print_steps) {
               tbox::plog << "    Breakoff box " << *bi << " -> " << give_mapped_box_in_transit
                          << std::endl;
            }
         }
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         for (std::vector<hier::Box>::const_iterator bi = leftover.begin();
              bi != leftover.end();
              ++bi) {
            BoxInTransit keep_mapped_box_in_transit(
               brk_mapped_box_in_transit,
               *bi,
               d_mpi.getRank(),
               next_available_index);
            keep_mapped_box_in_transit.load = (int)computeLoad(
                  keep_mapped_box_in_transit.orig_mapped_box,
                  keep_mapped_box_in_transit.getBox());
            next_available_index += 2 + d_degree;
            trial_src.insert(keep_mapped_box_in_transit);
            if (d_print_steps) {
               tbox::plog << "    Leftover box " << *bi << " -> " << keep_mapped_box_in_transit
                          << std::endl;
            }
         }
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         trial_src.erase(brk_mapped_box_in_transit);

         /*
          * Compute the new penalty to see if it improves our best result so far.
          */
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         double trial_balance_penalty = computeBalancePenalty(trial_src,
               trial_dst,
               (double)ideal_transfer - trial_actual_transfer);
         double trial_surface_penalty = computeSurfacePenalty(trial_src,
               trial_dst);
         double trial_slender_penalty = computeSlenderPenalty(trial_src,
               trial_dst);
         double trial_combined_penalty = combinedBreakingPenalty(
               trial_balance_penalty,
               trial_surface_penalty,
               trial_slender_penalty);
         if (d_print_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "    Trial's imbalance: "
                       << (ideal_transfer - trial_actual_transfer)
                       << " balance,surface,slender,combined penalty: "
                       << trial_balance_penalty << ' '
                       << trial_surface_penalty << ' '
                       << trial_slender_penalty << ' '
                       << trial_combined_penalty
                       << std::endl;
         }

         if (trial_combined_penalty < best_combined_penalty) {
            if (d_print_steps) {
               tbox::plog << "    Keeping this trial." << std::endl;
            }
            found_breakage = true;
            best_actual_transfer = (int)breakoff_amt;
            best_src = trial_src;
            best_dst = trial_dst;
            best_combined_penalty = trial_combined_penalty;
         } else {
            if (d_print_steps) {
               tbox::plog << "    Rejecting this trial." << std::endl;
            }
         }

      } else {
         if (d_print_steps) {
            tbox::plog << "    Break step could not break " << ideal_transfer
                       << std::endl;
         }
      }

   }

   if (found_breakage) {
      src.swap(best_src);
      dst.swap(best_dst);
      actual_transfer = best_actual_transfer;
   }

   t_shift_loads_by_breaking->stop();
   return found_breakage;
}



/*
 *************************************************************************
 * Attempt to swap some NodesInTransit between 2 sets of
 * NodesInTransit (src and dst) to shift ideal_transfer work units.
 *
 * Transfering a BoxInTransit from one TransitSet to another
 * is considered a degenerate "swap" (a BoxInTransit is
 * swapped for nothing) handled by this function.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of dst
 * ideal_transfer < 0 means to raise the load of src
 * The iterative do loop may overshoot the ideal_transfer
 * and may have to swap to shift some of the load
 * back.
 *
 * Return whether any changes were made.
 *************************************************************************
 */
bool TreeLoadBalancer::shiftLoadsBySwapping(
   TransitSet& src,
   TransitSet& dst,
   int ideal_transfer,
   int& actual_transfer) const
{
   t_shift_loads_by_swapping->start();

   if (d_print_steps) {
      tbox::plog << "  Attempting to swap " << ideal_transfer << std::endl;
   }

   bool return_value = false;
   bool found_swap;

   actual_transfer = 0;

   do {

      /*
       * Ammount we seek to transfer from hi to lo
       * (the "ideal" for this particular iteration).
       * Unlike ideal_transfer and actual_transfer, this quantity is positive.
       */
      int rem_transfer = ideal_transfer - actual_transfer;
      if (d_print_steps) {
         tbox::plog << "    Swap progress: " << actual_transfer
                    << " / " << ideal_transfer << " remaining transfer = "
                    << rem_transfer << std::endl;
      }

      found_swap = false;

      TransitSet::iterator isrc;
      TransitSet::iterator idst;
      int swap_transfer;
      found_swap = findLoadSwapPair(src,
            dst,
            rem_transfer,
            swap_transfer,
            isrc,
            idst);

      if (found_swap) {

         // We can improve balance_penalty by swapping isrc with idst.
         if (d_print_steps) {
            tbox::plog << "    Swapping " << swap_transfer << " units using ";
            if (isrc != src.end()) tbox::plog << *isrc;
            else tbox::plog << "X";
            tbox::plog << " <--> ";
            if (idst != dst.end()) tbox::plog << *idst;
            else tbox::plog << "X";
            tbox::plog << std::endl;
         }

         if (isrc != src.end()) {
            dst.insert(*isrc);
            src.erase(isrc);
         }
         if (idst != dst.end()) {
            src.insert(*idst);
            dst.erase(idst);
         }

         actual_transfer += swap_transfer;
         return_value = true;

      } else {
         if (d_print_steps) {
            tbox::plog << "    Cannot find swap pair for " << rem_transfer
                       << " units." << std::endl;
         }
      }

   } while (found_swap && (actual_transfer != ideal_transfer));

   if (d_print_steps) {
      tbox::plog << "  Final imbalance for swap " << ideal_transfer
      - actual_transfer << std::endl;
   }

   t_shift_loads_by_swapping->stop();

   return return_value;
}



/*
 *************************************************************************
 * Find a BoxInTransit in src and a BoxInTransit in dst which when swapped
 * results in shifting close to ideal_shift from src to dst.
 * Return whether a swap pair was found.
 *
 * If isrc is set but idst is not, it means that isrc should
 * be moved to dst, but no dst should be moved back.  This is
 * the degenerate case of swapping isrc for a BoxInTransit with zero
 * load.
 *************************************************************************
 */
bool TreeLoadBalancer::findLoadSwapPair(
   TransitSet& src,
   TransitSet& dst,
   int ideal_transfer,
   int& actual_transfer,
   TransitSet::iterator& isrc,
   TransitSet::iterator& idst) const
{
   if (ideal_transfer < 0) {
      // The logic below does not handle bi-directional so handle it here.
      bool rval = findLoadSwapPair(dst,
            src,
            -ideal_transfer,
            actual_transfer,
            idst,
            isrc);
      actual_transfer = -actual_transfer;
      return rval;
   }

   t_find_swap_pair->start();

   if (d_print_steps) {
      tbox::plog << "  findLoadSwapPair looking for transfer of "
                 << ideal_transfer
                 << " between " << src.size() << "-BoxInTransit src and "
                 << dst.size() << "-BoxInTransit dst." << std::endl;
   }
   if (d_print_swap_steps) {
      tbox::plog << "    src (" << src.size() << "):" << std::endl;
      for (TransitSet::iterator si = src.begin(); si != src.end(); ++si) {
         tbox::plog << *si << std::endl;
      }
      tbox::plog << "    dst (" << dst.size() << "):" << std::endl;
      for (TransitSet::iterator si = dst.begin(); si != dst.end(); ++si) {
         tbox::plog << *si << std::endl;
      }
   }

   /*
    * Look for a high and a low NodeInTransit to swap, subject to condition
    *
    * 0 < load(high) - load(low) < lim_transfer.
    *
    * Compute the
    * balance_penalty if high and low were swapped.  Keep looking
    * until we find the pair giving the lowest balance_penalty on swapping.
    *
    * isrc and idst point to the current best pair to swap.  new_balance_penalty
    * is the balance_penalty if we swap them.
    *
    * jsrc and jdst are trial pairs to check to see if we can improve on
    * new_balance_penalty.
    *
    * We will look for two "best" pairs:
    *
    * On the high side, swapping more than ideal_transfer:
    * (jsrc_hiside , jdst_hiside) are the best swap pair such that
    * jsrc_hiside->load - jdst_hiside->load >= ideal_transfer
    *
    * On the low side, swapping less than ideal_transfer:
    * (jsrc_loside, jdst_loside are the best swap pair such that
    * jsrc_loside->load - jdst_loside->load < ideal_transfer
    *
    * Then we decide on whether to choose the hig-hside swap or the
    * low-side swap based on balance_penalty of each candidate swap.
    */

   // Initialization indicating no swap pair found.
   TransitSet::iterator jsrc_hiside = src.end();
   TransitSet::iterator jdst_hiside = dst.end();
   TransitSet::iterator jsrc_loside = src.end();
   TransitSet::iterator jdst_loside = dst.end();

   // A dummy BoxInTransit for set searches.
   hier::Box dummy_box(d_dim);
   BoxInTransit dummy_search_target(d_dim);

   // Distance between ideal and candidate, >= 0
   int imbalance_loside = tbox::MathUtilities<int>::getMax();
   int imbalance_hiside = tbox::MathUtilities<int>::getMax();

   if (dst.empty()) {
      /*
       * There is no dst BoxInTransit, so the swap would
       * degnerate to moving a src BoxInTransit to dst.  Find
       * the best src BoxInTransit to move.
       */
      dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
      dummy_search_target.load = ideal_transfer;
      TransitSet::iterator jsrc = src.lower_bound(dummy_search_target);

      if (d_print_swap_steps) {
         tbox::plog << "  findLoadSwapPair with empty dst: ";
      }
      if (jsrc != src.begin()) {
         jsrc_hiside = jsrc;
         --jsrc_hiside;
         imbalance_hiside = jsrc_hiside->load - ideal_transfer;
         if (d_print_swap_steps) {
            tbox::plog << "  hi src: " << (*jsrc_hiside)
                       << " with transfer " << (*jsrc_hiside).load
                       << " missing by " << imbalance_hiside;
         }
      }

      if (jsrc != src.end()) {
         jsrc_loside = jsrc;
         imbalance_loside = ideal_transfer - jsrc_loside->load;
         if (d_print_swap_steps) {
            tbox::plog << "  lo src: " << (*jsrc_loside)
                       << " with transfer " << (*jsrc_loside).load
                       << " missing by " << imbalance_loside;
         }
      }
      if (d_print_swap_steps) {
         tbox::plog << std::endl;
      }

   } else {

      /*
       * Start search through src beginning with the NodeInTransit
       * whose load exceed the biggest dst BoxInTransit by at
       * least ideal_transfer.
       */
      dummy_search_target = *dst.begin();
      dummy_search_target.load += ideal_transfer;
      TransitSet::iterator jsrc_beg = src.lower_bound(dummy_search_target);

      for (TransitSet::iterator jsrc = jsrc_beg; jsrc != src.end(); ++jsrc) {

         /*
          * Set jdst pointing to where we should start looking in dst.
          * Look for a load less than the load of jsrc by
          * ideal_transfer.
          */
         dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
         dummy_search_target.load = tbox::MathUtilities<int>::Max(
               jsrc->load - ideal_transfer,
               0);
         TransitSet::iterator jdst = dst.lower_bound(dummy_search_target);

         if (jdst != dst.end()) {

            /*
             * lower_bound returns jdst that gives >= ideal_transfer
             * when swapped with jsrc.  Check transfererence between
             * jsrc and two dst NodeInTransit on either side of the
             * ideal dst BoxInTransit to swap.
             */
            int tmp_miss = (jsrc->load - jdst->load) - ideal_transfer;
            TBOX_ASSERT(tmp_miss >= 0);
            if ((tmp_miss < imbalance_hiside)) {
               jsrc_hiside = jsrc;
               jdst_hiside = jdst;
               imbalance_hiside = tmp_miss;
               if (d_print_swap_steps) {
                  tbox::plog << "    new hi-swap pair: " << (*jsrc_hiside)
                             << " & " << (*jdst_hiside) << " with transfer "
                             << (jsrc_hiside->load - jdst_hiside->load)
                             << " missing by " << imbalance_hiside
                             << std::endl;
               }
            }

            if (jdst != dst.begin()) {
               --jdst; // Now, jsrc and jdst transferer by *less* than ideal_transfer.
               tmp_miss = ideal_transfer - (jsrc->load - jdst->load);
               TBOX_ASSERT(tmp_miss >= 0);
               if (tmp_miss < imbalance_loside) {
                  jsrc_loside = jsrc;
                  jdst_loside = jdst;
                  imbalance_loside = tmp_miss;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap pair: " << (*jsrc_loside)
                                << " & " << (*jdst_loside) << " with transfer "
                                << (jsrc_loside->load - jdst_loside->load)
                                << " missing by " << imbalance_loside
                                << std::endl;
                  }
               }
            }

         } else {

            /*
             * The ideal dst to swap is smaller than the smallest dst
             * BoxInTransit.  Check the transfererence between
             * jsrc and smallest BoxInTransit, the case of
             * moving jsrc in whole.
             */
            if (jsrc->load > ideal_transfer) {
               // Moving jsrc to src is moving too much--hiside.
               int tmp_miss = jsrc->load - ideal_transfer;
               TBOX_ASSERT(tmp_miss >= 0);
               if (tmp_miss < imbalance_hiside) {
                  jsrc_hiside = jsrc;
                  jdst_hiside = dst.end();
                  imbalance_hiside = tmp_miss;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new hi-swap source: " << (*jsrc_hiside)
                                << " & " << "no dst" << " with transfer "
                                << (jsrc_hiside->load)
                                << " missing by " << imbalance_hiside
                                << std::endl;
                  }
               }
            } else {
               // Moving jsrc to src is moving (just right or) too little--loside.
               int tmp_miss = ideal_transfer - jsrc->load;
               TBOX_ASSERT(tmp_miss >= 0);
               if (tmp_miss < imbalance_loside) {
                  jsrc_loside = jsrc;
                  jdst_loside = dst.end();
                  imbalance_loside = tmp_miss;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap source: " << (*jsrc_loside)
                                << " & " << "no dst" << " with transfer "
                                << (jsrc_loside->load)
                                << " missing by " << imbalance_loside
                                << std::endl;
                  }
               }
               /*
                * Break out of the loop early because there is no
                * point checking smaller src NodeInTransit.
                */
               break;
            }
         }
      }

   }

   /*
    * Swapping does not produce new cuts, so it is ok to omit the penalties
    * arising from cutting.
    */
   double current_balance_penalty = (double)ideal_transfer;
   double balance_penalty_loside = (double)imbalance_loside;
   double balance_penalty_hiside = (double)imbalance_hiside;

   if (d_print_swap_steps) {
      tbox::plog << "    Swap candidates give penalties (unswap,lo,hi): "
                 << current_balance_penalty << " , " << balance_penalty_loside
                 << " , " << balance_penalty_hiside << std::endl;
   }

   bool return_value = false;
   if (balance_penalty_loside < current_balance_penalty &&
       balance_penalty_loside <= balance_penalty_hiside) {
      isrc = jsrc_loside;
      idst = jdst_loside;
      return_value = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking loside." << std::endl;
      }
   } else if (balance_penalty_hiside < current_balance_penalty &&
              balance_penalty_hiside <= balance_penalty_loside) {
      isrc = jsrc_hiside;
      idst = jdst_hiside;
      return_value = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking hiside." << std::endl;
      }
   } else {
      if (d_print_swap_steps) {
         tbox::plog << "    Keeping original (no swap)." << std::endl;
      }
   }

   if (return_value) {
      actual_transfer = (*isrc).load;
      if (idst != dst.end()) {
         actual_transfer -= (*idst).load;
      }
   }

   t_find_swap_pair->stop();
   return return_value;
}



/*
 *************************************************************************
 * Master method for breaking off a load.
 *
 * Try different heuristics and pick the "best" way to break off a load.
 * The best is defined as the one with the lowest combined penalty.
 *
 * Return whether a successful break was made.
 *************************************************************************
 */
bool TreeLoadBalancer::breakOffLoad(
   const hier::Box& mapped_box,
   double ideal_load_to_break,
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load) const
{
   TBOX_ASSERT(ideal_load_to_break > 0);

   /*
    * NOTE: We need in this method a way to weigh the
    * value of proximity to the ideal breakoff vs the
    * increased area of the cuts.  However, the weight
    * given to area-optimized cuts should be considered
    * only with real application performance data.
    *
    * NOTE: We can compute the amount of new box boundaries
    * generated by computing the box boundary before and
    * after, and subtracting.  Easier than reconstructing
    * the cuts from the box definitions.
    *
    * NOTE: We should weight negatively the production of
    * high surface-to-volume boxes.
    */

   t_break_off_load->start();

   breakoff.clear();
   leftover.clear();

   if (ideal_load_to_break < d_min_size.getProduct()) {
      /*
       * Assuming uniform load balancing, there is no hope
       * if breaking off a piece of the desired size.
       */
      if (d_print_break_steps) {
         tbox::plog << "      ideal_load_to_break " << ideal_load_to_break
                    << " < " << d_min_size.getProduct() << d_min_size
                    << ":  Cannot break Box " << mapped_box << std::endl;
      }
      t_break_off_load->stop();
      return false;
   }

   /*
    * To avoid repeated computations of bad cuts,
    * we precompute bad_cuts here to provide to
    * methods that actually use the information.
    */
   tbox::Array<tbox::Array<bool> > bad_cuts(d_dim.getValue());
   t_find_bad_cuts->start();
   hier::BoxUtilities::findBadCutPoints(bad_cuts,
      mapped_box,
      d_domain_boxes,
      d_bad_interval);
   t_find_bad_cuts->stop();

   // Penalty for not transfering ideal load.
   double best_balance_penalty = computeBalancePenalty(mapped_box,
         ideal_load_to_break);
   // Penalty for new surfaces generated (none generated yet).
   double best_surface_penalty = computeSurfacePenalty(mapped_box);
   // Penalty for slender boxes.
   double best_slender_penalty = computeSlenderPenalty(mapped_box);

   double best_combined_penalty = tbox::MathUtilities<double>::getMax();

   if (d_print_break_steps) {
      tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
      tbox::plog.precision(6);
      tbox::plog << "      pre-break imbalance: " << ideal_load_to_break
                 << " balance,surface,slender,combined penalties: "
                 << best_balance_penalty << ", " << best_surface_penalty
                 << ", "
                 << best_slender_penalty << ", " << best_combined_penalty
                 << std::endl;
   }

   brk_load = 0;
   bool found_any_break = false;

   {
      std::vector<hier::Box> planar_breakoff;
      std::vector<hier::Box> planar_leftover;
      double planar_brk_load;

      bool found_this_break = breakOffLoad_planar(
            mapped_box,
            ideal_load_to_break,
            bad_cuts,
            planar_breakoff,
            planar_leftover,
            planar_brk_load);

      if (found_this_break) {
         found_any_break = true;
         double planar_balance_penalty = computeBalancePenalty(planar_breakoff,
               planar_leftover,
               (double)(planar_brk_load - ideal_load_to_break));
         double planar_surface_penalty = computeSurfacePenalty(planar_breakoff,
               planar_leftover);
         double planar_slender_penalty = computeSlenderPenalty(planar_breakoff,
               planar_leftover);
         double planar_combined_penalty =
            combinedBreakingPenalty(planar_balance_penalty,
               planar_surface_penalty,
               planar_slender_penalty);
         if (d_print_break_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "      Planar-break broke off "
                       << planar_brk_load << " / " << ideal_load_to_break
                       << " from " << mapped_box << '|'
                       << mapped_box.numberCells() << '|'
                       << mapped_box.size() << " into "
                       << planar_breakoff.size()
                       << " breakoff: ";
            for (std::vector<hier::Box>::const_iterator bi =
                    planar_breakoff.begin();
                 bi != planar_breakoff.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        and " << planar_leftover.size()
                       << " leftover boxes:";
            for (std::vector<hier::Box>::const_iterator bi =
                    planar_leftover.begin();
                 bi != planar_leftover.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        imbalance: "
                       << (planar_brk_load - ideal_load_to_break)
                       << " balance,surface,slender,combined penalties: "
                       << planar_balance_penalty << ", "
                       << planar_surface_penalty << ", "
                       << planar_slender_penalty
                       << ", " << planar_combined_penalty << std::endl;
         }
         if (planar_combined_penalty < best_combined_penalty) {
            if (d_print_break_steps) {
               tbox::plog << "      Keeping planar cut result." << std::endl;
            }
            breakoff.swap(planar_breakoff);
            leftover.swap(planar_leftover);
            brk_load = planar_brk_load;
            best_balance_penalty = planar_balance_penalty;
            best_surface_penalty = planar_surface_penalty;
            best_slender_penalty = planar_slender_penalty;
            best_combined_penalty = planar_combined_penalty;
         } else {
            if (d_print_break_steps) {
               tbox::plog << "      Rejecting planar cut result." << std::endl;
            }
         }
      }
   }

   /*
    * If above cut algorithms fail to break or improve the penalty, try
    * more cutting algorithms.
    */
   {

      std::vector<hier::Box> cubic1_breakoff;
      std::vector<hier::Box> cubic1_leftover;
      double cubic1_brk_load;

      bool found_this_break = breakOffLoad_cubic1(
            mapped_box,
            ideal_load_to_break,
            bad_cuts,
            cubic1_breakoff,
            cubic1_leftover,
            cubic1_brk_load);

      if (found_this_break) {
         found_any_break = true;
         double cubic1_balance_penalty = computeBalancePenalty(
               cubic1_breakoff,
               cubic1_leftover,
               (double)(cubic1_brk_load - ideal_load_to_break));
         double cubic1_surface_penalty = computeSurfacePenalty(cubic1_breakoff,
               cubic1_leftover);
         double cubic1_slender_penalty = computeSlenderPenalty(cubic1_breakoff,
               cubic1_leftover);
         double cubic1_combined_penalty =
            combinedBreakingPenalty(cubic1_balance_penalty,
               cubic1_surface_penalty,
               cubic1_slender_penalty);
         if (d_print_break_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "      Cubic-break broke off "
                       << cubic1_brk_load << " / " << ideal_load_to_break
                       << " from " << mapped_box << '|'
                       << mapped_box.numberCells() << '|'
                       << mapped_box.size() << " into "
                       << cubic1_breakoff.size()
                       << " breakoff: ";
            for (std::vector<hier::Box>::const_iterator bi =
                    cubic1_breakoff.begin();
                 bi != cubic1_breakoff.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        and " << cubic1_leftover.size()
                       << " leftover boxes:";
            for (std::vector<hier::Box>::const_iterator bi =
                    cubic1_leftover.begin();
                 bi != cubic1_leftover.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        imbalance: "
                       << (cubic1_brk_load - ideal_load_to_break)
                       << " balance,surface,slender,combined penalties: "
                       << cubic1_balance_penalty << ", "
                       << cubic1_surface_penalty << ", "
                       << cubic1_slender_penalty
                       << ", " << cubic1_combined_penalty << std::endl;
         }
         if (cubic1_combined_penalty < best_combined_penalty) {
            if (d_print_break_steps) {
               tbox::plog << "      choosing breakOffLoad_cubic1 result."
                          << std::endl;
            }
            breakoff.swap(cubic1_breakoff);
            leftover.swap(cubic1_leftover);
            brk_load = cubic1_brk_load;
            best_balance_penalty = cubic1_balance_penalty;
            best_surface_penalty = cubic1_surface_penalty;
            best_slender_penalty = cubic1_slender_penalty;
            best_combined_penalty = cubic1_combined_penalty;
         } else {
            if (d_print_break_steps) {
               tbox::plog << "      Rejecting cubic1 cut result." << std::endl;
            }
         }
      } else {
         if (d_print_break_steps) {
            tbox::plog << "      breakOffLoad_cubic1 could not break "
                       << ideal_load_to_break << " from " << mapped_box
                       << '/' << mapped_box.numberCells()
                       << '/' << mapped_box.numberCells().getProduct()
                       << std::endl;
         }
      }

   }

   t_break_off_load->stop();

   return found_any_break;
}



/*
 *************************************************************************
 * Measuring surface area of boxes is used in penalizing
 * the creation of new surfaces.
 *************************************************************************
 */

double TreeLoadBalancer::computeBoxSurfaceArea(
   const std::vector<hier::Box>& boxes) const
{
   int surface_area = 0;
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        ++bi) {
      const hier::Box& box = *bi;
      int volume = box.size();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         surface_area += volume / box.numberCells(d);
      }
   }
   surface_area *= 2;
   return surface_area;
}



/*
 *************************************************************************
 * Measuring surface area of boxes is used in penalizing
 * the creation of new surfaces.
 *************************************************************************
 */

int TreeLoadBalancer::computeBoxSurfaceArea(
   const hier::Box& box) const
{
   int surface_area = 0;
   int volume = box.size();
   for (int d = 0; d < d_dim.getValue(); ++d) {
      surface_area += volume / box.numberCells(d);
   }
   surface_area *= 2;
   return surface_area;
}



/*
 *************************************************************************
 * Compute the total combined penalty associated with box cutting.
 * Any kind of binary function can be used.  I am just experimenting
 * with various types and weights.
 *
 * All input penalties should all be non-dimensional.
 *************************************************************************
 */

double TreeLoadBalancer::combinedBreakingPenalty(
   double balance_penalty,
   double surface_penalty,
   double slender_penalty) const
{
   double combined_penalty =
      d_balance_penalty_wt * balance_penalty * balance_penalty
      + d_surface_penalty_wt * surface_penalty * surface_penalty
      + d_slender_penalty_wt * slender_penalty * slender_penalty;
   return combined_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted balance penalty for
 * two containers of boxes and how imbalanced they are.
 *************************************************************************
 */

double TreeLoadBalancer::computeBalancePenalty(
   const std::vector<hier::Box>& a,
   const std::vector<hier::Box>& b,
   double imbalance) const
{
   NULL_USE(a);
   NULL_USE(b);
   return tbox::MathUtilities<double>::Abs(imbalance);
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted balance penalty for
 * two containers of boxes and how imbalanced they are.
 *************************************************************************
 */

double TreeLoadBalancer::computeBalancePenalty(
   const TransitSet& a,
   const TransitSet& b,
   double imbalance) const
{
   NULL_USE(a);
   NULL_USE(b);
   return tbox::MathUtilities<double>::Abs(imbalance);
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted balance penalty for
 * a box and how imbalanced it is.
 *************************************************************************
 */

double TreeLoadBalancer::computeBalancePenalty(
   const hier::Box& a,
   double imbalance) const
{
   NULL_USE(a);
   return tbox::MathUtilities<double>::Abs(imbalance);
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double TreeLoadBalancer::computeSurfacePenalty(
   const std::vector<hier::Box>& a,
   const std::vector<hier::Box>& b) const
{
   double surface_penalty = 0;
   for (std::vector<hier::Box>::const_iterator bi = a.begin();
        bi != a.end();
        ++bi) {
      surface_penalty += computeSurfacePenalty(*bi);
   }
   for (std::vector<hier::Box>::const_iterator bi = b.begin();
        bi != b.end();
        ++bi) {
      surface_penalty += computeSurfacePenalty(*bi);
   }
   return surface_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double TreeLoadBalancer::computeSurfacePenalty(
   const TransitSet& a,
   const TransitSet& b) const
{
   double surface_penalty = 0;
   for (TransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->mapped_box);
   }
   for (TransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->mapped_box);
   }
   return surface_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double TreeLoadBalancer::computeSurfacePenalty(
   const hier::Box& a) const
{
   int boxvol = a.size();
   double surface_area = computeBoxSurfaceArea(a);
   double best_surface = 2 * d_dim.getValue() * pow((double)boxvol,
         (double)(d_dim.getValue() - 1) / d_dim.getValue());
   double surface_penalty = surface_area / best_surface - 1.0;
   surface_penalty = surface_penalty * surface_penalty; // Make it blow up.
   return surface_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted slenderness penalty for two
 * containers of boxes.  The reference zero penalty refers to a box with
 * all sides the same length.
 *************************************************************************
 */

double TreeLoadBalancer::computeSlenderPenalty(
   const std::vector<hier::Box>& a,
   const std::vector<hier::Box>& b) const
{
   double slender_penalty = 0;
   for (std::vector<hier::Box>::const_iterator bi = a.begin();
        bi != a.end();
        ++bi) {
      slender_penalty += computeSlenderPenalty(*bi);
   }
   for (std::vector<hier::Box>::const_iterator bi = b.begin();
        bi != b.end();
        ++bi) {
      slender_penalty += computeSlenderPenalty(*bi);
   }
   return slender_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted slenderness penalty for two
 * containers of boxes.  The reference zero penalty refers to a box with
 * all sides the same length.
 *************************************************************************
 */

double TreeLoadBalancer::computeSlenderPenalty(
   const TransitSet& a,
   const TransitSet& b) const
{
   double slender_penalty = 0;
   for (TransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->mapped_box);
   }
   for (TransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->mapped_box);
   }
   return slender_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted slenderness penalty for two
 * containers of boxes.  The reference zero penalty refers to a box with
 * all sides the same length.
 *************************************************************************
 */

double TreeLoadBalancer::computeSlenderPenalty(
   const hier::Box& a) const
{
   const hier::IntVector boxdim = a.numberCells();
   double slender_penalty = (double)boxdim.max() / boxdim.min() - 1.0;
   slender_penalty = slender_penalty * slender_penalty; // Make it blow up.
   return slender_penalty;
}



/*
 *************************************************************************
 * Break up box bursty against box solid and adds the pieces to list.
 * This version differs from that in BoxList in that it tries to minimize
 * slivers.
 *************************************************************************
 */

void TreeLoadBalancer::burstBox(
   const hier::Box& bursty,
   const hier::Box& solid,
   std::vector<hier::Box>& boxes) const
{
   /*
    * This method lacks logic to handle the case of solid not being
    * completely inside bursty.  That feature is not currently needed.
    */
   TBOX_ASSERT(bursty.contains(solid));

   const hier::IntVector solid_size = solid.numberCells();

   boxes.clear();
   hier::Box cutme = bursty;
   while (!cutme.isSpatiallyEqual(solid)) {

      int cut_dir = 999999;
      bool cut_above_solid = false; // Whether to slice off the piece above solid (vs below).
      /*
       * Find direction and place to cut.  To minimize slivers, cut
       * from cutme the thickest slab (in direction normal to cut)
       * possible.
       */
      int slab_thickness = 0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (cutme.numberCells(d) > solid_size(d)) {
            const int thickness_from_upper_cut = cutme.upper() (d)
               - solid.upper() (d);
            if (thickness_from_upper_cut > slab_thickness) {
               slab_thickness = thickness_from_upper_cut;
               cut_dir = d;
               cut_above_solid = true;
            }
            const int thickness_from_lower_cut = solid.lower() (d)
               - cutme.lower() (d);
            if (thickness_from_lower_cut > slab_thickness) {
               slab_thickness = thickness_from_lower_cut;
               cut_dir = d;
               cut_above_solid = false;
            }
         }
      }
      TBOX_ASSERT(cut_dir >= 0 && cut_dir < d_dim.getValue());

      hier::Box removeme = cutme;
      if (cut_above_solid) {
         cutme.upper() (cut_dir) = solid.upper() (cut_dir);
         removeme.lower() (cut_dir) = solid.upper() (cut_dir) + 1;
      } else {
         cutme.lower() (cut_dir) = solid.lower() (cut_dir);
         removeme.upper() (cut_dir) = solid.lower() (cut_dir) - 1;
      }

      boxes.push_back(removeme);

   }

   if (d_print_break_steps) {
      tbox::plog << "      burstBox: " << bursty << " = " << solid;
      for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
           bi != boxes.end();
           bi++) {
         tbox::plog << " + " << *bi;
      }
      tbox::plog << std::endl;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        bi++) {
      for (std::vector<hier::Box>::const_iterator bj = boxes.begin();
           bj != boxes.end();
           bj++) {
         if (bi != bj) {
            TBOX_ASSERT(!bi->intersects(*bj));
         }
      }
   }
   hier::BoxList l1(bursty);
   hier::BoxList l2(solid);
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        bi++) {
      l2.pushFront(*bi);
   }
   l1.removeIntersections(l2);
   TBOX_ASSERT(l1.size() == 0);
   l2.removeIntersections(bursty);
   TBOX_ASSERT(l2.size() == 0);
#endif
}



/*
 *************************************************************************
 *************************************************************************
 */
double TreeLoadBalancer::computeLocalLoads(
   const hier::BoxLevel& mapped_box_level) const
{
   // Count up workload.
   double load = 0.0;
   const hier::BoxSet& mapped_boxes = mapped_box_level.getBoxes();
   for (hier::BoxSet::OrderedConstIterator ni = mapped_boxes.orderedBegin();
        ni != mapped_boxes.orderedEnd();
        ++ni) {
      double mapped_box_load = computeLoad(*ni);
      load += mapped_box_load;
   }
   return (double)load;
}



/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void TreeLoadBalancer::printClassData(
   std::ostream& os) const
{
   os << "\nTreeLoadBalancer::printClassData..." << std::endl;
   os << "\nTreeLoadBalancer: this = "
      << (TreeLoadBalancer *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;

   int ln;

   os << "d_workload_data_id..." << std::endl;
   for (ln = 0; ln < d_workload_data_id.getSize(); ln++) {
      os << "    d_workload_data_id[" << ln << "] = "
         << d_workload_data_id[ln] << std::endl;
   }
}



/*
 *************************************************************************
 *
 * Read values (described in the class header) from input database.
 *
 *************************************************************************
 */

void TreeLoadBalancer::getFromInput(
   tbox::Pointer<tbox::Database> db)
{

   if (!db.isNull()) {

      d_print_steps =
         db->getBoolWithDefault("print_steps",
            d_print_steps);
      d_print_break_steps =
         db->getBoolWithDefault("print_break_steps",
            d_print_break_steps);
      d_print_swap_steps =
         db->getBoolWithDefault("print_swap_steps",
            d_print_swap_steps);
      d_print_edge_steps =
         db->getBoolWithDefault("print_edge_steps",
            d_print_edge_steps);
      d_check_connectivity =
         db->getBoolWithDefault("check_connectivity",
            d_check_connectivity);
      d_check_map =
         db->getBoolWithDefault("check_map",
            d_check_map);

      d_report_load_balance = db->getBoolWithDefault("report_load_balance",
            d_report_load_balance);
      d_barrier_before = db->getBoolWithDefault("barrier_before",
            d_barrier_before);
      d_barrier_after = db->getBoolWithDefault("barrier_after",
            d_barrier_after);

      d_n_root_cycles = db->getIntegerWithDefault("n_root_cycles",
            d_n_root_cycles);

      d_balance_penalty_wt = db->getDoubleWithDefault("balance_penalty_wt",
            d_balance_penalty_wt);
      d_surface_penalty_wt = db->getDoubleWithDefault("surface_penalty_wt",
            d_surface_penalty_wt);
      d_slender_penalty_wt = db->getDoubleWithDefault("slender_penalty_wt",
            d_slender_penalty_wt);
      d_precut_penalty_wt = db->getDoubleWithDefault("precut_penalty_wt",
            d_precut_penalty_wt);

   }
}



/*
 *************************************************************************
 *
 * Private utility functions to determine parameter values for level.
 *
 *************************************************************************
 */

int TreeLoadBalancer::getWorkloadDataId(
   int level_number) const
{

   TBOX_ASSERT(level_number >= 0);

   int wrk_id = (level_number < d_workload_data_id.getSize() ?
                 d_workload_data_id[level_number] :
                 d_master_workload_data_id);

   return wrk_id;
}



/*
 ***************************************************************************
 *
 ***************************************************************************
 */
void TreeLoadBalancer::assertNoMessageForPrivateCommunicator() const
{
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if (d_mpi_dup.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      int flag;
      tbox::SAMRAI_MPI::Status mpi_status;
      int mpi_err = d_mpi_dup.Iprobe(MPI_ANY_SOURCE,
            MPI_ANY_TAG,
            &flag,
            &mpi_status);
      if (mpi_err != MPI_SUCCESS) {
         TBOX_ERROR("Error probing for possible lost messages.");
      }
      if (flag == true) {
         int count = -1;
         mpi_err = tbox::SAMRAI_MPI::Get_count(&mpi_status, MPI_INT, &count);
         TBOX_ERROR(
            "Library error!\n"
            << "TreeLoadBalancer detected before or\n"
            << "after using a private communicator that there\n"
            << "is a message yet to be received.  This is\n"
            << "an error because all messages using the\n"
            << "private communicator should have been\n"
            << "accounted for.  Message status:\n"
            << "source " << mpi_status.MPI_SOURCE << '\n'
            << "tag " << mpi_status.MPI_TAG << '\n'
            << "count " << count << " (assuming integers)\n"
            << "current tags: "
            << ' ' << TreeLoadBalancer_LOADTAG0 << ' '
            << TreeLoadBalancer_LOADTAG1
            );
      }
   }
}



/*
 *************************************************************************
 * Break off a load from a box by making a single planar cut across
 * the box's longest dimension.
 *
 * Currently assuming uniform load of one unit per cell.
 *
 * Return whether a successful break was made.
 *************************************************************************
 */
bool TreeLoadBalancer::breakOffLoad_planar(
   const hier::Box& mapped_box,
   double ideal_load_to_break,
   const tbox::Array<tbox::Array<bool> >& bad_cuts,
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load) const
{

   const tbox::Dimension dim(d_dim);

   if (d_print_break_steps) {
      tbox::plog << "      breakOffLoad_planar attempting to break "
                 << ideal_load_to_break << " from Box " << mapped_box
                 << " min_size=" << d_min_size << std::endl;
   }

   brk_load = 0;
   breakoff.clear();
   leftover.clear();

   const hier::IntVector& box_dims = mapped_box.numberCells();

   const int box_vol = box_dims.getProduct();

   if (box_vol <= ideal_load_to_break) {
      // Easy: break off everything.
      breakoff.push_back(mapped_box);
      brk_load = box_vol;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_planar broke off entire Box "
                    << mapped_box
                    << std::endl;
      }
      return true;
   }

   /*
    * Determine ordering of box_dims from shortest to longest.
    */
   hier::IntVector ordered_dims(dim);
   sortIntVector(ordered_dims, box_dims);

   /*
    * best_difference is the difference between the best cut found and
    * ideal_load_to_break.  Initialize it for zero breakoff.
    */
   double best_difference = (double)ideal_load_to_break;
   hier::Box best_breakoff_box(dim);
   hier::Box best_leftover_box(dim);

   for (int d = d_dim.getValue() - 1; d >= 0; --d) {

      /*
       * Search directions from longest to shortest
       * because we prefer to break across longest dir.
       */
      const int brk_dir = ordered_dims(d_dim.getValue() - 1);

      const int brk_area = box_vol / box_dims(brk_dir);

      const tbox::Array<bool>& bad = bad_cuts[brk_dir];

      /*
       * Try rounding the break length down (round==0) and up (round==1),
       * looking for the best place to cut.
       */
      for (int round = 0; round <= 1; ++round) {

         /*
          * Rounding up or down must heed d_cut_factor,
          * so brk_len must be an integer multiple of d_cut_factor.
          */
         const int brk_len =
            (((int)(ideal_load_to_break / brk_area))
             / d_cut_factor(brk_dir) + round)
            * d_cut_factor(brk_dir);

         if (brk_len < d_min_size(brk_dir)) {
            // Breakoff box violates minimum size.
            continue;
         }

         if (box_dims(brk_dir) - brk_len > 0 &&
             box_dims(brk_dir) - brk_len < d_min_size(brk_dir)) {
            // Leftover box violates minimum size.
            continue;
         }

         const int brk_volume = brk_area * brk_len;

         /*
          * Compute the difference between the current breakage
          * and the ideal.
          */
         const double difference = brk_volume <= ideal_load_to_break ?
            ((double)ideal_load_to_break - brk_volume) :
            ((double)brk_volume - ideal_load_to_break);

         if (difference < best_difference) {
            // This cut gives better difference, if it can be done.

#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(brk_len >= 0 && brk_len <= bad.size());
#endif
            if (brk_len == box_dims(brk_dir) ||
                bad[brk_len] == false) {
               // Cutting brk_len from low side is ok.
               best_difference = difference;
               best_breakoff_box = mapped_box;
               best_breakoff_box.upper() (brk_dir) =
                  best_breakoff_box.lower() (brk_dir) + brk_len - 1;
               best_leftover_box = mapped_box;
               best_leftover_box.lower() (brk_dir) =
                  best_breakoff_box.upper() (brk_dir) + 1;
               break;
            } else if (bad[box_dims(brk_dir) - brk_len] == false) {
               // Cutting brk_len from high side is ok.
               best_difference = difference;
               best_breakoff_box = mapped_box;
               best_breakoff_box.lower() (brk_dir) =
                  best_breakoff_box.upper() (brk_dir) - brk_len + 1;
               best_leftover_box = mapped_box;
               best_leftover_box.upper() (brk_dir) =
                  best_breakoff_box.lower() (brk_dir) - 1;
               break;
            }
         }
      }

   }

   bool successful_break = false;

   if (!best_breakoff_box.empty()) {
      breakoff.push_back(best_breakoff_box);
      brk_load = best_breakoff_box.size();
      successful_break = true;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_planar broke off box " << mapped_box
                    << " for breakoff box " << best_breakoff_box
                    << " and leftover " << best_leftover_box << std::endl;
      }
   } else {
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_planar could not break "
                    << ideal_load_to_break << " from Box " << mapped_box
                    << std::endl;
      }
   }
   if (!best_leftover_box.empty()) {
      leftover.push_back(best_leftover_box);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   for (std::vector<hier::Box>::iterator bi = breakoff.begin();
        bi != breakoff.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancer library error:\n"
               << "breakoff box " << b << ", size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "break box size " << best_breakoff_box.numberCells() << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
   for (std::vector<hier::Box>::iterator bi = leftover.begin();
        bi != leftover.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancer library error:\n"
               << "leftover box " << b << ", size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "break box size " << best_breakoff_box.numberCells() << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
#endif

   return successful_break;
}



/*
 *************************************************************************
 * Break off a load from a box by making a box cut that is as close
 * to cubic as feasible.
 *
 * Currently assuming uniform load of one unit per cell.
 *
 * Return whether a successful break was made.
 *
 * Differs from breakOffLoad in that it will always
 * performs a break and if needed, break off more than
 * the ideal.  The calling method should take this into
 * account.
 *************************************************************************
 */
bool TreeLoadBalancer::breakOffLoad_cubic1(
   const hier::Box& mapped_box,
   double ideal_load_to_break,
   const tbox::Array<tbox::Array<bool> >& bad_cuts,
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load) const
{

   const hier::IntVector box_dims(mapped_box.numberCells());

   const double mapped_box_load(box_dims.getProduct());

   if (ideal_load_to_break >= mapped_box_load) {
      // Easy: break off everything.
      leftover.clear();
      breakoff.clear();
      breakoff.push_back(mapped_box);
      brk_load = mapped_box_load;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_cubic1 broke off entire Box "
                    << mapped_box
                    << std::endl;
      }
      return true;
   }

   if (ideal_load_to_break > 0.5 * mapped_box_load) {
      /*
       * This algorithm is better when breaking off a small portion.
       * Since the ideal is a bigger portion, switch breakoff with leftover.
       */
      if (d_print_break_steps) {
         tbox::plog
         << "      breakOffload_cubic1 reversing direction to break "
         << (box_dims.getProduct() - ideal_load_to_break)
         << " instead of " << ideal_load_to_break << " / "
         << box_dims.getProduct() << std::endl;
      }
      bool success =
         breakOffLoad_cubic1(mapped_box,
            box_dims.getProduct() - ideal_load_to_break,
            bad_cuts,
            leftover,
            breakoff,
            brk_load);
      if (success) {
         brk_load = box_dims.getProduct() - brk_load;
      }
      return success;
   }

   if (d_print_break_steps) {
      tbox::plog << "      breakOffload_cubic1 attempting to break "
                 << ideal_load_to_break << " from Box " << mapped_box
                 << " min_size=" << d_min_size << std::endl;
   }

   breakoff.clear();
   leftover.clear();

   /*
    * brk_size is the size of the box we want to break off of
    * mapped_box.  We start with the smallest allowed brk_size that
    * will not create remainders that violate size constraints.
    *
    * In the do loop below, we increase brk_size to bring brk_load
    * closer to ideal_oad_to_break.
    */
   hier::IntVector brk_size(d_min_size);
   brk_size.max(d_cut_factor);
   brk_size.min(box_dims);
   /*
    * If remainder is too small, zero it out to avoid
    * having non-zero remainder smaller than d_min_size.
    */
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if ((box_dims(d) - brk_size(d) > 0) &&
          (box_dims(d) - brk_size(d) < d_min_size(d))) {
         brk_size(d) = box_dims(d);
      }
   }
   brk_load = brk_size.getProduct();

   if (d_print_break_steps) {
      tbox::plog << "      brk: " << std::flush;
      tbox::plog << "  " << brk_size << "->" << brk_load << std::flush;
   }

   /*
    * stop_growing: whether dimensions of brk_size is already
    * big engough so that it cannot not be grown without breaking
    * off too much.
    */
   hier::IntVector stop_growing(d_dim, 0);
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (brk_size[d] == box_dims[d]) stop_growing[d] = 1;
   }

   if (brk_load < ideal_load_to_break) {
      /*
       * The do loop gradually increases brk_size to bring brk_load
       * closer to ideal_load_to_break.
       *
       * Select dimension to increase in size, inc_dim.  Use the
       * smallest dimension that is still allowed to grow.
       *
       * Try a new break size that is bigger than the current brk_size
       * by the minimal allowed amount.  If this brings us closer to
       * the ideal break amount, mark it.
       *
       * Exit the loop when we cannot grow anymore or we are already
       * breaking off more than the ideal.
       */
      do {

         int inc_dim = -1;
         for (int d = 0; d < d_dim.getValue(); ++d) {
            if (!stop_growing(d) &&
                (inc_dim == -1 || brk_size(inc_dim) > brk_size(d))) inc_dim = d;
         }
         if (inc_dim == -1) break; // No growable dimension.

         TBOX_ASSERT(brk_size(inc_dim) < box_dims(inc_dim));

         hier::IntVector new_brk_size(brk_size);
         new_brk_size(inc_dim) += d_cut_factor(inc_dim);
         if (d_print_break_steps) {
            tbox::plog << "  " << brk_size << "=>" << brk_load << std::flush;
         }

         // Prevent remainder being smaller than d_min_size.
         if (box_dims(inc_dim) - new_brk_size(inc_dim) < d_min_size(inc_dim)) {
            new_brk_size(inc_dim) = box_dims(inc_dim);
            if (d_print_break_steps) {
               tbox::plog << "  " << brk_size << "==>" << brk_load
                          << std::flush;
            }
         }

         if (new_brk_size(inc_dim) == box_dims(inc_dim)) {
            stop_growing(inc_dim) = 1;
         }

         int new_brk_load = new_brk_size.getProduct();

         if (new_brk_load <= ideal_load_to_break) {
            /*
             * new_brk_load is closer to ideal than current brk_load.
             * Don't break out of the loop yet.  We will grow it again
             * and check.
             */
            brk_size = new_brk_size;
            brk_load = new_brk_load;
            if (d_print_break_steps) {
               tbox::plog << "  " << brk_size << "===>" << brk_load
                          << std::flush;
            }
         } else if ((new_brk_load - ideal_load_to_break) <
                    ((double)ideal_load_to_break - brk_load)) {
            /*
             * new_brk_load is bigger than ideal but is still an
             * improvement over brk_load.  Accept it, but break out of
             * the loop because further growing will not improve
             * result.
             */
            brk_size = new_brk_size;
            brk_load = new_brk_load;
            if (d_print_break_steps) {
               tbox::plog << "  " << brk_size << "====>" << brk_load
                          << std::flush;
            }
            break;
         } else {
            /*
             * Even though dimension inc_dim has not reached the box
             * dimension, stop growing it because any further growth
             * leads to too big a load.
             */
            stop_growing(inc_dim) = 1;
         }

      } while (true);
   }

   if (d_print_break_steps) {
      tbox::plog << std::endl;
   }

   /*
    * Find a place to put the break-off box so that it does not lie
    * across a bad cut.  If no such placement is found, set
    * placement_impossible to true.
    */
   hier::Box breakoff_box(d_dim);
   const hier::IntVector& lower(mapped_box.lower());
   const hier::IntVector& upper(mapped_box.upper());
   bool placement_impossible = false;
   if (d_print_break_steps) {
      tbox::plog << "      Placing " << brk_size
                 << " to avoid bad cut points" << std::flush;
   }
   for (int d = 0; d < d_dim.getValue(); ++d) {
      /*
       * To minimize the number of boxes remaining after breakoff_box
       * is removed, prefer to place breakoff_box against the upper or
       * lower side of the Box.  First try putting the breakoff_box
       * along the upper side of dimension d.  If it cuts the box at a
       * bad point, try the lower side.  If it still cuts at a bad
       * points, slide the box toward the upper side until it does not
       * cut at any bad points, being careful not to be so close to
       * the box boundaries that we generate remainder boxes violating
       * the minimum size.  If no place can be found to put
       * breakoff_box, set placement_impossible and give up.  (We
       * could go back and reshape the box at this point, but we won't
       * because there is probably another box that would work without
       * reshaping.)
       */
      const tbox::Array<bool>& bad = bad_cuts[d];

      int& gl = breakoff_box.lower()[d];
      int& gu = breakoff_box.upper()[d];

      gu = upper[d];
      gl = gu - (brk_size[d] - 1);
      if (!bad[gl - lower[d]]) {
         if (d_print_break_steps) {
            tbox::plog << "  d=" << d << " upper" << std::flush;
         }
         continue;
      }

      gl = lower[d];
      gu = gl + (brk_size[d] - 1);
      if (gu + 1 - lower[d] < bad.size() && !bad[gu + 1 - lower[d]]) {
         if (d_print_break_steps) {
            tbox::plog << "  d=" << d << " lower" << std::flush;
         }
         continue;
      }

      gl = lower[d] + d_min_size[d];
      gu = gl + (brk_size[d] - 1);
      while (gu <= upper[d] - d_min_size[d] &&
             (bad[gl - lower[d]] || bad[gu + 1 - lower[d]])) {
         ++gl;
         ++gu;
      }
      if (gu <= upper[d] - d_min_size[d]) {
         if (d_print_break_steps) {
            tbox::plog << "  d=" << d << " middle" << std::flush;
         }
         continue;
      }

      if (d_print_break_steps) {
         tbox::plog << "  cannot place dim " << d
                    << " without creating bad cuts."
                    << std::flush;
      }
      placement_impossible = true;
      /*
       * Cannot find place for breakoff_box along dimension d.
       * No point in looking at higher dimensions.
       */
      break;
   }
   if (d_print_break_steps) {
      tbox::plog << std::endl;
   }
   if (placement_impossible) {
      return false;
   }

   breakoff.clear();
   breakoff.push_back(breakoff_box);

   burstBox(mapped_box,
      breakoff_box,
      leftover);

#ifdef DEBUG_CHECK_ASSERTIONS
   for (std::vector<hier::Box>::iterator bi = breakoff.begin();
        bi != breakoff.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancer library error:\n"
               << "breakoff box " << b << ", with size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "orig box " << mapped_box << "\n"
               << "break box " << breakoff_box << "\n"
               << "break box size " << brk_size << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
   for (std::vector<hier::Box>::iterator bi = leftover.begin();
        bi != leftover.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancer library error:\n"
               << "leftover box " << b << ", with size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "orig box " << mapped_box << "\n"
               << "break box " << breakoff_box << "\n"
               << "break box size " << brk_size << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
#endif

   return true;
}



/*
**************************************************************************
**************************************************************************
*/

void TreeLoadBalancer::setShadowData(
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_mapped_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const hier::IntVector& refinement_ratio) const
{
   d_min_size = min_size;
   d_max_size = max_size;
   d_bad_interval = bad_interval;
   d_cut_factor = cut_factor;
   /*
    * Domain boxes are used by breakOffLoad to determine where
    * the bad cuts are.  Computing domain_boxes from domain_mapped_box_level
    * should be moved above the this method.
    */

   /*
    * We expect the domain mapped_box_level to be in globalized state.
    */
   TBOX_ASSERT(
      domain_mapped_box_level.getParallelState() ==
      hier::BoxLevel::GLOBALIZED);

   d_domain_boxes.clear();
   domain_mapped_box_level.getGlobalBoxes(d_domain_boxes);
   d_domain_boxes.refine(refinement_ratio);
}



/*
**************************************************************************
**************************************************************************
*/

void TreeLoadBalancer::unsetShadowData() const {
   d_min_size = hier::IntVector(d_dim, -1);
   d_max_size = hier::IntVector(d_dim, -1);
   d_domain_boxes.clear();
   d_bad_interval = hier::IntVector(d_dim, -1);
   d_cut_factor = hier::IntVector(d_dim, -1);
}



/*
**************************************************************************
* Move Boxes in balance_mapped_box_level from ranks outside of
* rank_group to ranks inside rank_group.  Modify the given connectors
* to make them correct following this moving of boxes.
**************************************************************************
*/

void TreeLoadBalancer::prebalanceBoxLevel(
   hier::BoxLevel& balance_mapped_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const tbox::RankGroup& rank_group) const
{
   hier::OverlapConnectorAlgorithm oca;

   if (balance_to_anchor.isInitialized()) {
      TBOX_ASSERT(anchor_to_balance.checkTransposeCorrectness(balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor.checkTransposeCorrectness(anchor_to_balance) == 0);
   }

   /*
    * tmp_mapped_box_level will contain the same boxes as
    * balance_mapped_box_level, but all will live on the processors
    * specified in rank_group.
    */
   hier::BoxLevel tmp_mapped_box_level(d_dim);
   tmp_mapped_box_level.initialize(
      balance_mapped_box_level.getRefinementRatio(),
      balance_mapped_box_level.getGridGeometry(),
      balance_mapped_box_level.getMPI());

   /*
    * If a rank is not in rank_group it is called a "sending" rank, as
    * it will send any Boxes it has to a rank in rank_group.
    */
   bool is_sending_rank = rank_group.isMember(d_mpi.getRank()) ? false : true;

   int output_nproc = rank_group.size();

   /*
    * the send and receive comm objects
    */
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* box_send = NULL;
   tbox::AsyncCommPeer<int>* box_recv = NULL;
   tbox::AsyncCommPeer<int>* id_send = NULL;
   tbox::AsyncCommPeer<int>* id_recv = NULL;

   /*
    * A sending rank will send its Boxes to a receiving rank, and
    * that receiving processor will add it to its local set of Boxes.
    * When the box is added on the receiving processor, it will receive
    * a new LocalId.  This LocalId value needs to be sent back to
    * the sending processor, in order to construct the mapping connectors.
    *
    * Therefore the sending ranks construct comm objects for sending boxes
    * and receiving LocalIdes.
    *
    * Sending processors send to ranks in the rank_group determined by
    * a modulo heuristic.
    */
   if (is_sending_rank) {
      box_send = new tbox::AsyncCommPeer<int>;
      box_send->initialize(&comm_stage);
      box_send->setPeerRank(rank_group.getMappedRank(d_mpi.getRank() % output_nproc));
      box_send->setMPI(d_mpi);
      box_send->setMPITag(TreeLoadBalancer_PREBALANCE0 + 2 * d_mpi.getRank(),
         TreeLoadBalancer_PREBALANCE1 + 2 * d_mpi.getRank());

      id_recv = new tbox::AsyncCommPeer<int>;
      id_recv->initialize(&comm_stage);
      id_recv->setPeerRank(rank_group.getMappedRank(d_mpi.getRank() % output_nproc));
      id_recv->setMPI(d_mpi);
      id_recv->setMPITag(TreeLoadBalancer_PREBALANCE0 + 2 * d_mpi.getRank(),
         TreeLoadBalancer_PREBALANCE1 + 2 * d_mpi.getRank());
   }

   /*
    * The receiving ranks construct comm objects for receiving boxes
    * and sending LocalIdes.
    */
   int num_recvs = 0;
   if (rank_group.isMember(d_mpi.getRank())) {
      tbox::List<int> recv_ranks;
      for (int i = 0; i < d_mpi.getSize(); i++) {
         if (!rank_group.isMember(i) &&
             rank_group.getMappedRank(i % output_nproc) == d_mpi.getRank()) {
            recv_ranks.appendItem(i);
         }
      }
      num_recvs = recv_ranks.size();
      if (num_recvs > 0) {
         box_recv = new tbox::AsyncCommPeer<int>[num_recvs];
         id_send = new tbox::AsyncCommPeer<int>[num_recvs];
         int recv_count = 0;
         for (tbox::List<int>::Iterator ri(recv_ranks); ri; ri++) {
            box_recv[recv_count].initialize(&comm_stage);
            box_recv[recv_count].setPeerRank(ri());
            box_recv[recv_count].setMPI(d_mpi);
            box_recv[recv_count].setMPITag(TreeLoadBalancer_PREBALANCE0 + 2 * ri(),
               TreeLoadBalancer_PREBALANCE1 + 2 * ri());

            id_send[recv_count].initialize(&comm_stage);
            id_send[recv_count].setPeerRank(ri());
            id_send[recv_count].setMPI(d_mpi);
            id_send[recv_count].setMPITag(TreeLoadBalancer_PREBALANCE0 + 2 * ri(),
               TreeLoadBalancer_PREBALANCE1 + 2 * ri());

            recv_count++;
         }
         TBOX_ASSERT(num_recvs == recv_count);
      }
   }

   /*
    * These NeighborhoodSets are used to create mapping Connectors which
    * describe the mapping from the box configuration of the given
    * balance_mapped_box_level, to the new configuration stored in
    * tmp_mapped_box_level.  These mapping Connectors are necessary to
    * modify the two Connectors given in the argument list, so that on
    * return from this method, they will be correct for the new
    * balance_mapped_box_level.
    */
   hier::NeighborhoodSet balance_to_tmp_edges(d_dim);
   hier::NeighborhoodSet tmp_to_balance_edges(d_dim);

   /*
    * Where Boxes already exist on ranks in rank_group,
    * move them directly to tmp_mapped_box_level.
    */
   if (!is_sending_rank) {
      const hier::BoxSet& unchanged_mapped_boxes =
         balance_mapped_box_level.getBoxes();

      for (hier::BoxSet::OrderedConstIterator ni =
              unchanged_mapped_boxes.orderedBegin();
           ni != unchanged_mapped_boxes.orderedEnd(); ++ni) {

         const hier::Box& mapped_box = *ni;
         tmp_mapped_box_level.addBox(mapped_box);
      }
   }

   const int buf_size = hier::Box::commBufferSize(d_dim);

   /*
    * On sending ranks, pack the Boxes into buffers and send.
    */
   if (is_sending_rank) {
      const hier::BoxSet& sending_mapped_boxes =
         balance_mapped_box_level.getBoxes();
      const int num_sending_boxes =
         static_cast<int>(sending_mapped_boxes.size());

      int* buffer = new int[buf_size * num_sending_boxes];
      int box_count = 0;
      for (hier::BoxSet::OrderedConstIterator ni = sending_mapped_boxes.orderedBegin();
           ni != sending_mapped_boxes.orderedEnd(); ++ni) {

         const hier::Box& mapped_box = *ni;

         mapped_box.putToIntBuffer(&buffer[box_count * buf_size]);
         box_count++;
      }
      box_send->beginSend(buffer, buf_size * num_sending_boxes);

      delete[] buffer;
   }

   /*
    * On receiving ranks, complete the receives, add the boxes to local
    * tmp_mapped_box_level, insert boxes into tmp_to_balance_edges, and then
    * send the new LocalIdes back to the sending processors.
    */
   if (!is_sending_rank && num_recvs > 0) {
      for (int i = 0; i < num_recvs; i++) {
         box_recv[i].beginRecv();
      }
      int num_completed_recvs = 0;
      tbox::Array<bool> completed(num_recvs, false);
      while (num_completed_recvs < num_recvs) {
         for (int i = 0; i < num_recvs; i++) {
            if (!completed[i] && box_recv[i].checkRecv()) {
               num_completed_recvs++;
               completed[i] = true;
               const int num_boxes = box_recv[i].getRecvSize() / buf_size;
               const int* buffer = box_recv[i].getRecvData();
               int* id_buffer = new int[num_boxes];

               for (int b = 0; b < num_boxes; b++) {
                  hier::Box mapped_box(d_dim);

                  mapped_box.getFromIntBuffer(&buffer[b * buf_size]);

                  hier::BoxSet::OrderedIterator tmp_iter =
                     tmp_mapped_box_level.addBox(mapped_box,
                        mapped_box.getBlockId());

                  hier::BoxId tmp_mapped_box_id = (*tmp_iter).getId();

                  tmp_to_balance_edges.insertNeighbor(tmp_mapped_box_id,
                                                      mapped_box);

                  id_buffer[b] = tmp_mapped_box_id.getLocalId().getValue();
               }
               id_send[i].beginSend(id_buffer, num_boxes);

               delete[] id_buffer;
            }
         }
      }
      for (int i = 0; i < num_recvs; i++) {
         if (!id_send[i].checkSend()) {
            id_send[i].completeCurrentOperation();
         }
      }
   }

   /*
    * On sending ranks, receive the LocalIds, and add the edges
    * to balance_to_tmp_edges.
    */
   if (is_sending_rank) {
      if (!box_send->checkSend()) {
         box_send->completeCurrentOperation();
      }

      id_recv->beginRecv();

      if (!id_recv->checkRecv()) {
         id_recv->completeCurrentOperation();
      }
      const int* buffer = id_recv->getRecvData();

      const hier::BoxSet& sending_mapped_boxes =
         balance_mapped_box_level.getBoxes();
      TBOX_ASSERT(static_cast<unsigned int>(id_recv->getRecvSize()) == sending_mapped_boxes.size());

      int box_count = 0;
      for (hier::BoxSet::OrderedConstIterator ni =
              sending_mapped_boxes.orderedBegin();
           ni != sending_mapped_boxes.orderedEnd(); ++ni) {

         hier::Box new_mapped_box(
            (*ni),
            (hier::LocalId)buffer[box_count],
            rank_group.getMappedRank(d_mpi.getRank() % output_nproc),
            (*ni).getBlockId());

         balance_to_tmp_edges.insertNeighbor((*ni).getId(), new_mapped_box);
         box_count++;
      }
   }

   /*
    * Construct the mapping Connectors with the appropriate NeighborhoodSets.
    */
   hier::Connector balance_to_tmp(balance_mapped_box_level,
                                  tmp_mapped_box_level,
                                  hier::IntVector::getZero(d_dim),
                                  balance_to_tmp_edges);

   hier::Connector tmp_to_balance(tmp_mapped_box_level,
                                  balance_mapped_box_level,
                                  hier::IntVector::getZero(d_dim),
                                  tmp_to_balance_edges);

   if (anchor_to_balance.isInitialized()) {
      /*
       * This modify operation copies tmp_mapped_box_level to
       * balance_mapped_box_level, and changes anchor_to_balance and
       * balance_to_anchor such that they are correct for the new state
       * of balance_mapped_box_level.
       */
      const hier::MappingConnectorAlgorithm mca;
      mca.modify(anchor_to_balance,
         balance_to_anchor,
         balance_to_tmp,
         tmp_to_balance,
         &balance_mapped_box_level,
         &tmp_mapped_box_level);

      TBOX_ASSERT(anchor_to_balance.checkTransposeCorrectness(balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor.checkTransposeCorrectness(anchor_to_balance) == 0);
   } else {
      hier::BoxLevel::swap(balance_mapped_box_level, tmp_mapped_box_level);
   }

   /*
    * Clean up raw pointer allocation.
    */
   if (is_sending_rank) {
      delete box_send;
      delete id_recv;
   }
   if (num_recvs) {
      delete[] box_recv;
      delete[] id_send;
   }
}



/*
**************************************************************************
**************************************************************************
*/

std::ostream& operator << (
   std::ostream& co,
   const BoxInTransit& r)
{
   co << r.mapped_box
   << r.mapped_box.numberCells() << '|'
   << r.mapped_box.numberCells().getProduct();
   for (int i = static_cast<int>(r.proc_hist.size()) - 1; i >= 0; --i) {
      co << '-' << r.proc_hist[i];
   }
   co << '-' << r.orig_mapped_box
   << r.mapped_box.numberCells() << '|'
   << r.mapped_box.numberCells().getProduct();
   return co;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void TreeLoadBalancer::setTimers()
{
   /*
    * The first constructor gets timers from the TimerManager.
    * and sets up their deallocation.
    */
   if (t_load_balance_mapped_box_level.isNull()) {
      t_load_balance_mapped_box_level = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::loadBalanceBoxLevel()");
      t_get_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_map");
      t_use_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::use_map");
      t_constrain_size = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constrain_size");
      t_map_big_boxes = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::mapOversizedBoxes()");
      t_compute_local_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_local_load");
      t_compute_global_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_global_load");
      t_compute_tree_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load");
      t_compute_tree_load0 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load0");
      t_compute_tree_load1 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load1");
      t_compute_tree_load2 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load2");
      t_break_off_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::breakOffLoad()");
      t_find_bad_cuts = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::find_bad_cuts");
      t_reassign_loads = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::reassignLoads()");
      t_shift_loads_by_swapping = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::shiftLoadsBySwapping()");
      t_shift_loads_by_breaking = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::shiftLoadsByBreaking()");
      t_find_swap_pair = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::findLoadSwapPair()");
      t_send_load_to_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_children");
      t_send_load_to_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_parent");
      t_get_load_from_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_children");
      t_get_load_from_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_parent");
      t_send_edge_to_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_edge_to_children");
      t_send_edge_to_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_edge_to_parent");
      t_get_edge_from_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_edge_from_children");
      t_get_edge_from_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_edge_from_parent");
      t_report_loads = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::report_loads");
      t_misc1 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::misc1");
      t_misc2 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::misc2");
      t_misc3 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::misc3");
      t_finish_comms = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::finish_comms");
      t_local_balancing = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::local_balancing");
      t_local_edges = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::local_edges");
      t_pack_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::pack_load");
      t_unpack_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::unpack_load");
      t_pack_edge = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::pack_edge");
      t_unpack_edge = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::unpack_edge");
      t_parent_load_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_load_comm");
      t_children_load_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::children_load_comm");
      t_parent_edge_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_edge_comm");
      t_children_edge_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::children_edge_comm");
      t_barrier_before = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::barrier_before");
      t_barrier_after = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::barrier_after");
      t_MPI_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::MPI_wait");
   }
}



/*
 *************************************************************************
 *************************************************************************
 */
void TreeLoadBalancer::printStatistics(
   std::ostream& s) const
{
   gatherAndReportLoadBalance(d_load_stat,
      tbox::SAMRAI_MPI::getSAMRAIWorld(),
      s);
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
