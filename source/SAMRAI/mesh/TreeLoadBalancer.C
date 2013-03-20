/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer_C
#define included_mesh_TreeLoadBalancer_C

#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/CenteredRankTree.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cmath>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int TreeLoadBalancer::TreeLoadBalancer_LOADTAG0;
const int TreeLoadBalancer::TreeLoadBalancer_LOADTAG1;
const int TreeLoadBalancer::TreeLoadBalancer_EDGETAG0;
const int TreeLoadBalancer::TreeLoadBalancer_EDGETAG1;
const int TreeLoadBalancer::TreeLoadBalancer_PREBALANCE0;
const int TreeLoadBalancer::TreeLoadBalancer_PREBALANCE1;
const int TreeLoadBalancer::TreeLoadBalancer_FIRSTDATALEN;
const int TreeLoadBalancer::TreeLoadBalancer_MIN_NPROC_FOR_AUTOMATIC_MULTICYCLE;

const int TreeLoadBalancer::d_default_data_id = -1;


// Round a to the nearest higher integer divisible by b.  This should work even for a < 0.
#define ROUND_TO_HI(a,b) ((a)-((((a)%(b))-(b))%(b)))
// Round a to the nearest lower integer divisible by b.  This should work even for a < 0.
#define ROUND_TO_LO(a,b) ((a)-((((a)%(b))+(b))%(b)))


/*
 *************************************************************************
 * TreeLoadBalancer constructor.
 *************************************************************************
 */

TreeLoadBalancer::TreeLoadBalancer(
   const tbox::Dimension& dim,
   const std::string& name,
   const boost::shared_ptr<tbox::Database>& input_db,
   const boost::shared_ptr<tbox::RankTreeStrategy> &rank_tree):
   d_dim(dim),
   d_object_name(name),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mpi_is_dupe(false),
   d_max_cycle_spread_ratio(1000000),
   d_rank_tree(rank_tree ? rank_tree : boost::shared_ptr<tbox::RankTreeStrategy>(new tbox::CenteredRankTree) ),
   d_comm_graph_writer(),
   d_master_workload_data_id(d_default_data_id),
   d_flexible_load_tol(0.0),
   d_min_load_fraction_per_box(0.03),
   d_balance_penalty_wt(1.0),
   d_surface_penalty_wt(1.0),
   d_slender_penalty_wt(1.0),
   d_slender_penalty_threshold(3.0),
   d_precut_penalty_wt(1.0),
   // Data shared during balancing.
   d_min_size(d_dim),
   d_max_size(d_dim),
   d_bad_interval(d_dim),
   d_cut_factor(d_dim),
   // Output control.
   d_report_load_balance(false),
   d_summarize_map(false),
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

bool
TreeLoadBalancer::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}



/*
**************************************************************************
**************************************************************************
*/
void
TreeLoadBalancer::setWorkloadPatchDataIndex(
   int data_id,
   int level_number)
{
   boost::shared_ptr<pdat::CellDataFactory<double> > datafact(
      hier::VariableDatabase::getDatabase()->getPatchDescriptor()->
      getPatchDataFactory(data_id),
      BOOST_CAST_TAG);

   TBOX_ASSERT(datafact);

   if (level_number >= 0) {
      int asize = static_cast<int>(d_workload_data_id.size());
      if (asize < level_number + 1) {
         d_workload_data_id.resize(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
            d_workload_data_id[i] = d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = data_id;
      }
   } else {
      d_master_workload_data_id = data_id;
      for (int ln = 0; ln < static_cast<int>(d_workload_data_id.size()); ln++) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}



/*
 *************************************************************************
 * This method implements the abstract LoadBalanceStrategy interface,
 * but it is not where the tree load balancer algorithm is implemented.
 *
 * This method does some preliminary setup then calls
 * loadBalanceWithinRankGroup to compute the new balanced
 * BoxLevel and the mapping Connectors between the old and the new.
 * Then it applies the mapping to update the balance<==>anchor
 * Connectors.  It may do this multiple times, as specified by the
 * cycling parameter.
 *
 * After load balancing, it enforces the maximum size restriction
 * by breaking up large boxes and update balance<==>anchor again.
 *************************************************************************
 */
void
TreeLoadBalancer::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_anchor,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   TBOX_ASSERT(!balance_to_anchor || balance_to_anchor->hasTranspose());
   TBOX_ASSERT(!balance_to_anchor || balance_to_anchor->isTransposeOf(balance_to_anchor->getTranspose()));
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY6(d_dim,
      balance_box_level,
      min_size,
      max_size,
      domain_box_level,
      bad_interval,
      cut_factor);
   if (hierarchy) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);
   }


   if ( d_mpi_is_dupe ) {
      /*
       * If user has set the duplicate communicator, make sure it is
       * compatible with the BoxLevel involved.
       */
      TBOX_ASSERT(d_mpi.getSize() == balance_box_level.getMPI().getSize());
      TBOX_ASSERT(d_mpi.getRank() == balance_box_level.getMPI().getRank());
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_mpi.getSize() > 1) {
         int compare_result;
         tbox::SAMRAI_MPI::Comm_compare(
            d_mpi.getCommunicator(),
            balance_box_level.getMPI().getCommunicator(),
            &compare_result);
         if (compare_result != MPI_CONGRUENT) {
            TBOX_ERROR("TreeLoadBalancer::loadBalanceBoxLevel:\n"
               << "The input balance_box_level has a SAMRAI_MPI that is\n"
               << "not congruent with the one set with setSAMRAI_MPI().\n"
               << "You must use freeMPICommunicator() before balancing\n"
               << "a BoxLevel with an incongruent SAMRAI_MPI.");
         }
      }
#endif
   }
   else {
      d_mpi = balance_box_level.getMPI();
   }

   if (d_print_steps ||
       d_print_break_steps) {
      tbox::plog << "TreeLoadBalancer::loadBalanceBoxLevel called with:"
                 << "\n  min_size = " << min_size
                 << "\n  max_size = " << max_size
                 << "\n  bad_interval = " << bad_interval
                 << "\n  cut_factor = " << cut_factor
                 << "\n  prebalance:\n"
                 << balance_box_level.format("  ", 2);
   }


   /*
    * Periodic image Box should be ignored during load balancing
    * because they have no real work.  The load-balanced results
    * should contain no periodic images.
    *
    * To avoid need for special logic to skip periodic images while
    * load balancing, we just remove periodic images in the
    * balance_box_level and all periodic edges in
    * anchor<==>balance.
    */

   balance_box_level.removePeriodicImageBoxes();
   if (balance_to_anchor) {

      balance_to_anchor->getTranspose().removePeriodicRelationships();
      balance_to_anchor->getTranspose().setHead(balance_box_level, true);

      balance_to_anchor->removePeriodicRelationships();
      balance_to_anchor->setBase(balance_box_level, true);

   }


   if (d_barrier_before) {
      t_barrier_before->start();
      d_mpi.Barrier();
      t_barrier_before->stop();
   }

   if (!rank_group.containsAllRanks()) {
      prebalanceBoxLevel(balance_box_level,
         balance_to_anchor,
         rank_group);
   }

   t_load_balance_box_level->start();

   d_min_size = min_size;
   d_max_size = max_size;
   d_bad_interval = bad_interval;
   d_cut_factor = cut_factor;
   d_min_load = d_min_size.getProduct();
   /*
    * Domain boxes are used by breakOffLoad to determine where
    * the bad cuts are.  Computing domain_boxes from domain_box_level
    * should be moved above the this method.
    */

   /*
    * We expect the domain box_level to be in globalized state.
    */
   TBOX_ASSERT(
      domain_box_level.getParallelState() ==
      hier::BoxLevel::GLOBALIZED);

   d_block_domain_boxes.clear();
   int nblocks =
      domain_box_level.getGridGeometry()->getNumberBlocks();
   d_block_domain_boxes.resize(nblocks);

   if (nblocks == 1) {
      domain_box_level.getGlobalBoxes(d_block_domain_boxes[0]);
      d_block_domain_boxes[0].refine(balance_box_level.getRefinementRatio());
   } else {
      for (int b = 0; b < nblocks; ++b) {
         d_block_domain_boxes[b] = hier::BoxContainer(
            domain_box_level.getGlobalBoxes(), hier::BlockId(b));

         d_block_domain_boxes[b].refine(balance_box_level.getRefinementRatio());
      }
   }
   /*
    * TODO: d_block_domain_boxes should be made a search tree and
    * findBadCutPoints should take a search tree form.  We do a lot
    * of searches through d_block_domain_boxes.  This would only make
    * a difference when the domain is described using many boxes.
    */


   t_compute_local_load->start();
   double local_load = computeLocalLoads(balance_box_level);
   t_compute_local_load->stop();

   LoadType max_local_load = local_load;

   LoadType global_sum_load = local_load;

   size_t nproc_with_initial_load =
      balance_box_level.getLocalNumberOfBoxes() > 0;

   {
      /*
       * Determine the total load and number of processes that has any
       * initial load.
       */
      t_compute_global_load->start();
      if (d_mpi.getSize() > 1) {
         double dtmp[2], dtmp_sum[2], dtmp_max[2];

         dtmp[0] = local_load;
         dtmp[1] = static_cast<double>(nproc_with_initial_load);
         d_mpi.Allreduce(dtmp,
            dtmp_sum,
            2,
            MPI_DOUBLE,
            MPI_SUM);
         global_sum_load = dtmp_sum[0];
         nproc_with_initial_load = (size_t)dtmp_sum[1];

         d_mpi.Allreduce(dtmp,
                         dtmp_max,
                         1,
                         MPI_DOUBLE,
                         MPI_MAX);
         max_local_load = dtmp_max[0];

      }
      t_compute_global_load->stop();
      if (d_print_steps) {
         tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
         tbox::plog.precision(6);
         tbox::plog << "TreeLoadBalancer::loadBalanceBoxLevel"
                    << " max_local_load=" << max_local_load
                    << " global_sum_load=" << global_sum_load
                    << " (initially born on "
                    << nproc_with_initial_load << " procs) across all "
                    << d_mpi.getSize()
                    << " procs, averaging " << global_sum_load / d_mpi.getSize()
                    << " or " << pow(global_sum_load / d_mpi.getSize(), 1.0 / d_dim.getValue())
                    << "^" << d_dim << " per proc." << std::endl;
      }
   }


   d_global_avg_load = global_sum_load / rank_group.size();


   /*
    * Add additional minimum box size restriction based on
    * d_min_load_fraction_per_box: Should be no smaller than a cubic
    * box that satisfies this work load.
    */
   if ( d_min_load_fraction_per_box > 0.0 ) {
      const hier::IntVector tmp_vec(d_min_size);

      int box_size_for_min_load_restriction =
         static_cast<int>(pow(d_global_avg_load*d_min_load_fraction_per_box,
                              1.0/d_dim.getValue())+ 0.5);
      d_min_size.max( hier::IntVector( d_dim, box_size_for_min_load_restriction ) );
      d_min_size.ceilingDivide(cut_factor);
      d_min_size *= cut_factor;
      d_min_load = d_min_size.getProduct();

      if (d_print_steps) {
         tbox::plog << "min_load_fraction_per_box changed min_size from " << tmp_vec;
         tbox::plog << " to " << d_min_size << '\n';
      }
   }


   /*
    * Compute how many root cycles to use based on severity of imbalance
    * using formula d_max_cycle_spread_ratio^number_of_cycles >= fanout_size.
    */
   const double fanout_size = max_local_load/d_global_avg_load;
   const int number_of_cycles =
      int(ceil( log(fanout_size)/log(d_max_cycle_spread_ratio) ));
      if (d_print_steps) {
         tbox::plog << "TreeLoadBalancer::loadBalanceBoxLevel"
                    << " max_cycle_spread_ratio=" << d_max_cycle_spread_ratio
                    << " fanout_size=" << fanout_size
                    << " number_of_cycles=" << number_of_cycles
                    << std::endl;
      }



   /*
    * The icycle loop spreads out the work each time through.  If
    * using more than one cycle, only the last one tries to balance
    * across all of d_mpi.
    */

   for (int icycle = 0; icycle < number_of_cycles; ++icycle) {

      // If not the first cycle, local_load needs updating.
      if (icycle > 0) {
         t_compute_local_load->start();
         local_load = computeLocalLoads(balance_box_level);
         t_compute_local_load->stop();
      }

      if (d_report_load_balance) {
         // Debugging: check overall load balance at intermediate cycles.
         tbox::plog
         << "TreeLoadBalancer::loadBalanceBoxLevel results before cycle "
         << icycle << ":" << std::endl;
         BalanceUtilities::gatherAndReportLoadBalance(
            local_load,
            balance_box_level.getMPI());
      }


      const bool last_cycle = (icycle == number_of_cycles-1);

      /*
       * Determine whether to use rank_group as is or subgroup it based
       * on cycles.
       */

      int number_of_groups = 1;
      int group_num = 0;

      tbox::RankGroup cycle_rank_group(d_mpi);
      if ( !last_cycle && rank_group.containsAllRanks() ) {
         createBalanceRankGroupBasedOnCycles(
            cycle_rank_group,
            number_of_groups,
            group_num,
            icycle,
            number_of_cycles);
      }


      /*
       * Compute the load for the group.  If this is the last cycle,
       * the group must include all processes, and the group's load
       * is the global sum load.  Else, use all-reduce to get the
       * group load.
       */
      t_compute_tree_load->start();

      double group_sum_load;

      if (icycle == number_of_cycles - 1) {

         group_sum_load = global_sum_load;

      } else {

         t_compute_tree_load_for_cycle[icycle]->start();

         /*
          * Use MPI's vector all-reduce to get individual group loads.
          * This gives more info than the process needs, but because the
          * number of groups << number of procs, it is still faster
          * (probably) than hand coded conmunication.
          */
         std::vector<double> group_loads(number_of_groups, 0.0);
         group_loads[group_num] = local_load;
         if (d_mpi.getSize() > 1) {
            d_mpi.AllReduce(&group_loads[0],
                            static_cast<int>(group_loads.size()),
                            MPI_SUM);
         }
         group_sum_load = group_loads[group_num];

         t_compute_tree_load_for_cycle[icycle]->stop();

      }

      t_compute_tree_load->stop();

      if (d_print_steps) {
         tbox::plog << "TreeLoadBalancer::loadBalanceBoxLevel"
                    << " cycle number=" << icycle
                    << " number_of_groups=" << number_of_groups
                    << " my group_num=" << group_num
                    << " my group size=" << cycle_rank_group.size()
                    << " my group_sum_load=" << group_sum_load
                    << std::endl;
      }

      /*
       * Compute the load-balancing map.
       */

      loadBalanceWithinRankGroup(
         balance_box_level,
         balance_to_anchor,
         rank_group,
         group_sum_load );

      if (d_barrier_after) {
         t_barrier_after->start();
         d_mpi.Barrier();
         t_barrier_after->stop();
      }

   }


   /*
    * Undo effects of min_load_fraction_per_box on d_min_size.  We do
    * not want that constraint during the remaining load balance
    * steps.
    */
   d_min_size = min_size;
   d_min_load = d_min_size.getProduct();


   /*
    * If max_size is given (positive), constrain boxes to the given
    * max_size.  If not given, skip the enforcement step to save some
    * communications.
    */

   if (max_size > hier::IntVector::getZero(d_dim)) {

      t_constrain_size->start();
      if (balance_to_anchor) {
         constrainMaxBoxSizes(
            balance_box_level,
            &balance_to_anchor->getTranspose());
      }
      else {
         constrainMaxBoxSizes(
            balance_box_level,
            balance_to_anchor);
      }
      t_constrain_size->stop();

      if (d_print_steps) {
         tbox::plog << " TreeLoadBalancer completed constraining box sizes."
                    << "\n";
      }

   }


   /*
    * Finished load balancing.  Clean up and wrap up.
    */

   d_min_size = hier::IntVector(d_dim, -1);
   d_max_size = hier::IntVector(d_dim, -1);
   d_block_domain_boxes.clear();
   d_bad_interval = hier::IntVector(d_dim, -1);
   d_cut_factor = hier::IntVector(d_dim, -1);

   t_load_balance_box_level->stop();

   local_load = computeLocalLoads(balance_box_level);
   d_load_stat.push_back(local_load);
   d_box_count_stat.push_back(
      static_cast<int>(balance_box_level.getBoxes().size()));

   if (d_print_steps) {
      tbox::plog << "Post balanced:\n" << balance_box_level.format("", 2);
   }

   if (d_report_load_balance) {
      t_report_loads->start();
      tbox::plog
      << "TreeLoadBalancer::loadBalanceBoxLevel results after "
      << number_of_cycles << " cycles:" << std::endl;
      BalanceUtilities::gatherAndReportLoadBalance(local_load,
         balance_box_level.getMPI());
      t_report_loads->stop();
   }

   if (d_check_connectivity && balance_to_anchor) {
      hier::Connector& anchor_to_balance = balance_to_anchor->getTranspose();
      tbox::plog << "TreeLoadBalancer checking balance-anchor connectivity."
                 << std::endl;
      int errs = 0;
      if (anchor_to_balance.checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in anchor_to_balance!\n";
      }
      if (balance_to_anchor->checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in balance_to_anchor!\n";
      }
      if (anchor_to_balance.checkTransposeCorrectness(*balance_to_anchor)) {
         ++errs;
         tbox::perr << "Error found in balance-anchor transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "anchor_box_level:\n" << anchor_to_balance.getBase().format("", 2)
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "anchor_to_balance:\n" << anchor_to_balance.format("", 2)
            << "balance_to_anchor:\n" << balance_to_anchor->format("", 2));
      }
      tbox::plog << "TreeLoadBalancer checked balance-anchor connectivity."
                 << std::endl;
   }

   if (d_barrier_after) {
      t_barrier_after->start();
      d_mpi.Barrier();
      t_barrier_after->stop();
   }

}



/*
 *************************************************************************
 * Constrain maximum box sizes in the given BoxLevel and
 * update given Connectors to the changed BoxLevel.
 *************************************************************************
 */
void
TreeLoadBalancer::constrainMaxBoxSizes(
   hier::BoxLevel& box_level,
   hier::Connector* anchor_to_level) const
{
   TBOX_ASSERT(!anchor_to_level || anchor_to_level->hasTranspose());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, box_level);

   t_map_big_boxes->start();

   if (d_print_break_steps) {
      tbox::plog << "Mapping oversized boxes starting with "
                 << box_level.getBoxes().size() << " boxes."
                 << std::endl;
   }

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   hier::BoxLevel constrained(box_level.getRefinementRatio(),
      box_level.getGridGeometry(),
      box_level.getMPI());
   hier::MappingConnector unconstrained_to_constrained(box_level,
      constrained,
      zero_vector);

   const hier::BoxContainer& unconstrained_boxes = box_level.getBoxes();

   hier::LocalId next_available_index = box_level.getLastLocalId() + 1;

   for (hier::BoxContainer::const_iterator ni = unconstrained_boxes.begin();
        ni != unconstrained_boxes.end(); ++ni) {

      const hier::Box& box = *ni;

      const hier::IntVector box_size = box.numberCells();

      /*
       * If box already conform to max size constraint, keep it.
       * Else chop it up and keep the parts.
       */

      if (box_size <= d_max_size) {

         if (d_print_break_steps) {
            tbox::plog << "    Not oversized: " << box
                       << box.numberCells() << "\n";
         }
         constrained.addBox(box);

      } else {

         if (d_print_break_steps) {
            tbox::plog << "    Breaking oversized " << box
                       << box.numberCells() << " ->";
         }
         hier::BoxContainer chopped(box);
         hier::BoxUtilities::chopBoxes(
            chopped,
            d_max_size,
            d_min_size,
            d_cut_factor,
            d_bad_interval,
            d_block_domain_boxes[box.getBlockId().getBlockValue()]);
         TBOX_ASSERT( !chopped.isEmpty() );

         if (chopped.size() != 1) {

            hier::Connector::NeighborhoodIterator base_box_itr =
               unconstrained_to_constrained.makeEmptyLocalNeighborhood(
                  box.getBoxId());

            for (hier::BoxContainer::iterator li = chopped.begin();
                 li != chopped.end(); ++li) {

               const hier::Box fragment = *li;

               const hier::Box new_box(fragment,
                                       next_available_index++,
                                       d_mpi.getRank());
               TBOX_ASSERT(new_box.getBlockId() == ni->getBlockId());

               if (d_print_break_steps) {
                  tbox::plog << "  " << new_box
                             << new_box.numberCells();
               }

               constrained.addBox(new_box);

               unconstrained_to_constrained.insertLocalNeighbor(
                  new_box,
                  base_box_itr);

            }

            if (d_print_break_steps) {
               tbox::plog << "\n";
            }

         } else {
            TBOX_ASSERT( box.isSpatiallyEqual( chopped.front() ) );
            if (d_print_break_steps) {
               tbox::plog << " Unbreakable!" << "\n";
            }
            constrained.addBox(box);
         }

      }

   }

   if (d_print_steps) {
      tbox::plog
      << " TreeLoadBalancer::constrainMaxBoxSizes completed building unconstrained_to_constrained"
      << "\n";
   }

   if (anchor_to_level && anchor_to_level->isFinalized()) {
      // Modify anchor<==>level Connectors and swap box_level with constrained.
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix(d_object_name);
      mca.modify(*anchor_to_level,
                 unconstrained_to_constrained,
                 &box_level,
                 &constrained);
   } else {
      // Swap box_level and constrained without touching anchor<==>level.
      hier::BoxLevel::swap(box_level, constrained);
   }

   t_map_big_boxes->stop();
}



/*
 *************************************************************************
 * Given an "unbalanced" BoxLevel, compute the BoxLevel that is
 * load-balanced and compute the mapping between the unbalanced and
 * balanced BoxLevels.
 *
 * If given a RankGroup with less than all ranks, we treat it as a
 * specific user request to balance only within the RankGroup and just
 * use the RankGroup as is.  Otherwise, we may generate sub-groups
 * based on the cycle number and balance within the generated
 * sub-group.
 *
 * The objective of balancing over multiple cycles is to avoid
 * unscalable performance in the cases where just a few processes own
 * all the initial load.  By slowly spreading out the load, no process
 * has to set up Connector unbalanced_to_balanced with number of
 * relationships that scales with the machine size.
 *
 * If the local process is not a member of the RankGroup, it does not
 * participate in the work and just sets the output objects to be
 * locally empty.
 *************************************************************************
 */
void
TreeLoadBalancer::loadBalanceWithinRankGroup(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_anchor,
   const tbox::RankGroup& rank_group,
   const double group_sum_load ) const
{
   TBOX_ASSERT(!balance_to_anchor || balance_to_anchor->hasTranspose());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, balance_box_level);

   double group_avg_load = group_sum_load / rank_group.size();

   /*
    * Initialize empty balanced_box_level and mappings.
    */
   hier::BoxLevel balanced_box_level(
      balance_box_level.getRefinementRatio(),
      balance_box_level.getGridGeometry(),
      balance_box_level.getMPI());
   hier::MappingConnector balanced_to_unbalanced(balanced_box_level,
         balance_box_level,
         hier::IntVector::getZero(d_dim));
   hier::MappingConnector unbalanced_to_balanced(balance_box_level,
         balanced_box_level,
         hier::IntVector::getZero(d_dim));
   unbalanced_to_balanced.setTranspose(&balanced_to_unbalanced, false);


   if ( !rank_group.isMember(d_mpi.getRank()) ) {
      /*
       * The following assert should be guaranteed by an earlier call
       * to prebalanceBoxLevel.  Having boxes without being in the
       * given rank group leads to undefined results.
       */
      TBOX_ASSERT( balance_box_level.getLocalNumberOfBoxes() == 0 );

t_post_load_distribution_barrier->start();
d_mpi.Barrier(); // Temporary barrier to determine if the follow communication slows down unfinished communications in load distribution phase.
t_post_load_distribution_barrier->stop();
      if (balance_to_anchor && balance_to_anchor->hasTranspose()) {
         hier::MappingConnectorAlgorithm mca;
         mca.setTimerPrefix(d_object_name);
         t_use_map->start();
         mca.modify(
            balance_to_anchor->getTranspose(),
            unbalanced_to_balanced,
            &balance_box_level,
            &balanced_box_level);
         t_use_map->stop();
      } else {
         hier::BoxLevel::swap(balance_box_level, balanced_box_level);
      }
      return;
   }


   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << "TreeLoadBalancer::LoadBalanceWithinRankGroup balancing "
                 << group_sum_load << " units in group of "
                 << d_mpi.getSize() << " procs, averaging " << group_avg_load
                 << " or " << pow(group_avg_load, 1.0 / d_dim.getValue())
                 << "^" << d_dim << " per proc."
                 << "  Avg is " << group_avg_load/d_min_size.getProduct()
                 << " times min size of " << d_min_size
                 << std::endl;
   }


   t_get_map->start();

   t_load_distribution->start();

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


   /*
    * Arrange the group ranks in a rank tree in order to get
    * the parent/children in the group.
    */
   d_rank_tree->setupTree( rank_group, d_mpi.getRank() );
   if ( d_print_steps ) {
      // Write local part of tree to log.
      tbox::plog << "TreeLoadBalancer tree:\n"
                 << "  Root rank: " << d_rank_tree->getRootRank() << '\n'
                 << "  Child number: " << d_rank_tree->getChildNumber() << '\n'
                 << "  Generation number: " << d_rank_tree->getGenerationNumber() << '\n'
                 << "  Number of children: " << d_rank_tree->getNumberOfChildren() << '\n'
                 << "  Local relatives: "
                 << "  " << d_rank_tree->getParentRank() << " <- [" << d_rank_tree->getRank() << "] -> {";
      for ( unsigned int i=0; i<d_rank_tree->getNumberOfChildren(); ++i ) {
         tbox::plog << ' ' << d_rank_tree->getChildRank(i);
      }
      tbox::plog << " }" << std::endl;
   }

   const int num_children = d_rank_tree->getNumberOfChildren();


   /*
    * Communication objects for sending to/receiving from
    * parent/children: We could combine all of these AsyncCommStages
    * and most of the AsyncCommPeers, but we intentionally keep them
    * separate to aid performance analysis.
    */

   tbox::AsyncCommStage child_send_stage;
   tbox::AsyncCommPeer<char>* child_sends = 0;
   tbox::AsyncCommStage parent_send_stage;
   tbox::AsyncCommPeer<char>* parent_send = 0;

   setupAsyncCommObjects(
      child_send_stage,
      child_sends,
      parent_send_stage,
      parent_send,
      rank_group );
   child_send_stage.setCommunicationWaitTimer(t_child_send_wait);
   parent_send_stage.setCommunicationWaitTimer(t_parent_send_wait);

   tbox::AsyncCommStage child_recv_stage;
   tbox::AsyncCommPeer<char>* child_recvs = 0;
   tbox::AsyncCommStage parent_recv_stage;
   tbox::AsyncCommPeer<char>* parent_recv = 0;

   setupAsyncCommObjects(
      child_recv_stage,
      child_recvs,
      parent_recv_stage,
      parent_recv,
      rank_group );
   child_recv_stage.setCommunicationWaitTimer(t_child_recv_wait);
   parent_recv_stage.setCommunicationWaitTimer(t_parent_recv_wait);



   /*
    * Outline of the tree load balancing algorithm as implemented:
    *
    * 1. For each child of the local process:
    * Receive data from subtree rooted at child (number in
    * subtree, excess work, remaining work in subtree, etc.).
    *
    * 2. Compute data for subtree rooted at self by combining
    * local data with children subtree data.
    *
    * 3. If parent exists:
    * Send subtree info (number in subtree, excess work,
    * remaining work in subtree, etc.) to parent.
    *
    * 4. If parent exists and we need more work:
    * Receive additional work from parent.
    *
    * 5. Partition additional work among children and self.
    *
    * 6. For each child:
    * Send additional work (if any).
    */

   /*
    * Step 1:
    *
    * Post receive for data from subtree rooted at children.
    * We have to do a few local setups, but post the receive
    * now to overlap communication.
    */
   t_get_load_from_children->start();
   for (int c = 0; c < num_children; ++c) {
      child_recvs[c].setRecvTimer(t_child_recv_wait);
      child_recvs[c].setWaitTimer(t_child_recv_wait);
      child_recvs[c].beginRecv();
      if (child_recvs[c].isDone()) {
         child_recvs[c].pushToCompletionQueue();
      }
   }
   t_get_load_from_children->stop();


   /*
    * Step 2, local part:
    *
    * The local process must generate indices for new and imported
    * boxes.  To do it deterministically, no generated index should
    * depend on message arrival order.  To achieve this, we maintain
    * 2+deg values in next_available_index: one for the local
    * process, one for the parent and one for each child.  deg is the
    * degree of tree d_rank_tree.  The first
    * index given to a locally generated Box is some index unused by
    * balance_box_level.  The first index given to a Box
    * from the parent is the same value plus 1.  The first index given
    * to a box from child 0 is the same value plus 2.  And so on.
    * Each time a value from next_available_index is used, we
    * increment it by 2+deg so that the 2+deg available
    * values can never be the same.  Moreover, boxes from a certain
    * source always take indices from its own set, independent of when
    * boxes from other sources arrive.
    */
   std::vector<hier::LocalId> next_available_index(2 + d_rank_tree->getDegree());
   next_available_index[0] = balance_box_level.getLastLocalId() + 1;

   /*
    * The next line makes next_available_index[0] divisible by 2+deg.
    * It is not strictly necessary but makes debugging much easier because
    * we can quickly associate any value with the source of its Box.
    */
   next_available_index[0] +=
      hier::LocalId(2+d_rank_tree->getDegree()) - (next_available_index[0] % (2 + d_rank_tree->getDegree()));

   for (unsigned int c = 1; c < d_rank_tree->getDegree() + 2; ++c) {
      next_available_index[c] = next_available_index[0] + c;
   }


   /*
    * Data for storing and transfering subtree info.
    */
   SubtreeData my_subtree;
   std::vector<SubtreeData> child_subtrees(num_children);


   /*
    * unassigned is a container of BoxInTransit that has been released by
    * a process and has not yet been assigned to another.  First, put
    * excess local work (if any) in unassigned.  Imported
    * BoxInTransits are placed here before determining whether to keep
    * them or send them to another part of the tree.
    */
   TransitSet unassigned(balance_box_level.getBoxes().begin(),
                         balance_box_level.getBoxes().end());


   /*
    * Compute local proc's load and store in
    * my_subtree.  This will eventually include data for the subtree.
    * We will add the rest of the subtree's work when we receive that
    * data from the children.
    */
   my_subtree.d_num_procs = 1;
   my_subtree.d_subtree_load_ideal = group_avg_load;
   my_subtree.d_subtree_load_current = unassigned.getSumLoad();
   my_subtree.d_subtree_load_upperlimit = group_avg_load*(1+d_flexible_load_tol);

   my_subtree.d_eff_num_procs = my_subtree.d_num_procs;
   my_subtree.d_eff_load_ideal = my_subtree.d_subtree_load_ideal;
   my_subtree.d_eff_load_current = my_subtree.d_subtree_load_current;
   my_subtree.d_eff_load_upperlimit = my_subtree.d_subtree_load_upperlimit;

   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << "Initial local load is " << unassigned.getSumLoad()
                 << " (" << (unassigned.getSumLoad()/group_avg_load)
                 << ") in " << unassigned.size() << " boxes"
                 << ", surplus = " << my_subtree.surplus()
                 << ", excess = " << my_subtree.excess()
                 << std::endl;
   }

   /*
    * The incoming unbalanced boxes need a mapping to describe their
    * change, but we don't know what they will become, so create empty
    * maps for now.  Should any not change, we'll erase their
    * neighborhood later.
    */
   for (TransitSet::const_iterator ni=unassigned.begin(); ni!=unassigned.end(); ++ni ) {
      unbalanced_to_balanced.makeEmptyLocalNeighborhood(ni->d_orig_box.getBoxId());
   }



   /*
    * Step 2, remote part:
    *
    * Finish getting tree and load data from children.
    * Add imported BoxInTransit to unassigned bin.
    */
   t_get_load_from_children->start();
   while ( child_recv_stage.numberOfCompletedMembers() > 0 ||
           child_recv_stage.advanceSome() ) {

      tbox::AsyncCommPeer<char>* child_recv =
         CPP_CAST<tbox::AsyncCommPeer<char> *>(child_recv_stage.popCompletionQueue());

      TBOX_ASSERT(child_recv != 0);
      TBOX_ASSERT(child_recv >= child_recvs);
      TBOX_ASSERT(child_recv < child_recvs + num_children);

      const int cindex = static_cast<int>(child_recv - child_recvs);

      TBOX_ASSERT(cindex >= 0 && cindex < num_children);

      /*
       * Extract data from the child cindex, storing it in
       * child_subtrees[cindex].
       */
      tbox::MessageStream mstream(child_recv->getRecvSize(),
                                  tbox::MessageStream::Read,
                                  child_recv->getRecvData(),
                                  false);
      unpackSubtreeDataUp(
         child_subtrees[cindex],
         next_available_index[cindex],
         mstream);

      unassigned.insert( child_subtrees[cindex].d_work_traded.begin(),
                         child_subtrees[cindex].d_work_traded.end() );

      if (d_print_steps) {
         tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
         tbox::plog.precision(6);
         tbox::plog << "Received from child "
                    << cindex << ':' << d_rank_tree->getChildRank(cindex) << ": "
                    << child_subtrees[cindex].d_work_traded.size() << " boxes ("
                    << child_subtrees[cindex].d_work_traded.getSumLoad() << " units):";
         for ( TransitSet::const_iterator ni=child_subtrees[cindex].d_work_traded.begin();
               ni!=child_subtrees[cindex].d_work_traded.end(); ++ni ) {
            const BoxInTransit& box_in_transit = *ni;
            tbox::plog << "  " << box_in_transit;
         }
         tbox::plog << std::endl;
      }

      my_subtree.addChild( child_subtrees[cindex] );

   }

   size_t unassigned_highwater = unassigned.size();


   // We should have received everything at this point.
   TBOX_ASSERT(!child_recv_stage.hasPendingRequests());


   if ( my_subtree.effDeficit() > 0 && !d_rank_tree->isRoot() ) {
      my_subtree.d_wants_work_from_parent = true;
   }


   if (d_print_steps) {
      tbox::plog << "Received children subtree data." << std::endl;
      for (int c = 0; c < num_children; ++c) {
         tbox::plog << "Child "
                    << c << ':' << d_rank_tree->getChildRank(c)
                    << " subtree:\n";
         child_subtrees[c].printClassData( "  ", tbox::plog );
      }
      tbox::plog << "Initial subtree:\n";
      my_subtree.printClassData( "  ", tbox::plog );
      tbox::plog << "unassigned has: " << unassigned.size()
                 << " boxes (" << unassigned.getSumLoad() << " units)."
                 << std::endl;
   }

   t_get_load_from_children->stop();


   /*
    * Step 3:
    *
    * Send subtree info and excess work (if any) up to parent.
    */
   t_send_load_to_parent->start();
   if (parent_send != 0) {

      if ( my_subtree.effExcess() > 0 ) {
         /*
          * Don't send more than the surplus because that would
          * overload the complement of the subtree.  Don't send less
          * than effective excess because that would overload the
          * subtree.  Sometimes underloaded children subtrees cause,
          * effective excess > surplus.  In these cases, don't send
          * the effective excess, because that would overload the
          * complement of the subtree.  It is better to overload the
          * subtree than to progressively push surplus up, making the
          * root extremely overloaded.  Keeping the overload in the
          * subtree is better for data locality.
          */
         const LoadType export_load_low = tbox::MathUtilities<double>::Min(my_subtree.effExcess(), my_subtree.surplus());
         const LoadType export_load_high = my_subtree.surplus();
         const LoadType export_load_ideal = export_load_low;

         if ( export_load_low > 0 ) {

            if (d_print_steps) {
               tbox::plog << "Attempting to reassign " << export_load_ideal
                          << " [" << export_load_low << ", " << export_load_high
                          << "] of unassigned load to parent.\n";
            }

            LoadType export_load_actual = adjustLoad(
               my_subtree.d_work_traded /* to parent */,
               unassigned,
               next_available_index[d_rank_tree->getDegree()],
               export_load_ideal,
               export_load_low,
               export_load_high );
            TBOX_ASSERT( export_load_actual >= 0 );

            my_subtree.d_subtree_load_current -= my_subtree.d_work_traded.getSumLoad();
            my_subtree.d_eff_load_current -= my_subtree.d_work_traded.getSumLoad();

            if (d_print_steps) {
               tbox::plog << "Assigned " << my_subtree.d_work_traded.size()
                          << " boxes (" << export_load_actual << " / ["
                          << export_load_low << ", " << export_load_high << "] "
                          << " units) to parent's export bin:";
               for (TransitSet::const_iterator ni = my_subtree.d_work_traded.begin();
                    ni!=my_subtree.d_work_traded.end(); ++ni) {
                  tbox::plog << "  " << *ni;
               }
               tbox::plog << std::endl;
            }

         }

      }

      /*
       * Send local subtree info, along with any exported work,
       * up to parent.
       */
      tbox::MessageStream mstream;
      packSubtreeDataUp(mstream, my_subtree);
      if (d_print_steps) {
         tbox::plog << "Sending to parent " << d_rank_tree->getParentRank() << ": "
                    << my_subtree.d_work_traded.size() << " boxes ("
                    << my_subtree.d_work_traded.getSumLoad() << " units)."
                    << "  message length = " << mstream.getCurrentSize() << " bytes"
                    << std::endl;
         tbox::plog << "unassigned has: " << unassigned.size()
                    << " boxes (" << unassigned.getSumLoad() << " units)."
                    << std::endl;
      }
      parent_send->setSendTimer(t_parent_send_wait);
      parent_send->setWaitTimer(t_parent_send_wait);
      parent_send->beginSend(static_cast<const char*>(mstream.getBufferStart()),
                             static_cast<int>(mstream.getCurrentSize()));

   }
   t_send_load_to_parent->stop();



   /*
    * Step 4:
    *
    * Finish the send-up and begin the send-down.
    */
   if (my_subtree.d_wants_work_from_parent) {
      t_parent_load_comm->start();
      t_get_load_from_parent->start();

      parent_recv->setRecvTimer(t_parent_recv_wait);
      parent_recv->setWaitTimer(t_parent_recv_wait);
      parent_recv->beginRecv();

      t_get_load_from_parent->stop();
      t_parent_load_comm->stop();
   }


   /*
    * May do some things here that do not depend on message from
    * parents.  Steve Smith suggested sending work down to underloaded
    * children at this point, even if it makes the local process
    * underloaded as a result.  The local process can recover the
    * correct amount of work when it comes down from the parent.  This
    * would allow some children subtrees to wait less, but it has
    * other consequences.
    */


   if (my_subtree.d_wants_work_from_parent) {

      /*
       * Receive and unpack message from parent.
       */
      t_get_load_from_parent->start();

      parent_recv->completeCurrentOperation();

      tbox::MessageStream mstream(parent_recv->getRecvSize(),
                                  tbox::MessageStream::Read,
                                  parent_recv->getRecvData(),
                                  false);
      unpackSubtreeDataDown(
         my_subtree,
         next_available_index[1 + d_rank_tree->getDegree()],
         mstream);

      unassigned.insert( my_subtree.d_work_traded.begin(),
                         my_subtree.d_work_traded.end() );
      my_subtree.d_subtree_load_current += my_subtree.d_work_traded.getSumLoad();
      my_subtree.d_eff_load_current += my_subtree.d_work_traded.getSumLoad();

      if ( unassigned_highwater < unassigned.size() ) {
         unassigned_highwater = unassigned.size();
      }

      if (d_print_steps) {
         if (d_print_steps) {
            tbox::plog << "Received from parent " << d_rank_tree->getParentRank() << ":"
                       << my_subtree.d_work_traded.size() << " boxes ("
                       << my_subtree.d_work_traded.getSumLoad() << " units):";
            for ( TransitSet::const_iterator ni=my_subtree.d_work_traded.begin();
                  ni!=my_subtree.d_work_traded.end(); ++ni ) {
               const BoxInTransit& box_in_transit = *ni;
               tbox::plog << "  " << box_in_transit;
            }
            tbox::plog << std::endl;
         }
      }

      t_get_load_from_parent->stop();
   }
   else {
      if (d_print_steps) {
         tbox::plog << "Did not request work from parent.\n";
      }
   }

   if (d_print_steps) {
      tbox::plog << "Postparent subtree:\n";
      my_subtree.printClassData( "  ", tbox::plog );
      tbox::plog << "unassigned has: " << unassigned.size()
                 << " boxes (" << unassigned.getSumLoad() << " units)."
                 << std::endl;
   }



   /*
    * Step 5 and 6:
    *
    * Reassign and send work to each child that requested work.
    */

   t_send_load_to_children->start();

   for (int ichild = 0; ichild < num_children; ++ichild) {

      SubtreeData& recip_subtree = child_subtrees[ichild];

      if (recip_subtree.d_wants_work_from_parent) {

         const LoadType surplus_per_eff_des =
            computeSurplusPerEffectiveDescendent(
               unassigned,
               group_avg_load,
               child_subtrees,
               ichild );

         const LoadType export_load_ideal = recip_subtree.effDeficit()
            + (surplus_per_eff_des < 0.0 ? 0.0 :
               surplus_per_eff_des*recip_subtree.d_eff_num_procs);

         const LoadType export_load_low = recip_subtree.effDeficit()
            + surplus_per_eff_des*recip_subtree.d_eff_num_procs;

         const LoadType export_load_high =
            tbox::MathUtilities<double>::Max(export_load_ideal,
                                             recip_subtree.effMargin());

         TBOX_ASSERT( export_load_high >= export_load_ideal );
         TBOX_ASSERT( export_load_ideal >= export_load_low );

         if ( export_load_low > 0.0 ) {

            if (d_print_steps) {
               tbox::plog << "Adjusting export bin for child "
                          << ichild << ':' << d_rank_tree->getChildRank(ichild)
                          << " to " << export_load_ideal
                          << " [" << export_load_low << ", " << export_load_high << "]\n";
            }

            const LoadType export_load_actual = adjustLoad(
               recip_subtree.d_work_traded,
               unassigned,
               next_available_index[d_rank_tree->getDegree()],
               export_load_ideal,
               export_load_low,
               export_load_high );
            TBOX_ASSERT(export_load_actual >= 0);
            recip_subtree.d_subtree_load_current += export_load_actual;
            recip_subtree.d_eff_load_current += export_load_actual;

            if (d_print_steps) {
               tbox::plog << "Assigned " << recip_subtree.d_work_traded.size()
                          << " boxes (" << export_load_actual << " / " << export_load_ideal
                          << " [" << export_load_low << ", " << export_load_high << "] "
                          << " units) to child "
                          << ichild << ':' << d_rank_tree->getChildRank(ichild)
                          << " for " << recip_subtree.d_num_procs << " procs:";
               for (TransitSet::const_iterator ni = recip_subtree.d_work_traded.begin();
                    ni != recip_subtree.d_work_traded.end(); ++ni) {
                  tbox::plog << "  " << *ni;
               }
               tbox::plog << std::endl;
            }

         }

         tbox::MessageStream mstream;
         packSubtreeDataDown(mstream, recip_subtree);
         if (d_print_steps) {
            tbox::plog << "Sending to child "
                       << ichild << ':' << d_rank_tree->getChildRank(ichild)
                       << ' ' << recip_subtree.d_work_traded.size() << " boxes ("
                       << recip_subtree.d_work_traded.getSumLoad() << " units)."
                       << "  message length = " << mstream.getCurrentSize() << " bytes"
                       << std::endl;
            tbox::plog << "unassigned has: " << unassigned.size()
                       << " boxes (" << unassigned.getSumLoad() << " units)."
                       << std::endl;
         }
         child_sends[ichild].setSendTimer(t_child_send_wait);
         child_sends[ichild].setWaitTimer(t_child_send_wait);
         child_sends[ichild].beginSend(static_cast<const char*>(mstream.getBufferStart()),
                                       static_cast<int>(mstream.getCurrentSize()));

      }

   }

   t_send_load_to_children->stop();


   if (d_print_steps) {
      tbox::plog << "After settling parent and children, unassigned is "
                 << unassigned.getSumLoad() << " ("
                 << (unassigned.getSumLoad()/group_avg_load)
                 << ") units in " << unassigned.size() << " boxes:\n";
      for ( TransitSet::const_iterator bi=unassigned.begin();
            bi!=unassigned.end(); ++bi ) {
         tbox::plog << "    " << *bi << std::endl;
      }
      tbox::plog << "Final local load is " << unassigned.getSumLoad()
                 << " (" << (unassigned.getSumLoad()/group_avg_load)
                 << ") in " << unassigned.size() << " boxes"
                 << ", surplus = " << (unassigned.getSumLoad()-group_avg_load)
                 << ", excess = " << (unassigned.getSumLoad()-(1+d_flexible_load_tol)*group_avg_load)
                 << std::endl;

      tbox::plog << "Final subtree:\n";
      my_subtree.printClassData( "  ", tbox::plog );

      for (int ichild = 0; ichild < num_children; ++ichild) {
         tbox::plog << "Final child "
                    << ichild << ":" << d_rank_tree->getChildRank(ichild)
                    << " subtree:\n";
         child_subtrees[ichild].printClassData("  ", tbox::plog );
      }
   }

   /*
    * All unassigned boxes should go into balanced_box_level.  Put
    * them there and generate relationships in balanced<==>unbalanced
    * mapping Connectors where required.
    *
    * We can generate balanced--->unbalanced edges for all unassigned
    * boxes because we have their origin info.  If the unassigned box
    * originated locally, we can generate the unbalanced--->balanced
    * edge for them as well.  However, we can't generate these edges
    * for boxes originating remotely.  For these boxes, leave them in
    * unassigned for the step of notifying their origin owners that we
    * have them.  Otherwise, remove boxes from unassigned.
    */
   t_local_balancing->start();

   for (TransitSet::iterator
        ni = unassigned.begin();
        ni != unassigned.end(); /* incremented in loop */) {

      const BoxInTransit& box_in_transit = *ni;
      balanced_box_level.addBox(box_in_transit.d_box);

      if (box_in_transit.d_box.isIdEqual(box_in_transit.d_orig_box)) {
         // Unchanged box requires no mapping.  Nothing else need to be done.
         TBOX_ASSERT( unbalanced_to_balanced.isEmptyNeighborhood(ni->d_box.getBoxId()) );
         unbalanced_to_balanced.eraseLocalNeighborhood(ni->d_box.getBoxId());
         unassigned.erase(ni++);
      } else {

         balanced_to_unbalanced.insertLocalNeighbor(
            box_in_transit.d_orig_box,
            box_in_transit.d_box.getBoxId());

         if (box_in_transit.d_orig_box.getOwnerRank() == d_mpi.getRank()) {
            unbalanced_to_balanced.insertLocalNeighbor(
               box_in_transit.d_box,
               box_in_transit.d_orig_box.getBoxId());
            unassigned.erase(ni++);
         }
         else {
            // Leave this box in unassigned for notifying originating
            // process of where it landed.
            ++ni;
         }
      }

   }

   t_local_balancing->stop();

   t_load_distribution->stop();


   /*
    * Finish messages before starting edge info exchange.
    * We have only sends to complete, so it should not take
    * long to advance them all to completion.
    */
   t_finish_sends->start();
   child_send_stage.advanceAll();
   parent_send_stage.advanceAll();
   t_finish_sends->stop();
   child_send_stage.clearCompletionQueue();
   parent_send_stage.clearCompletionQueue();
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < num_children; ++i) {
      TBOX_ASSERT(child_sends[i].isDone());
      TBOX_ASSERT(child_recvs[i].isDone());
   }
   if (parent_send != 0) {
      TBOX_ASSERT(parent_send->isDone());
      TBOX_ASSERT(parent_recv->isDone());
   }
#endif


   constructSemilocalUnbalancedToBalanced(
      unbalanced_to_balanced,
      unassigned );

   if ( d_summarize_map ) {
      tbox::plog << "TreeLoadBalancer::loadBalanceWithinRankGroup unbalanced--->balanced map:\n"
                 << unbalanced_to_balanced.format("\t",0)
                 << "Map statistics:\n" << unbalanced_to_balanced.formatStatistics("\t")
                 << "TreeLoadBalancer::loadBalanceWithinRankGroup balanced--->unbalanced map:\n"
                 << balanced_to_unbalanced.format("\t",0)
                 << "Map statistics:\n" << balanced_to_unbalanced.formatStatistics("\t")
                 << '\n';
   }


   if (d_check_connectivity) {
      tbox::plog
      << "TreeLoadBalancer checking unbalanced-balanced connectivity."
      << std::endl;
      int errs = 0;
      if (unbalanced_to_balanced.checkOverlapCorrectness(true, true)) {
         ++errs;
         tbox::perr << "Error found in unbalanced_to_balanced!\n";
      }
      if (balanced_to_unbalanced.checkOverlapCorrectness(true, true)) {
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
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "balanced_box_level:\n" << balanced_box_level.format("", 2)
            << "unbalanced_to_balanced:\n" << unbalanced_to_balanced.format("", 2)
            << "balanced_to_unbalanced:\n" << balanced_to_unbalanced.format("", 2));
      }
   }

   if (d_check_map) {
      if (unbalanced_to_balanced.findMappingErrors() != 0) {
         TBOX_ERROR(
            "TreeLoadBalancer::loadBalanceWithinRankGroup Mapping errors found in unbalanced_to_balanced!");
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         TBOX_ERROR(
            "TreeLoadBalancer::loadBalanceWithinRankGroup Transpose errors found!");
      }
   }


   t_get_map->stop();

   if (balance_to_anchor && balance_to_anchor->hasTranspose()) {
      t_use_map->start();
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix(d_object_name);
      mca.modify(
         balance_to_anchor->getTranspose(),
         unbalanced_to_balanced,
         &balance_box_level,
         &balanced_box_level);
      t_use_map->stop();
   } else {
      hier::BoxLevel::swap(balance_box_level, balanced_box_level);
   }


   if ( d_comm_graph_writer ) {
      /*
       * Record these edges:
       * - Two edges to parent: load traded, and boxes traded.
       * - Same two edges from parent, plus one more for message size.
       * - Two timer for edges from children, one from parent.
       * Record these nodes:
       * - Number of final boxes:
       */
      d_comm_graph_writer->addRecord( d_mpi, int(0), size_t(8), size_t(7) );

      const int prank = (d_rank_tree->isRoot() ? -1 : d_rank_tree->getParentRank());

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(0),
         "load up",
         double(my_subtree.d_wants_work_from_parent ? 0 : my_subtree.d_work_traded.getSumLoad()),
         tbox::CommGraphWriter::TO,
         prank );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(1),
         "boxes up",
         double(my_subtree.d_wants_work_from_parent ? 0 : my_subtree.d_work_traded.size()),
         tbox::CommGraphWriter::TO,
         prank );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(2),
         "load down",
         double(my_subtree.d_wants_work_from_parent ? my_subtree.d_work_traded.getSumLoad() : 0),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(3),
         "boxes down",
         double(my_subtree.d_wants_work_from_parent ? my_subtree.d_work_traded.size() : 0),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(4),
         "bytes down",
         double(my_subtree.d_wants_work_from_parent ? parent_recv->getRecvSize() : int(0)),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(5),
         "child wait",
         t_child_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         d_rank_tree->getChildRank(0) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(6),
         "child wait",
         t_child_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         d_rank_tree->getChildRank(1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(7),
         "parent wait",
         t_parent_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(0),
         "initial box count",
         double(balance_box_level.getLocalNumberOfBoxes()) );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(1),
         "initial load",
         double(balance_box_level.getLocalNumberOfCells())/group_avg_load );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(2),
         "final box count",
         double(balanced_box_level.getLocalNumberOfBoxes()) );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(3),
         "final surplus",
         double(balanced_box_level.getLocalNumberOfCells())-group_avg_load );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(4),
         "final load",
         double(balanced_box_level.getLocalNumberOfCells())/group_avg_load );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(5),
         "subtree surplus",
         double(my_subtree.surplus()) );

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(6),
         "unassigned highwater",
         double(unassigned_highwater) );

   }


   if (d_print_steps) {
      tbox::plog << "TreeLoadBalancer::LoadBalanceWithinRankGroup returning"
                 << std::endl;
   }

   destroyAsyncCommObjects(child_sends, parent_send);
   destroyAsyncCommObjects(child_recvs, parent_recv);

   return;
}



/*
 *************************************************************************
 *
 * This method adjusts the load in a main_bin of BoxInTransits by
 * moving work between it and a holding bin.  It tries to bring
 * main_bin's load to the specified ideal_load.
 *
 * The high_load and low_load define an acceptable range around the
 * ideal_load.  As soon as the main load falls in this range, no
 * further change is tried, even if it may bring the load closer to
 * the ideal.
 *
 * This method makes a best effort and returns the amount of load
 * moved.  It can move BoxInTransit between given sets and, if needed,
 * break some BoxInTransit up to move part of the work.
 *
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 * Return amount of load moved from main_bin to hold_bin.  Negative
 * amount means load moved from hold_bin to main_bin.
 *
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::adjustLoad(
   TransitSet& main_bin,
   TransitSet& hold_bin,
   hier::LocalId& next_available_index,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load ) const
{
   if (d_print_steps) {
      tbox::plog << "  adjustLoad attempting to bring main load from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " or within [" << low_load << ", " << high_load << "]."
                 << std::endl;
   }
   TBOX_ASSERT( low_load <= ideal_load );
   TBOX_ASSERT( high_load >= ideal_load );


   LoadType actual_transfer = 0;

   if ((main_bin.empty() && ideal_load <= 0 ) ||
       (hold_bin.empty() && main_bin.getSumLoad() < ideal_load )) {
      return actual_transfer;
   }

   t_adjust_load->start();

   /*
    * The algorithm cycles through a do-loop.  Each time around, we
    * try to swap some BoxInTransit between main_bin and hold_bin
    * until we have main_bin's load in [low_load,high_load] or we
    * cannot improve the actual_transfer any further.  Then, we try
    * breaking up a BoxInTransit to improve the results.  If we break
    * some BoxInTransit, we generate some more swapping options that
    * were not there before, so we loop back to try swapping again.
    *
    * If a break phase does not break any Box (and does not generate
    * more swap options), the loop will stop making changes.  We exit
    * the loop at that point (and whenever we get main_bin's load in
    * the correct range).
    */
   do {

      /*
       * Try to balance load through swapping.
       */
      LoadType swap_transfer = adjustLoadBySwapping(
         main_bin,
         hold_bin,
         ideal_load,
         low_load,
         high_load);

      actual_transfer += swap_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            main_bin,
            hold_bin,
            (main_bin.getSumLoad() - ideal_load));
         tbox::plog << "  Balance penalty after adjustLoadBySwapping = "
                    << balance_penalty
                    << ", needs " << (ideal_load-main_bin.getSumLoad())
                    << " more with " << main_bin.size() << " main_bin and "
                    << hold_bin.size() << " hold_bin Boxes remaining."
                    << "\n  main_bin now has " << main_bin.getSumLoad()
                    << " in " << main_bin.size() << " boxes."
                    << std::endl;
      }

      // Skip breaking if already in range.
      if (main_bin.getSumLoad() <= high_load && main_bin.getSumLoad() >= low_load ) break;

      /*
       * Skip breaking if adding/subtracting d_min_load overshoots the range and worsens distance to range.
       */
      if ( tbox::MathUtilities<double>::Abs(main_bin.getSumLoad() - 0.5*(high_load+low_load)) <= 0.5*d_min_load ) {
         break;
      }

      /*
       * Assuming that we did the best we could, swapping
       * some BoxInTransit without breaking any, we now break up a Box
       * in the overloaded side for partial transfer to the
       * underloaded side.
       */
      LoadType brk_transfer = -adjustLoadByBreaking(
         main_bin,
         hold_bin,
         next_available_index,
         ideal_load,
         low_load,
         high_load );
      actual_transfer += brk_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            main_bin,
            hold_bin,
            (main_bin.getSumLoad() - ideal_load));
         tbox::plog << "  Balance penalty after adjustLoadByBreaking = "
                    << balance_penalty
                    << ", needs " << (ideal_load-main_bin.getSumLoad())
                    << " more with " << main_bin.size() << " main_bin and "
                    << hold_bin.size() << " hold_bin Boxes remaining."
                    << "\n  main_bin now has " << main_bin.getSumLoad()
                    << " in " << main_bin.size() << " boxes."
                    << std::endl;
      }
      if (brk_transfer == 0) {
         /*
          * If no box can be broken to improve the actual_transfer,
          * there is nothing further we can do.  The swap phase, tried
          * before the break phase, also generated no transfer, so
          * there's no point trying again.  Break out now to save
          * retrying the swap phase.
          */
         if (d_print_steps) {
            tbox::plog << "  adjustLoad stopping due to unsuccessful break."
                       << std::endl;
         }
         break;
      }

      /*
       * Now that we have broken up a Box, redo this loop to
       * see if swapping can produce a better result.
       */
   } while ( ( main_bin.getSumLoad() >= high_load ) ||
             ( main_bin.getSumLoad() <= low_load ) );

   if ( d_print_steps ) {
      const LoadType point_miss = main_bin.getSumLoad() - ideal_load;
      const LoadType range_miss =
         main_bin.getSumLoad() > high_load ? main_bin.getSumLoad() - high_load :
         main_bin.getSumLoad() < low_load ? low_load - main_bin.getSumLoad() : 0;
      tbox::plog << "  adjustLoad point_miss=" << point_miss
                 << "  range_miss="
                 << (range_miss > 0 ? " ":"") // Add space if missed range
                 << (range_miss > 0.5*d_min_size.getProduct() ? " ":"") // Add space if missed range by a lot
                 << range_miss
                 << "  " << main_bin.getSumLoad() << '/'
                 << ideal_load << " [" << low_load << ',' << high_load << ']'
                 << std::endl;
   }

   t_adjust_load->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::packSubtreeDataUp(
   tbox::MessageStream& msg,
   const SubtreeData& subtree_data) const
{
   t_pack_load->start();
   msg << subtree_data.d_num_procs;
   msg << subtree_data.d_subtree_load_current;
   msg << subtree_data.d_subtree_load_ideal;
   msg << subtree_data.d_subtree_load_upperlimit;
   msg << subtree_data.d_eff_num_procs;
   msg << subtree_data.d_eff_load_current;
   msg << subtree_data.d_eff_load_ideal;
   msg << subtree_data.d_eff_load_upperlimit;
   const TransitSet& for_export = subtree_data.d_work_traded;
   msg << static_cast<int>(for_export.size());
   for (TransitSet::const_iterator
        ni = for_export.begin(); ni != for_export.end(); ++ni) {
      const BoxInTransit& box_in_transit = *ni;
      box_in_transit.putToMessageStream(msg);
   }
   msg << subtree_data.d_wants_work_from_parent;
   t_pack_load->stop();
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::unpackSubtreeDataUp(
   SubtreeData& subtree_data,
   hier::LocalId& next_available_index,
   tbox::MessageStream &msg ) const
{
   t_unpack_load->start();
   int num_boxes = 0;
   msg >> subtree_data.d_num_procs;
   msg >> subtree_data.d_subtree_load_current;
   msg >> subtree_data.d_subtree_load_ideal;
   msg >> subtree_data.d_subtree_load_upperlimit;
   msg >> subtree_data.d_eff_num_procs;
   msg >> subtree_data.d_eff_load_current;
   msg >> subtree_data.d_eff_load_ideal;
   msg >> subtree_data.d_eff_load_upperlimit;
   msg >> num_boxes;
   /*
    * As we pull each BoxInTransit out, give it a new id that reflects
    * its new owner.
    */
   BoxInTransit received_box(d_dim);
   for (int i = 0; i < num_boxes; ++i) {
      received_box.getFromMessageStream(msg);
      BoxInTransit renamed_box(received_box,
                               received_box.getBox(),
                               d_mpi.getRank(),
                               next_available_index);
      next_available_index += 2 + d_rank_tree->getDegree();
      subtree_data.d_work_traded.insert(renamed_box);
   }
   msg >> subtree_data.d_wants_work_from_parent;
   t_unpack_load->stop();
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::packSubtreeDataDown(
   tbox::MessageStream& msg,
   const SubtreeData& subtree_data) const
{
   t_pack_load->start();
   const TransitSet& for_export = subtree_data.d_work_traded;
   msg << static_cast<int>(for_export.size());
   for (TransitSet::const_iterator
        ni = for_export.begin(); ni != for_export.end(); ++ni) {
      const BoxInTransit& box_in_transit = *ni;
      box_in_transit.putToMessageStream(msg);
   }
   t_pack_load->stop();
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::unpackSubtreeDataDown(
   SubtreeData& subtree_data,
   hier::LocalId& next_available_index,
   tbox::MessageStream &msg ) const
{
   t_unpack_load->start();
   int num_boxes = 0;
   msg >> num_boxes;
   /*
    * As we pull each BoxInTransit out, give it a new id that reflects
    * its new owner.
    */
   BoxInTransit received_box(d_dim);
   for (int i = 0; i < num_boxes; ++i) {
      received_box.getFromMessageStream(msg);
      BoxInTransit renamed_box(received_box,
                               received_box.getBox(),
                               d_mpi.getRank(),
                               next_available_index);
      next_available_index += 2 + d_rank_tree->getDegree();
      subtree_data.d_work_traded.insert(renamed_box);
   }
   t_unpack_load->stop();
}



/*
 *************************************************************************
 * Compute the surplus per descendent still waiting for load from parent.
 * This surplus, if any, is the difference between the available
 * load_for_descendents, spread over those descendents.
 *
 * Surplus per descendent will be zero if we don't need to push surplus
 * work to children.
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::computeSurplusPerEffectiveDescendent(
   const TransitSet &unassigned,
   const LoadType group_avg_load,
   const std::vector<SubtreeData> &child_subtrees,
   int first_child ) const
{
   const LoadType load_for_me = group_avg_load*(1 + d_flexible_load_tol);

   // Available amount for descendents after removing load_for_me:
   const LoadType load_for_descendents = unassigned.getSumLoad() - load_for_me;

   // Total of ideal exports to children:
   LoadType ideal_export_to_children = 0.0;
   int num_effective_des = 0;
   for (size_t ichild = first_child; ichild < child_subtrees.size(); ++ichild) {
      if ( child_subtrees[ichild].d_wants_work_from_parent ) {
         ideal_export_to_children += child_subtrees[ichild].effDeficit();
         num_effective_des += child_subtrees[ichild].d_eff_num_procs;
      }
   }

   // Amount of surplus per effective descendent.
   const LoadType surplus_per_effective_descendent =
      (load_for_descendents - ideal_export_to_children) / num_effective_des;

   if ( d_print_steps ) {
      tbox::plog << "load_for_me = " << load_for_me
                 << ",  load_for_descendents = " << load_for_descendents
                 << ",  num_effective_des = " << num_effective_des
                 << ",  ideal_export_to_children = " << ideal_export_to_children
                 << ",  surplus_per_effective_descendent " << surplus_per_effective_descendent
                 << std::endl;
   }

   return surplus_per_effective_descendent;
}





/*
 *************************************************************************
 * Construct semilocal relationships in unbalanced--->balanced
 * Connector.
 *
 * Determine edges in unbalanced_to_balanced by sending balanced
 * BoxInTransit back to the owners of the unbalanced Boxes that
 * originated them.  We don't know what ranks will send back the
 * balanced boxes, so we keep receiving messages from any rank until
 * we have accounted for all the cells in the unbalanced BoxLevel.
 *************************************************************************
 */
void
TreeLoadBalancer::constructSemilocalUnbalancedToBalanced(
   hier::MappingConnector &unbalanced_to_balanced,
   const TreeLoadBalancer::TransitSet &kept_imports ) const
{
   t_construct_semilocal->start();

   // Stuff the imported BoxInTransits into buffers by their original owners.
   t_pack_edge->start();
   std::map<int,boost::shared_ptr<tbox::MessageStream> > outgoing_messages;
   for ( TransitSet::const_iterator bi=kept_imports.begin();
         bi!=kept_imports.end(); ++bi ) {
      const BoxInTransit &bit = *bi;
      boost::shared_ptr<tbox::MessageStream> &mstream = outgoing_messages[bit.d_orig_box.getOwnerRank()];
      if ( !mstream ) {
         mstream.reset(new tbox::MessageStream);
      }
      bit.putToMessageStream(*mstream);
   }
   t_pack_edge->stop();


   /*
    * Send outgoing_messages.  Optimization for mitigating contention:
    * Start by sending to the first recipient with a rank higher than
    * the local rank.
    */

   std::map<int,boost::shared_ptr<tbox::MessageStream> >::iterator recip_itr =
      outgoing_messages.upper_bound(d_mpi.getRank());
   if ( recip_itr == outgoing_messages.end() ) {
      recip_itr = outgoing_messages.begin();
   }

   int outgoing_messages_size = static_cast<int>(outgoing_messages.size());
   std::vector<tbox::SAMRAI_MPI::Request>
      send_requests( outgoing_messages_size, MPI_REQUEST_NULL );

   t_post_load_distribution_barrier->start();
   d_mpi.Barrier(); // This barrier seems to speed up the load balancing, maybe by allowing one communication phase to finish before beginning another.
   t_post_load_distribution_barrier->stop();

   for ( int send_number = 0; send_number < outgoing_messages_size; ++send_number ) {

      int recipient = recip_itr->first;
      tbox::MessageStream &mstream = *recip_itr->second;

      if ( d_print_edge_steps ) {
         tbox::plog << "Accounting for cells on proc " << recipient << '\n';
      }

      d_mpi.Isend(
         (void*)(mstream.getBufferStart()),
         static_cast<int>(mstream.getCurrentSize()),
         MPI_CHAR,
         recipient,
         TreeLoadBalancer_EDGETAG0,
         &send_requests[send_number]);

      ++recip_itr;
      if ( recip_itr == outgoing_messages.end() ) {
         recip_itr = outgoing_messages.begin();
      }

   }


   /*
    * Determine number of cells in unbalanced that are not yet accounted
    * for in balanced.
    */
   int num_unaccounted_cells = static_cast<int>(
      unbalanced_to_balanced.getBase().getLocalNumberOfCells());
   if ( d_print_edge_steps ) {
      tbox::plog << num_unaccounted_cells << " unbalanced cells\n";
   }

   const hier::BoxContainer &unbalanced_boxes = unbalanced_to_balanced.getBase().getBoxes();
   for ( hier::BoxContainer::const_iterator bi=unbalanced_boxes.begin();
         bi!=unbalanced_boxes.end(); ++bi ) {

      const hier::Box &unbalanced_box = *bi;

      hier::Connector::ConstNeighborhoodIterator neighborhood_itr =
         unbalanced_to_balanced.findLocal(unbalanced_box.getBoxId());

      if ( neighborhood_itr != unbalanced_to_balanced.end() ) {
         // unbalanced_box has neighborhood.
         for ( hier::Connector::ConstNeighborIterator ni=unbalanced_to_balanced.begin(neighborhood_itr);
               ni!=unbalanced_to_balanced.end(neighborhood_itr); ++ni ) {
            TBOX_ASSERT( ni->getOwnerRank() == d_mpi.getRank() );
            num_unaccounted_cells -= ni->size();
         }

      }
      else {
         // unbalanced_box has no neighborhood.
         num_unaccounted_cells -= unbalanced_box.size();
      }

   }
   if ( d_print_edge_steps ) {
      tbox::plog << num_unaccounted_cells << " unaccounted cells\n";
   }


   /*
    * Receive info about exported cells from processes that now own
    * those cells.  Receive until all cells are accounted for.
    */

   std::vector<char> incoming_message; // Keep outside loop to avoid reconstructions.
   BoxInTransit balanced_box_in_transit(d_dim);
   while ( num_unaccounted_cells > 0 ) {

      t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Status status;
      d_mpi.Probe( MPI_ANY_SOURCE,
                   TreeLoadBalancer_EDGETAG0,
                   &status );

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count( &status, MPI_CHAR, &count );
      incoming_message.resize( count, -1 );

      d_mpi.Recv(
         static_cast<void*>(&incoming_message[0]),
         count,
         MPI_CHAR,
         source,
         TreeLoadBalancer_EDGETAG0,
         &status );
      t_construct_semilocal_comm_wait->stop();

      tbox::MessageStream msg( incoming_message.size(),
                               tbox::MessageStream::Read,
                               static_cast<void*>(&incoming_message[0]),
                               false );
      const int old_count = num_unaccounted_cells;
      t_unpack_edge->start();
      while ( !msg.endOfData() ) {

         balanced_box_in_transit.getFromMessageStream(msg);
         TBOX_ASSERT( balanced_box_in_transit.d_orig_box.getOwnerRank() == d_mpi.getRank() );
         unbalanced_to_balanced.insertLocalNeighbor(
            balanced_box_in_transit.d_box,
            balanced_box_in_transit.d_orig_box.getBoxId() );
         num_unaccounted_cells -= balanced_box_in_transit.d_box.size();

      }
      t_unpack_edge->stop();

      if ( d_print_edge_steps ) {
         tbox::plog << "Process " << source << " accounted for "
                    << (old_count-num_unaccounted_cells) << " cells, leaving "
                    << num_unaccounted_cells << " unaccounted.\n";
      }

   }
   TBOX_ASSERT( num_unaccounted_cells == 0 );
   incoming_message.clear();


   // Wait for the sends to complete before clearing outgoing_messages.
   if (send_requests.size() > 0) {
      std::vector<tbox::SAMRAI_MPI::Status> status(send_requests.size());
      t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Waitall(
         static_cast<int>(send_requests.size()),
         &send_requests[0],
         &status[0]);
      t_construct_semilocal_comm_wait->stop();
      outgoing_messages.clear();
   }
   else {
      t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Waitall(0, NULL, NULL);
      t_construct_semilocal_comm_wait->stop();
   }

   t_construct_semilocal->stop();

   return;
}



/*
 *************************************************************************
 * Set the MPI commuicator.  If there's a private communicator, free
 * it first.  It's safe to free the private communicator because no
 * other code have access to it.
 *************************************************************************
 */
void
TreeLoadBalancer::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& samrai_mpi)
{
   if (samrai_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      TBOX_ERROR("TreeLoadBalancer::setSAMRAI_MPI error: Given\n"
         << "communicator is invalid.");
   }

   if ( d_mpi_is_dupe ) {
      d_mpi.freeCommunicator();
   }

   // Enable private communicator.
   d_mpi.dupCommunicator(samrai_mpi);
   d_mpi_is_dupe = true;
}



/*
 *************************************************************************
 * Set the MPI commuicator.
 *************************************************************************
 */
void
TreeLoadBalancer::freeMPICommunicator()
{
   if ( d_mpi_is_dupe ) {
      // Free the private communicator (if MPI has not been finalized).
      int flag;
      tbox::SAMRAI_MPI::Finalized(&flag);
      if (!flag) {
         d_mpi.freeCommunicator();
      }
   }
   d_mpi.setCommunicator(tbox::SAMRAI_MPI::commNull);
   d_mpi_is_dupe = false;
}



/*
 *************************************************************************
 * RankGroups load-balances within their membership and ignore other
 * groups.  When we balance over multiple cycles, the RankGroup for
 * the local process depends on the cycle, as computed by this method.
 *
 * The RankGroup size increases exponentially with the cycle number
 * such that for the last cycle the rank group includes all processes
 * in d_mpi.  It's a heuristic formula, as follows:
 *
 * Partition all ranks into similar sized groups.  With each cycle,
 * the group size grows exponentially while the number of groups
 * shrinks.  The last cycle_number has a single group of
 * d_mpi.getSize() processors.
 *
 * Let m = number of cycles
 * i = cycle number, [0,m)
 * p = communicator size
 *
 * The group size is p^((i+1)/m)
 *************************************************************************
 */
void
TreeLoadBalancer::createBalanceRankGroupBasedOnCycles(
   tbox::RankGroup &rank_group,
   int &number_of_groups,
   int &group_num,
   const int cycle_number,
   const int number_of_cycles) const
{

   /*
    * Compute the number of group and, implicitly, the group sizes.
    * Tiny groups tend to leave the members with possibly large
    * overloads.  In order to make all groups similar in size we round
    * down the number of groups (and round up the group size).
    */
   number_of_groups =
      static_cast<int>(pow(static_cast<double>(d_mpi.getSize()),
                           1.0 - double(cycle_number + 1) / number_of_cycles));

   /*
    * All groups will have a base population count of
    * d_mpi.getSize()/number_of_groups.  The remainder from the
    * integer division is distributed to a subset of groups, starting
    * from group 0, so these groups will have one more than the base.
    */
   const int base_group_size = d_mpi.getSize() / number_of_groups;
   const int first_base_sized_group = d_mpi.getSize() % number_of_groups;
   const int first_rank_in_base_sized_group =
      first_base_sized_group * (1 + base_group_size);

   if (d_mpi.getRank() < first_rank_in_base_sized_group) {
      group_num = d_mpi.getRank() / (1 + base_group_size);
      const int group_first_rank = group_num * (1 +base_group_size);
      rank_group.setMinMax( group_first_rank,
                            group_first_rank + base_group_size );
   } else {
      group_num = first_base_sized_group
         + (d_mpi.getRank() - first_rank_in_base_sized_group) / base_group_size;
      const int group_first_rank = first_rank_in_base_sized_group +
         (group_num - first_base_sized_group)*(1+base_group_size);
      rank_group.setMinMax( group_first_rank,
                            group_first_rank + base_group_size - 1 );
   }

   return;
}



/*
 *************************************************************************
 * Set up the asynchronous communication objects for the process tree
 * containing ranks defined by the RankGroup.
 *
 * The process tree lay-out is defined by the BalancedDepthFirstTree
 * class, thus defining parent and children of the local process.
 * This method sets the AsyncCommPeer objects for communication with
 * children and parent.
 *************************************************************************
 */
void
TreeLoadBalancer::setupAsyncCommObjects(
   tbox::AsyncCommStage& child_stage,
   tbox::AsyncCommPeer<char> *& child_comms,
   tbox::AsyncCommStage& parent_stage,
   tbox::AsyncCommPeer<char> *& parent_comm,
   const tbox::RankGroup &rank_group ) const
{

   child_comms = parent_comm = 0;

   const int num_children = d_rank_tree->getNumberOfChildren();

   if ( num_children > 0 ) {

      child_comms = new tbox::AsyncCommPeer<char>[num_children];

      for (int child_num = 0; child_num < num_children; ++child_num) {

         const int child_rank_in_grp = d_rank_tree->getChildRank(child_num);
         const int child_true_rank = rank_group.getMappedRank(child_rank_in_grp);

         child_comms[child_num].initialize(&child_stage);
         child_comms[child_num].setPeerRank(child_true_rank);
         child_comms[child_num].setMPI(d_mpi);
         child_comms[child_num].setMPITag(TreeLoadBalancer_LOADTAG0,
                                          TreeLoadBalancer_LOADTAG1);
         child_comms[child_num].limitFirstDataLength(
            sizeof(BoxInTransit)*TreeLoadBalancer_FIRSTDATALEN);
      }
   }

   if (d_rank_tree->getParentRank() != tbox::RankTreeStrategy::getInvalidRank()) {

      const int parent_rank_in_grp = d_rank_tree->getParentRank();
      int parent_true_rank = rank_group.getMappedRank(parent_rank_in_grp);

      parent_comm = new tbox::AsyncCommPeer<char>;
      parent_comm->initialize(&parent_stage);
      parent_comm->setPeerRank(parent_true_rank);
      parent_comm->setMPI(d_mpi);
      parent_comm->setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      parent_comm->limitFirstDataLength(
         sizeof(BoxInTransit)*TreeLoadBalancer_FIRSTDATALEN);

   }

   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::destroyAsyncCommObjects(
   tbox::AsyncCommPeer<char> *& child_comms,
   tbox::AsyncCommPeer<char> *& parent_comm) const
{
   if (d_mpi.getSize() == 1) {
      TBOX_ASSERT(child_comms == 0);
      TBOX_ASSERT(parent_comm == 0);
   } else {
      if ( child_comms != 0 ) {
         delete[] child_comms;
      }
      if ( parent_comm != 0 ) {
         delete parent_comm;
      }
      child_comms = parent_comm = 0;
   }
}



/*
 *************************************************************************
 * Attempt bring main_bin to within a specific load range by moving
 * one box to/from it from/to hold_bin.  This method is allowed to break
 * a box and move parts of it.
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::adjustLoadByBreaking(
   TransitSet& main_bin,
   TransitSet& hold_bin,
   hier::LocalId& next_available_index,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load ) const
{
   LoadType actual_transfer = 0;

   if (main_bin.getSumLoad() < low_load) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      actual_transfer = -adjustLoadByBreaking(
         hold_bin,
         main_bin,
         next_available_index,
         hold_bin.getSumLoad()-(ideal_load-main_bin.getSumLoad()),
         hold_bin.getSumLoad()-(high_load-main_bin.getSumLoad()),
         hold_bin.getSumLoad()-(low_load-main_bin.getSumLoad()) );
      return actual_transfer;
   }

   TBOX_ASSERT(low_load <= ideal_load);
   TBOX_ASSERT(ideal_load <= high_load);
   TBOX_ASSERT(main_bin.getSumLoad() > high_load);

   TBOX_ASSERT(main_bin.size() + hold_bin.size() > 0);

   t_shift_loads_by_breaking->start();

   const LoadType ideal_transfer = main_bin.getSumLoad() - ideal_load;
   const LoadType high_transfer = main_bin.getSumLoad() - low_load;
   const LoadType low_transfer = main_bin.getSumLoad() - high_load;

   if (d_print_steps) {
      tbox::plog << "    adjustLoadByBreaking asked to break off "
                 << ideal_transfer << " [" << low_transfer << ','
                 << high_transfer << "] from one of " << main_bin.size()
                 << " Boxes to add to set of " << hold_bin.size()
                 << " Boxes."
                 << std::endl;
   }


   // Data for the best cutting results so far:
   std::vector<hier::Box> breakoff;
   std::vector<hier::Box> leftover;
   double breakoff_amt = 0.0;
   BoxInTransit breakbox(d_dim);

   int break_acceptance_flags[3] = {0,0,0};
   int &found_breakage = break_acceptance_flags[2];

   for (TransitSet::iterator si = main_bin.begin(); si != main_bin.end(); ++si) {

      const BoxInTransit& candidate = *si;

      if (d_print_steps) {
         tbox::plog << "    Considering break candidate " << candidate
                    << std::endl;
      }

      std::vector<hier::Box> trial_breakoff;
      std::vector<hier::Box> trial_leftover;
      double trial_breakoff_amt;

      breakOffLoad(
         trial_breakoff,
         trial_leftover,
         trial_breakoff_amt,
         candidate.d_box,
         ideal_transfer,
         low_transfer,
         high_transfer );

      if (!trial_breakoff.empty()) {

         const bool accept_break = evaluateBreak(
            break_acceptance_flags, breakoff_amt, trial_breakoff_amt,
            ideal_transfer, low_transfer, high_transfer );
         if (d_print_break_steps) {
            tbox::plog << "      Break evaluation:"
                       << "  " << break_acceptance_flags[0]
                       << "  " << break_acceptance_flags[1]
                       << "  " << break_acceptance_flags[2]
                       << std::endl;
         }

         if (d_print_break_steps) {
            tbox::plog << "    Potential to replace " << candidate << " with "
                       << trial_breakoff.size() << " breakoff Boxes and "
                       << trial_leftover.size() << " leftover Boxes."
                       << "  break amount = " << trial_breakoff_amt
                       << "  in-range imp = " << break_acceptance_flags[0]
                       << "  balance imp = " << break_acceptance_flags[1]
                       << "  overal imp = " << break_acceptance_flags[2]
                       << "  accept_break = " << accept_break
                       << std::endl;
         }

         if (accept_break) {
            breakbox = candidate;
            breakoff_amt = trial_breakoff_amt;
            breakoff.swap(trial_breakoff);
            leftover.swap(trial_leftover);
            if ( break_acceptance_flags[0] == 1 ) {
               // We are in the [low,high] range.  That is sufficient.
               break;
            }
         }

      } else {
         if (d_print_break_steps) {
            tbox::plog << "    Break step could not break " << ideal_transfer
                       << " from main_bin box " << candidate
                       << std::endl;
         }
      }

   }


   if ( found_breakage == 1 ) {
      /*
       * Remove the chosen candidate.  Put its breakoff parts
       * in hold_bin and its leftover parts back into main_bin.
       */
      main_bin.erase(breakbox);
      for (std::vector<hier::Box>::const_iterator bi = breakoff.begin();
           bi != breakoff.end();
           ++bi) {
         BoxInTransit give_box_in_transit(
            breakbox,
            *bi,
            d_mpi.getRank(),
            next_available_index);
         give_box_in_transit.d_boxload = static_cast<int>(computeLoad(
                                                             give_box_in_transit.d_orig_box,
                                                             give_box_in_transit.getBox()));
         next_available_index += 2 + d_rank_tree->getDegree();
         hold_bin.insert(give_box_in_transit);
         actual_transfer += give_box_in_transit.d_boxload;
         if (d_print_break_steps) {
            tbox::plog << "    Breakoff box " << *bi << bi->numberCells()
                       << '|' << bi->size()
                       << " -> " << give_box_in_transit
                       << std::endl;
         }
      }
      for (std::vector<hier::Box>::const_iterator bi = leftover.begin();
           bi != leftover.end();
           ++bi) {
         BoxInTransit keep_box_in_transit(
            breakbox,
            *bi,
            d_mpi.getRank(),
            next_available_index);
         keep_box_in_transit.d_boxload = static_cast<int>(computeLoad(
                                                             keep_box_in_transit.d_orig_box,
                                                             keep_box_in_transit.getBox()));
         next_available_index += 2 + d_rank_tree->getDegree();
         main_bin.insert(keep_box_in_transit);
         if (d_print_break_steps) {
            tbox::plog << "    Leftover box " << *bi << bi->numberCells()
                       << '|' << bi->size()
                       << " -> " << keep_box_in_transit
                       << std::endl;
         }
      }
   }

   t_shift_loads_by_breaking->stop();
   return actual_transfer;
}



/*
 *************************************************************************
 * Attempt to adjust the load of a main_bin by swapping boxes with
 * a hold_bin.
 *
 * Transfering a BoxInTransit from one TransitSet to another
 * is considered a degenerate "swap" (a BoxInTransit is
 * swapped for nothing) handled by this function.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of main_bin
 * ideal_transfer < 0 means to raise the load of hold_bin
 * The iterative do loop may overshoot the ideal_transfer
 * and may have to swap to shift some of the load
 * back.
 *
 * Return whether any changes were made.
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::adjustLoadBySwapping(
   TransitSet& main_bin,
   TransitSet& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load ) const
{
   TBOX_ASSERT( high_load >= ideal_load );
   TBOX_ASSERT( low_load <= ideal_load );

   t_adjust_load_by_swapping->start();

   if (d_print_steps) {
      tbox::plog << "  Attempting to bring main_bin from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " [" << low_load << ',' << high_load
                 << "] by swapping."
                 << std::endl;
   }

   bool found_swap;

   LoadType actual_transfer = 0;

   do {

      /*
       * Ammount we seek to transfer from hi to lo
       * (the "ideal" for this particular iteration).
       * Unlike ideal_transfer and actual_transfer, this quantity is positive.
       */
      LoadType rem_transfer = main_bin.getSumLoad() - ideal_load;
      LoadType low_transfer = main_bin.getSumLoad() - high_load;
      LoadType high_transfer = main_bin.getSumLoad() - low_load;
      if (d_print_swap_steps) {
         tbox::plog << "    Swap progress: " << main_bin.getSumLoad()
                    << " / " << ideal_load << " remaining transfer = "
                    << rem_transfer << " [" << low_transfer << ','
                    << high_transfer << ']' << std::endl;
      }

      found_swap = false;

      LoadType swap_transfer;
      found_swap = swapLoadPair(
         main_bin,
         hold_bin,
         swap_transfer,
         rem_transfer,
         low_transfer,
         high_transfer);
      swap_transfer = -swap_transfer;

      if (found_swap) {
         actual_transfer += swap_transfer;
      }

   } while (found_swap &&
            (main_bin.getSumLoad() < low_load || main_bin.getSumLoad() > high_load ));

   if (d_print_swap_steps) {
      tbox::plog << "  Final balance for adjustLoadBySwapping: "
                 << main_bin.getSumLoad() << " / " << ideal_load
                 << "  Off by " << (main_bin.getSumLoad()-ideal_load)
                 << std::endl;
   }

   t_adjust_load_by_swapping->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 * Find a BoxInTransit in src and a BoxInTransit in dst which when
 * swapped results in shifting close to ideal_shift from src to dst.
 * Make the swap.  Return whether a swap pair was found.
 *************************************************************************
 */
bool
TreeLoadBalancer::swapLoadPair(
   TransitSet& src,
   TransitSet& dst,
   LoadType& actual_transfer,
   LoadType ideal_transfer,
   LoadType low_transfer,
   LoadType high_transfer ) const
{
   if (ideal_transfer < 0) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      bool rval = swapLoadPair(
         dst,
         src,
         actual_transfer,
         -ideal_transfer,
         -high_transfer,
         -low_transfer);
      actual_transfer = -actual_transfer;
      return rval;
   }

   t_find_swap_pair->start();

   if (d_print_swap_steps) {
      tbox::plog << "    swapLoadPair looking for transfer of "
                 << ideal_transfer
                 << " between " << src.size() << "-box src and "
                 << dst.size() << "-box dst." << std::endl;
      tbox::plog << "      src (" << src.size() << "):" << std::endl;
      for (TransitSet::iterator si = src.begin(); si != src.end(); ++si) {
         tbox::plog << "        " << *si << std::endl;
      }
      tbox::plog << "      dst (" << dst.size() << "):" << std::endl;
      for (TransitSet::iterator si = dst.begin(); si != dst.end(); ++si) {
         tbox::plog << "        " << *si << std::endl;
      }
   }

   /*
    * Look for two swap options.  The "high side" option would
    * transfer at least ideal_transfer.  The "low side" option would
    * transfer up to ideal_transfer.
    *
    * Each option is defined by a box from src and a box from dst,
    * designated by the iterators src_hiside, dst_hiside, src_loside
    * and dst_loside.  src_hiside points to the box in the src for the
    * high-side transfer, and similarly for dst_hiside.  src_loside
    * points to the box in the src for the low-side transfer, and
    * similarly for dst_loside.
    *
    * Note that in the degenerate case, the dst box does not exist,
    * and the swap degenerates to moving a box from the src to the
    * dst.
    *
    * Compute the balance_penalty if high and low were swapped.  Keep
    * looking until we find the pair giving the lowest balance_penalty
    * on swapping.
    *
    * isrc and idst point to the current best pair to swap.  new_balance_penalty
    * is the balance_penalty if we swap them.
    *
    * src_test and dst_test are trial pairs to check to see if we can improve on
    * new_balance_penalty.
    *
    * We will look for two "best" pairs:
    *
    * TODO: This method was originally written to compute the best
    * hiside and loside options separately and compare them at the
    * end.  That separation may not be needded anymore.  It may be
    * possible to simplify this method by keeping only the best option
    * at any time.
    */

   // Initialization indicating no swap pair found yet.
   TransitSet::iterator src_hiside = src.end();
   TransitSet::iterator dst_hiside = dst.end();
   TransitSet::iterator src_loside = src.end();
   TransitSet::iterator dst_loside = dst.end();

   // A dummy BoxInTransit for set searches.
   hier::Box dummy_box(d_dim);
   BoxInTransit dummy_search_target(d_dim);

   // Difference between swap results and ideal, >= 0
   LoadType hiside_transfer = 0.0;
   LoadType loside_transfer = 0.0;


   int loside_acceptance_flags[3] = {0,0,0};
   int hiside_acceptance_flags[3] = {0,0,0};

   if (dst.empty()) {
      /*
       * There is no dst BoxInTransit, so the swap would
       * degnerate to moving a box from src to dst.  Find
       * the best src BoxInTransit to move.
       */
      dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
      dummy_search_target.d_boxload = ideal_transfer;
      const TransitSet::iterator src_test = src.lower_bound(dummy_search_target);

      if (d_print_swap_steps) {
         tbox::plog << "  swapLoadPair with empty dst: ";
      }

      if (src_test != src.begin()) {
         TransitSet::iterator src_test1 = src_test;
         --src_test1;
         if ( evaluateBreak( hiside_acceptance_flags, hiside_transfer, src_test1->d_boxload,
                             ideal_transfer, low_transfer, high_transfer ) ) {
            src_hiside = src_test1;
            hiside_transfer = src_hiside->d_boxload;
            if (d_print_swap_steps) {
               tbox::plog << "  hi src: " << (*src_hiside)
                          << " with transfer " << src_hiside->d_boxload
                          << ", off by " << hiside_transfer-ideal_transfer
                          << ", acceptance_flags=" << hiside_acceptance_flags[0]
                          << ',' << hiside_acceptance_flags[1]
                          << ',' << hiside_acceptance_flags[2];
            }
         }
      }
      if (src_test != src.end()) {
         if ( evaluateBreak( loside_acceptance_flags, loside_transfer, src_test->d_boxload,
                             ideal_transfer, low_transfer, high_transfer ) ) {
            src_loside = src_test;
            loside_transfer = src_loside->d_boxload;
            if (d_print_swap_steps) {
               tbox::plog << "  lo src: " << (*src_loside)
                          << " with transfer " << src_loside->d_boxload
                          << ", off by " << loside_transfer-ideal_transfer
                          << ", acceptance_flags=" << loside_acceptance_flags[0]
                          << ',' << loside_acceptance_flags[1]
                          << ',' << loside_acceptance_flags[2];
            }
         }
      }
      if (d_print_swap_steps) {
         tbox::plog << std::endl;
      }

   } else {

      /*
       * Start search through src beginning with the box whose load
       * exceeds the biggest dst box by at least ideal_transfer.
       */
      dummy_search_target = *dst.begin();
      dummy_search_target.d_boxload += ideal_transfer;
      TransitSet::iterator src_beg = src.lower_bound(dummy_search_target);

      for (TransitSet::iterator src_test = src_beg; src_test != src.end(); ++src_test) {

         /*
          * Set dst_test pointing to where we should start looking in dst.
          * Look for a load less than the load of src_test by
          * ideal_transfer.
          */
         dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
         dummy_search_target.d_boxload = tbox::MathUtilities<LoadType>::Max(
               src_test->d_boxload - ideal_transfer,
               0);
         TransitSet::iterator dst_test = dst.lower_bound(dummy_search_target);

         if (dst_test != dst.end()) {

            /*
             * lower_bound returned dst_test that would transfer >=
             * ideal_transfer when swapped with src_test.  Check
             * transfererence between src_test and dst_test for the
             * high-side transfer.  Also check the next smaller box in
             * dst for the low-side transfer.
             */

            evaluateBreak( hiside_acceptance_flags, hiside_transfer,
                           src_test->d_boxload - dst_test->d_boxload,
                           ideal_transfer, low_transfer, high_transfer );

            if ( hiside_acceptance_flags[2] == 1 ) {
               src_hiside = src_test;
               dst_hiside = dst_test;
               hiside_transfer = src_hiside->d_boxload - dst_hiside->d_boxload;
               if (d_print_swap_steps) {
                  tbox::plog << "    new hi-swap pair: " << (*src_hiside)
                             << " & " << (*dst_hiside) << " with transfer "
                             << hiside_transfer
                             << " missing by " << hiside_transfer-ideal_transfer
                             << std::endl;
               }
            }

            if (dst_test != dst.begin()) {
               --dst_test; // Now, src_test and dst_test transferer by *less* than ideal_transfer.

               evaluateBreak( loside_acceptance_flags, loside_transfer,
                              src_test->d_boxload - dst_test->d_boxload,
                              ideal_transfer, low_transfer, high_transfer );

               if ( loside_acceptance_flags[2] == 1 ) {
                  src_loside = src_test;
                  dst_loside = dst_test;
                  loside_transfer = src_loside->d_boxload - dst_loside->d_boxload;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap pair: " << (*src_loside)
                                << " & " << (*dst_loside) << " with transfer "
                                << loside_transfer
                                << " missing by " << loside_transfer-ideal_transfer
                                << std::endl;
                  }
               }
            }

         } else {

            /*
             * The ideal dst to swap is smaller than the smallest dst
             * box.  So the only choice is swapping src_test for nothing.
             * Chech this against the current high- and low-side choices.
             */
            if (src_test->d_boxload > ideal_transfer) {
               // Moving src_test to src is moving too much--hiside.

               evaluateBreak( hiside_acceptance_flags, hiside_transfer,
                              src_test->d_boxload,
                              ideal_transfer, low_transfer, high_transfer );

               if ( hiside_acceptance_flags[2] == 1 ) {
                  src_hiside = src_test;
                  dst_hiside = dst.end();
                  hiside_transfer = src_hiside->d_boxload;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new hi-swap source: " << (*src_hiside)
                                << " & " << "no dst" << " with transfer "
                                << (src_hiside->d_boxload)
                                << " missing by " << hiside_transfer-ideal_transfer
                                << std::endl;
                  }
               }
            } else {
               // Moving src_test to src is moving (just right or) too little--loside.

               evaluateBreak( loside_acceptance_flags, loside_transfer,
                              src_test->d_boxload,
                              ideal_transfer, low_transfer, high_transfer );

               if ( loside_acceptance_flags[2] == 1 ) {
                  src_loside = src_test;
                  dst_loside = dst.end();
                  loside_transfer = src_loside->d_boxload;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap source: " << (*src_loside)
                                << " & " << "no dst" << " with transfer "
                                << (src_loside->d_boxload)
                                << " missing by " << loside_transfer-ideal_transfer
                                << std::endl;
                  }
               }
               /*
                * Break out of the loop early because there is no
                * point checking smaller src boxes.
                */
               break;
            }
         }

         if ( ( low_transfer <= loside_transfer && loside_transfer <= high_transfer ) ||
              ( low_transfer <= hiside_transfer && hiside_transfer <= high_transfer ) ) {
            // Found a transfer satisfying the range.  Stop searching.
            break;
         }

      }

   }

   /*
    * Swapping does not produce new cuts, so it is ok to omit the penalties
    * arising from cutting.
    */
   double current_balance_penalty = static_cast<double>(ideal_transfer);
   double balance_penalty_loside = static_cast<double>(loside_transfer-ideal_transfer);
   double balance_penalty_hiside = static_cast<double>(hiside_transfer-ideal_transfer);

   if (d_print_swap_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
      tbox::plog.precision(8);
      tbox::plog << "    Swap candidates give penalties (unswap,lo,hi): "
                 << current_balance_penalty << " , " << balance_penalty_loside
                 << " , " << balance_penalty_hiside << std::endl;
   }

   bool found_swap = false;
   TransitSet::iterator isrc = src.end();
   TransitSet::iterator idst = dst.end();
   actual_transfer = 0;

   if ( evaluateBreak( hiside_acceptance_flags, 0, hiside_transfer,
                       ideal_transfer, low_transfer, high_transfer ) ) {
      isrc = src_hiside;
      idst = dst_hiside;
      actual_transfer = hiside_transfer;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking hiside." << std::endl;
      }
   }

   if ( evaluateBreak( loside_acceptance_flags, actual_transfer, loside_transfer,
                       ideal_transfer, low_transfer, high_transfer ) ) {
      isrc = src_loside;
      idst = dst_loside;
      actual_transfer = loside_transfer;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking loside." << std::endl;
      }
   }


   if (found_swap) {

      // We can improve balance_penalty by swapping isrc with idst.
      if (d_print_swap_steps) {
         tbox::plog << "    Swapping " << actual_transfer << " units using ";
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


   } else {
      if (d_print_swap_steps) {
         if ( isrc == src.end() ) {
            tbox::plog << "    Cannot find swap pair for " << ideal_transfer
                       << " units." << std::endl;
         }
         else {
            tbox::plog << "    Keeping original (no swap)." << std::endl;
         }
      }
   }

   t_find_swap_pair->stop();
   return found_swap;
}



/*
 *************************************************************************
 * Master method for breaking off a load.
 *
 * Try different heuristics and pick the "best" way to break off a
 * load.  The best is defined as the one with the lowest combined
 * penalty.
 *
 * This method always return a breakage if at all possible, without
 * considering whether the break should be used.  For example,
 * requesting breakage of 1 cell in a 100x100 box might return a
 * breakage of a 100-cells sliver!
 *
 * Return whether a successful break was made.
 *************************************************************************
 */
bool
TreeLoadBalancer::breakOffLoad(
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load,
   const hier::Box& box,
   double ideal_load_to_break,
   double low_load,
   double high_load ) const
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


   /*
    * To avoid repeated computations of bad cuts,
    * we precompute bad_cuts here to provide to
    * methods that actually use the information.
    */
   std::vector<std::vector<bool> > bad_cuts(d_dim.getValue());
   t_find_bad_cuts->start();
   hier::BoxUtilities::findBadCutPoints(bad_cuts,
      box,
      d_block_domain_boxes[box.getBlockId().getBlockValue()],
      d_bad_interval);
   t_find_bad_cuts->stop();

   // Penalty for not transfering ideal load.
   double best_balance_penalty = computeBalancePenalty(box,
         ideal_load_to_break);

   if (d_print_break_steps) {
      tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
      tbox::plog.precision(6);
      tbox::plog << "      pre-break imbalance: " << ideal_load_to_break
                 << " balance penalty: " << best_balance_penalty
                 << std::endl;
   }

   brk_load = 0;
   bool found_any_break = false;

   {
      std::vector<hier::Box> planar_breakoff;
      std::vector<hier::Box> planar_leftover;
      double planar_brk_load;

      bool found_this_break = breakOffLoad_planar(
            planar_breakoff,
            planar_leftover,
            planar_brk_load,
            box,
            ideal_load_to_break,
            low_load,
            high_load,
            bad_cuts );

      if (found_this_break) {

         found_any_break = true;

         double planar_balance_penalty = computeBalancePenalty(planar_breakoff,
               planar_leftover,
               static_cast<double>(planar_brk_load - ideal_load_to_break));

         if (d_print_break_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "      breakOffLoad_planar broke off "
                       << planar_brk_load << " / " << ideal_load_to_break
                       << " from " << box << '|'
                       << box.numberCells() << '|'
                       << box.size() << " into "
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
                       << " balance penalties: "
                       << planar_balance_penalty
                       << std::endl;
         }


         int break_acceptance_flags[3] = {0,0,0};

         const bool accept_break = evaluateBreak(
            break_acceptance_flags, brk_load, planar_brk_load,
            ideal_load_to_break, low_load, high_load );
         if (d_print_break_steps) {
            tbox::plog << "      Break evaluation:"
                       << "  " << break_acceptance_flags[0]
                       << "  " << break_acceptance_flags[1]
                       << "  " << break_acceptance_flags[2]
                       << std::endl;
         }

         if (accept_break) {
            if (d_print_break_steps) {
               tbox::plog << "      Keeping planar cut result."
                          << "  " << planar_breakoff.size() << " boxes broken off."
                          << "  " << planar_leftover.size() << " boxes leftover."
                          << std::endl;
            }
            breakoff.swap(planar_breakoff);
            leftover.swap(planar_leftover);
            brk_load = planar_brk_load;
            best_balance_penalty = planar_balance_penalty;
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

      std::vector<hier::Box> cubic_breakoff;
      std::vector<hier::Box> cubic_leftover;
      double cubic_brk_load;

      bool found_this_break = breakOffLoad_cubic(
            cubic_breakoff,
            cubic_leftover,
            cubic_brk_load,
            box,
            ideal_load_to_break,
            low_load,
            high_load,
            bad_cuts );

      if (found_this_break) {

         found_any_break = true;

         double cubic_balance_penalty = computeBalancePenalty(
               cubic_breakoff,
               cubic_leftover,
               static_cast<double>(cubic_brk_load - ideal_load_to_break));

         if (d_print_break_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "      breakOffLoad_cubic broke off "
                       << cubic_brk_load << " / " << ideal_load_to_break
                       << " from " << box << '|'
                       << box.numberCells() << '|'
                       << box.size() << " into "
                       << cubic_breakoff.size()
                       << " breakoff: ";
            for (std::vector<hier::Box>::const_iterator bi =
                    cubic_breakoff.begin();
                 bi != cubic_breakoff.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        and " << cubic_leftover.size()
                       << " leftover boxes:";
            for (std::vector<hier::Box>::const_iterator bi =
                    cubic_leftover.begin();
                 bi != cubic_leftover.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        imbalance: "
                       << (cubic_brk_load - ideal_load_to_break)
                       << " balance penalties: "
                       << cubic_balance_penalty
                       << std::endl;
         }

         int break_acceptance_flags[3] = {0,0,0};

         const bool accept_break = evaluateBreak(
            break_acceptance_flags, brk_load, cubic_brk_load,
            ideal_load_to_break, low_load, high_load );
         if (d_print_break_steps) {
            tbox::plog << "      Break evaluation:"
                       << "  " << break_acceptance_flags[0]
                       << "  " << break_acceptance_flags[1]
                       << "  " << break_acceptance_flags[2]
                       << std::endl;
         }

         if (accept_break) {
            if (d_print_break_steps) {
               tbox::plog << "      choosing breakOffLoad_cubic result."
                          << "  " << cubic_breakoff.size() << " boxes broken off."
                          << "  " << cubic_leftover.size() << " boxes leftover."
                          << std::endl;
            }
            breakoff.swap(cubic_breakoff);
            leftover.swap(cubic_leftover);
            brk_load = cubic_brk_load;
            best_balance_penalty = cubic_balance_penalty;
         } else {
            if (d_print_break_steps) {
               tbox::plog << "      Rejecting cubic cut result." << std::endl;
            }
         }
      } else {
         if (d_print_break_steps) {
            tbox::plog << "      breakOffLoad_cubic could not break "
                       << ideal_load_to_break << " from " << box
                       << '/' << box.numberCells()
                       << '/' << box.numberCells().getProduct()
                       << std::endl;
         }
      }

   }

   t_break_off_load->stop();

   return found_any_break;
}



/*
 *************************************************************************
 * Determine whether a proposed break should be accepted based on
 * closeness to ideal and being within a given range.
 *
 * Return values in flags:
 * - [0]: -1, 0 or 1: degrades, leave-alone or improves in-range
 * - [1]: -1, 0 or 1: degrades, leave-alone or improves balance
 * - [2]: 0 or 1: whether new is an overall improvement over current
 *
 * Return whether new_load is an improvement over current_load.
 *************************************************************************
 */

bool
TreeLoadBalancer::evaluateBreak(
   int flags[],
   LoadType cur_load,
   LoadType new_load,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load ) const
{
   LoadType cur_range_miss = cur_load >= high_load ? cur_load-high_load :
      cur_load <= low_load ? low_load-cur_load : 0.0;
   LoadType new_range_miss = new_load >= high_load ? new_load-high_load :
      new_load <= low_load ? low_load-new_load : 0.0;
   flags[0] = new_range_miss < cur_range_miss ? 1 : new_range_miss > cur_range_miss ? -1 : 0;

   LoadType cur_diff = tbox::MathUtilities<double>::Abs(cur_load-ideal_load);
   LoadType new_diff = tbox::MathUtilities<double>::Abs(new_load-ideal_load);

   flags[1] = new_diff < cur_diff ? 1 : new_diff > cur_diff ? -1 : 0;

   /*
    * Combined evaluation gives preference to in-range improvement.
    * If in-range is the same, use balance improvement.
    */
   // flags[2] = flags[0] == 1 ? 1 : flags[1] == 1 ? 1 : 0;
   flags[2] = flags[0] != 0 ? flags[0] : flags[1] != 0 ? flags[1] : 0;

   return flags[2] == 1;
}



/*
 *************************************************************************
 * Measuring surface area of boxes is used in penalizing
 * the creation of new surfaces.
 *************************************************************************
 */

double
TreeLoadBalancer::computeBoxSurfaceArea(
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

int
TreeLoadBalancer::computeBoxSurfaceArea(
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
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double
TreeLoadBalancer::computeSurfacePenalty(
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

double
TreeLoadBalancer::computeSurfacePenalty(
   const TransitSet& a,
   const TransitSet& b) const
{
   double surface_penalty = 0;
   for (TransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->d_box);
   }
   for (TransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->d_box);
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

double
TreeLoadBalancer::computeSurfacePenalty(
   const hier::Box& a) const
{
   int boxvol = a.size();
   double surface_area = computeBoxSurfaceArea(a);
   double best_surface = 2 * d_dim.getValue() * pow(static_cast<double>(boxvol),
         static_cast<double>(d_dim.getValue() - 1) / d_dim.getValue());
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

double
TreeLoadBalancer::computeSlenderPenalty(
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

double
TreeLoadBalancer::computeSlenderPenalty(
   const TransitSet& a,
   const TransitSet& b) const
{
   double slender_penalty = 0;
   for (TransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->d_box);
   }
   for (TransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->d_box);
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

double
TreeLoadBalancer::computeSlenderPenalty(
   const hier::Box& a) const
{
   const hier::IntVector boxdim = a.numberCells();
   double slender_penalty = static_cast<double>(boxdim.max()) / boxdim.min() - 1.0;
   slender_penalty = slender_penalty * slender_penalty; // Make it blow up.
   return slender_penalty;
}



/*
 *************************************************************************
 * Break up box bursty against box solid and adds the pieces to list.
 * This version differs from that in BoxContainer in that it tries to minimize
 * slivers.
 *************************************************************************
 */

void
TreeLoadBalancer::burstBox(
   std::vector<hier::Box>& boxes,
   const hier::Box& bursty,
   const hier::Box& solid ) const
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
   hier::BoxContainer l1(bursty);
   hier::BoxContainer l2(solid);
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        bi++) {
      l2.pushFront(*bi);
   }
   l1.removeIntersections(l2);
   TBOX_ASSERT(l1.isEmpty());
   l2.removeIntersections(bursty);
   TBOX_ASSERT(l2.isEmpty());
#endif
}



/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::computeLocalLoads(
   const hier::BoxLevel& box_level) const
{
   // Count up workload.
   double load = 0.0;
   const hier::BoxContainer& boxes = box_level.getBoxes();
   for (hier::BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end();
        ++ni) {
      double box_load = computeLoad(*ni);
      load += box_load;
   }
   return static_cast<double>(load);
}



/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void
TreeLoadBalancer::printClassData(
   std::ostream& os) const
{
   os << "\nTreeLoadBalancer::printClassData..." << std::endl;
   os << "\nTreeLoadBalancer: this = "
      << (TreeLoadBalancer *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;

   int ln;

   os << "d_workload_data_id..." << std::endl;
   for (ln = 0; ln < static_cast<int>(d_workload_data_id.size()); ln++) {
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

void
TreeLoadBalancer::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{

   if (input_db) {

      d_print_steps = input_db->getBoolWithDefault("DEV_print_steps", false);
      d_print_break_steps =
         input_db->getBoolWithDefault("DEV_print_break_steps", false);
      d_print_swap_steps =
         input_db->getBoolWithDefault("DEV_print_swap_steps", false);
      d_print_edge_steps =
         input_db->getBoolWithDefault("DEV_print_edge_steps", false);
      d_check_connectivity =
         input_db->getBoolWithDefault("DEV_check_connectivity",
            d_check_connectivity);
      d_check_map =
         input_db->getBoolWithDefault("DEV_check_map",
            d_check_map);

      d_summarize_map = input_db->getBoolWithDefault("DEV_summarize_map",
         d_summarize_map);

      d_report_load_balance = input_db->getBoolWithDefault(
         "DEV_report_load_balance",
         d_report_load_balance);
      d_barrier_before = input_db->getBoolWithDefault("DEV_barrier_before",
         d_barrier_before);
      d_barrier_after = input_db->getBoolWithDefault("DEV_barrier_after",
         d_barrier_after);

      d_max_cycle_spread_ratio =
         input_db->getIntegerWithDefault("max_cycle_spread_ratio",
            d_max_cycle_spread_ratio);

      d_flexible_load_tol =
         input_db->getDoubleWithDefault("flexible_load_tolerance",
            d_flexible_load_tol);

      d_min_load_fraction_per_box =
         input_db->getDoubleWithDefault("min_load_fraction_per_box",
            d_min_load_fraction_per_box);
      if ( d_min_load_fraction_per_box < 0 ||
           d_min_load_fraction_per_box >= 1.0 ) {
         TBOX_ERROR("TreeLoadBalancer::getFromInput: min_load_fraction_per_box value of "
                    << d_min_load_fraction_per_box
                    << " is out of range.  It should be >= 0 and < 1 and on the order of 0.01.");
      }

      d_balance_penalty_wt =
         input_db->getDoubleWithDefault("DEV_balance_penalty_wt", 1.0);
      d_surface_penalty_wt =
         input_db->getDoubleWithDefault("DEV_surface_penalty_wt", 1.0);
      d_slender_penalty_wt =
         input_db->getDoubleWithDefault("DEV_slender_penalty_wt", 1.0);
      d_precut_penalty_wt =
         input_db->getDoubleWithDefault("DEV_precut_penalty_wt", 1.0);

   }
}



/*
 ***************************************************************************
 *
 ***************************************************************************
 */
void
TreeLoadBalancer::assertNoMessageForPrivateCommunicator() const
{
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      int flag;
      tbox::SAMRAI_MPI::Status mpi_status;
      int mpi_err = d_mpi.Iprobe(MPI_ANY_SOURCE,
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
 * the box's longest direction.
 *
 * Currently assuming uniform load of one unit per cell.
 *
 * Return whether a successful break was made.
 *************************************************************************
 */
bool
TreeLoadBalancer::breakOffLoad_planar(
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load,
   const hier::Box& box,
   double ideal_brk_load,
   double low_load,
   double high_load,
   const std::vector<std::vector<bool> >& bad_cuts ) const
{

   const tbox::Dimension dim(d_dim);

   if (d_print_break_steps) {
      tbox::plog << "      breakOffLoad_planar attempting to break "
                 << ideal_brk_load << " from Box "
                 << box << box.numberCells() << '|' << box.size()
                 << " min_size=" << d_min_size << std::endl;
   }

   breakoff.clear();
   leftover.clear();

   const hier::IntVector& box_dims = box.numberCells();

   const int box_vol = box_dims.getProduct();

   if (box_vol <= ideal_brk_load) {
      // Easy: break off everything.
      breakoff.push_back(box);
      brk_load = box_vol;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffLoad_planar broke off entire Box "
                    << box
                    << std::endl;
      }
      return true;
   }

   /*
    * Determine ordering of box_dims from shortest to longest.
    */
   hier::IntVector sorted_dirs(dim);
   sorted_dirs.sortIntVector(box_dims);

   hier::Box best_breakoff_box(dim);
   hier::Box best_leftover_box(dim);
   best_breakoff_box.setBlockId(box.getBlockId());
   best_leftover_box.setBlockId(box.getBlockId());
   brk_load = 0;

   int break_acceptance_flags[3] = {0,0,0};
   bool sufficient_brk_load = false;

   for (int d = d_dim.getValue() - 1; d >= 0 && !sufficient_brk_load; --d) {

      /*
       * Search directions from longest to shortest because we prefer
       * to break across longer directions.
       */
      const int brk_dir = sorted_dirs(d);

      const int brk_area = box_vol / box_dims(brk_dir);

      const std::vector<bool>& bad = bad_cuts[brk_dir];

      const double ideal_cut_length = double(ideal_brk_load)/brk_area;

      /*
       * Try 4 different cuts for direction brk_dir:
       * upper/lower: should we break off the upper end or lower end of the box.
       * hi/lo: should we round the break plane to the high or low side.
       *
       * plane refers to the index of the mesh line where the cut is.
       */

      // Ideal cut planes, not necessarily coincident with a grid line.
      const double ideal_upper_cut_plane = box.upper()(brk_dir) + 1 - ideal_cut_length;
      const double ideal_lower_cut_plane = box.lower()(brk_dir) + ideal_cut_length;

      // Compute valid cut planes on high and low sides of upper cut plane.
      int lo_upper_cut_plane = int(ideal_upper_cut_plane);
      int hi_upper_cut_plane = int(ideal_upper_cut_plane) + 1;
      lo_upper_cut_plane = ROUND_TO_LO(lo_upper_cut_plane, d_cut_factor(brk_dir));
      hi_upper_cut_plane = ROUND_TO_HI(hi_upper_cut_plane, d_cut_factor(brk_dir));
      while ( lo_upper_cut_plane > box.lower()(brk_dir)   && bad[lo_upper_cut_plane-box.lower()(brk_dir)] ) { lo_upper_cut_plane -= d_cut_factor(brk_dir); }
      while ( hi_upper_cut_plane < box.upper()(brk_dir)+1 && bad[hi_upper_cut_plane-box.lower()(brk_dir)] ) { hi_upper_cut_plane += d_cut_factor(brk_dir); }

      // Compute valid cut planes on high and low sides of lower cut plane.
      int lo_lower_cut_plane = int(ideal_lower_cut_plane);
      int hi_lower_cut_plane = int(ideal_lower_cut_plane) + 1;
      lo_lower_cut_plane = ROUND_TO_LO(lo_lower_cut_plane, d_cut_factor(brk_dir));
      hi_lower_cut_plane = ROUND_TO_HI(hi_lower_cut_plane, d_cut_factor(brk_dir));
      while ( lo_lower_cut_plane > box.lower()(brk_dir)   && bad[lo_lower_cut_plane-box.lower()(brk_dir)] ) { lo_lower_cut_plane -= d_cut_factor(brk_dir); }
      while ( hi_lower_cut_plane < box.upper()(brk_dir)+1 && bad[hi_lower_cut_plane-box.lower()(brk_dir)] ) { hi_lower_cut_plane += d_cut_factor(brk_dir); }


      if ( lo_lower_cut_plane - box.lower()(brk_dir) > d_min_size(brk_dir) &&
           box.upper()(brk_dir)+1 - lo_lower_cut_plane > d_min_size(brk_dir) ) {

         const int lo_lower_cut_vol = brk_area*( lo_lower_cut_plane - box.lower()(brk_dir) );

         if ( evaluateBreak( break_acceptance_flags, brk_load, lo_lower_cut_vol,
                             ideal_brk_load, low_load, high_load ) ) {
            brk_load = lo_lower_cut_vol;
            best_breakoff_box = best_leftover_box = box;
            best_breakoff_box.upper()(brk_dir) = lo_lower_cut_plane - 1;
            best_leftover_box.lower()(brk_dir) = lo_lower_cut_plane;
            TBOX_ASSERT( best_breakoff_box.size() == lo_lower_cut_vol );
         }
      }

      if ( ( hi_lower_cut_plane - box.lower()(brk_dir) > d_min_size(brk_dir) &&
           box.upper()(brk_dir)+1 - hi_lower_cut_plane > d_min_size(brk_dir) ) ||
           hi_lower_cut_plane >= box.upper()(brk_dir)+1 ) {

         const int hi_lower_cut_vol = brk_area*( hi_lower_cut_plane - box.lower()(brk_dir) );

         if ( evaluateBreak( break_acceptance_flags, brk_load, hi_lower_cut_vol,
                             ideal_brk_load, low_load, high_load ) ) {
            brk_load = hi_lower_cut_vol;
            best_breakoff_box = best_leftover_box = box;
            best_breakoff_box.upper()(brk_dir) = hi_lower_cut_plane - 1;
            best_leftover_box.lower()(brk_dir) = hi_lower_cut_plane;
            TBOX_ASSERT( best_breakoff_box.size() == hi_lower_cut_vol );
         }
      }

      if ( ( box.upper()(brk_dir)+1 - lo_upper_cut_plane > d_min_size(brk_dir) &&
           lo_upper_cut_plane - box.lower()(brk_dir) > d_min_size(brk_dir) ) ||
           lo_upper_cut_plane <= box.lower()(brk_dir) ) {

         const int lo_upper_cut_vol = brk_area*( box.upper()(brk_dir)+1 - lo_upper_cut_plane );

         if ( evaluateBreak( break_acceptance_flags, brk_load, lo_upper_cut_vol,
                          ideal_brk_load, low_load, high_load ) ) {
            brk_load = lo_upper_cut_vol;
            best_breakoff_box = best_leftover_box = box;
            best_breakoff_box.lower()(brk_dir) = lo_upper_cut_plane;
            best_leftover_box.upper()(brk_dir) = lo_upper_cut_plane - 1;
            TBOX_ASSERT( best_breakoff_box.size() == lo_upper_cut_vol );
         }
      }

      if ( box.upper()(brk_dir)+1 - hi_upper_cut_plane > d_min_size(brk_dir) &&
           hi_upper_cut_plane - box.lower()(brk_dir) > d_min_size(brk_dir) ) {

         const int hi_upper_cut_vol = brk_area*( box.upper()(brk_dir)+1 - hi_upper_cut_plane );

         if ( evaluateBreak( break_acceptance_flags, brk_load, hi_upper_cut_vol,
                             ideal_brk_load, low_load, high_load ) ) {
            brk_load = hi_upper_cut_vol;
            best_breakoff_box = best_leftover_box = box;
            best_breakoff_box.lower()(brk_dir) = hi_upper_cut_plane;
            best_leftover_box.upper()(brk_dir) = hi_upper_cut_plane - 1;
            TBOX_ASSERT( best_breakoff_box.size() == hi_upper_cut_vol );
         }
      }

      sufficient_brk_load = (brk_load >= low_load) && (brk_load <= high_load);

   }

   bool successful_break = false;
   if (!best_breakoff_box.empty()) {
      breakoff.push_back(best_breakoff_box);
      TBOX_ASSERT( brk_load == best_breakoff_box.size() );
      successful_break = true;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffLoad_planar broke box " << box << box.numberCells()
                    << " for breakoff box " << best_breakoff_box << best_breakoff_box.numberCells()
                    << " and leftover " << best_leftover_box << best_leftover_box.numberCells()
                    << std::endl;
      }
   } else {
      if (d_print_break_steps) {
         tbox::plog << "      breakOffLoad_planar could not break "
                    << ideal_brk_load << " from Box " << box
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
               << "ideal brk load " << ideal_brk_load);
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
               << "ideal brk load " << ideal_brk_load);
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
bool
TreeLoadBalancer::breakOffLoad_cubic(
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load,
   const hier::Box& box,
   double ideal_brk_load,
   double low_load,
   double high_load,
   const std::vector<std::vector<bool> >& bad_cuts ) const
{

   const hier::IntVector box_dims(box.numberCells());

   const double box_load(box_dims.getProduct());

   if (ideal_brk_load >= box_load) {
      // Easy: break off everything.
      leftover.clear();
      breakoff.clear();
      breakoff.push_back(box);
      brk_load = box_load;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffLoad_cubic broke off entire Box "
                    << box
                    << std::endl;
      }
      return true;
   }

   if (ideal_brk_load > 0.5 * box_load) {
      /*
       * This algorithm is better when breaking off a small portion.
       * Since the ideal is a bigger portion, switch breakoff with leftover.
       */
      if (d_print_break_steps) {
         tbox::plog
         << "      breakOffLoad_cubic reversing direction to break "
         << (box_dims.getProduct() - ideal_brk_load)
         << " instead of " << ideal_brk_load << " / "
         << box_dims.getProduct() << std::endl;
      }
      bool success =
         breakOffLoad_cubic(
            leftover,
            breakoff,
            brk_load,
            box,
            box_dims.getProduct() - ideal_brk_load,
            box_dims.getProduct() - high_load,
            box_dims.getProduct() - low_load,
            bad_cuts );
      if (success) {
         brk_load = box_dims.getProduct() - brk_load;
      }
      return success;
   }

   if (d_print_break_steps) {
      tbox::plog << "      breakOffLoad_cubic attempting to break "
                 << ideal_brk_load << " from Box "
                 << box << box.numberCells() << '|' << box.size()
                 << " min_size=" << d_min_size << std::endl;
   }

   breakoff.clear();
   leftover.clear();
   brk_load = 0.0;

   const hier::IntVector &one_vec = hier::IntVector::getOne(d_dim);
   const hier::IntVector &zero_vec = hier::IntVector::getZero(d_dim);

   hier::Box best_breakoff_box(d_dim);
   hier::IntVector best_breakoff_size = zero_vec;
   LoadType best_breakoff_load = 0;

   /*
    * We consider 2^dim boxes grown from the incoming box's corners.
    * Imagine 2 cutting planes per dimension, perpendicular to each
    * axis, an upper cut and a lower cut.  The box is cut into 3 parts
    * in each direction, so the 2*dim planes divide the box into 3^dim
    * boxes.
    *
    * 2D example:
    *
    *       +-----------------+
    *       |    |        |   |
    *       |----+--------+---|
    *       |    |        |   |
    *       |    |        |   |
    *       |    |        |   |
    *       |    |        |   |
    *       |----+--------+---|
    *       |    |        |   |
    *       +-----------------+
    *
    * upper_intersection is the point where all upper cuts intersect.
    * lower_intersection is the point where all lower cuts intersect.
    * We only consider the but boxes that touch the incoming box's
    * corners.  Using the other boxes result in too much fragmentation
    * of the incoming box.
    */
   hier::IntVector brk_size(d_min_size);
   brk_size.max(d_cut_factor);
   brk_size.min(box_dims);

   /*
    * Make sure brk_size is a multiple of d_cut_factor.
    */
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (brk_size(d) % d_cut_factor(d) != 0) {
         brk_size(d) = ((brk_size(d) / d_cut_factor(d)) + 1) * d_cut_factor(d);
      }
   }

   /*
    * Initialize the cut plane intersections.  We will grow the
    * corner boxes by gradually moving the intersections away
    * from their initial location.
    */
   hier::IntVector lower_intersection(box.lower() + d_min_size);
   hier::IntVector upper_intersection(box.upper() - d_min_size + one_vec);
   for ( int d=0; d<d_dim.getValue(); ++d ) {
      lower_intersection(d) = ROUND_TO_HI( lower_intersection(d),
                                           d_cut_factor(d) );
      upper_intersection(d) = ROUND_TO_LO( upper_intersection(d),
                                           d_cut_factor(d) );
   }


   const int num_corners = 1 << d_dim.getValue();

   for ( int bn=0; bn<num_corners; ++bn ) {

      if ( d_print_break_steps ) {
         tbox::plog << "Examining corner box " << bn << std::endl;
      }

      /*
       * Compute initial box at corner bn and its expansion rate.
       */
      hier::Box corner_box(box);
      hier::IntVector corner_box_size = zero_vec;
      LoadType corner_box_load = 0;
      hier::IntVector expansion_rate(d_dim);

      for ( int d=0; d<d_dim.getValue(); ++d ) {

         // In direction d, does corner_box touch the upper (vs lower) side of box:
         int touches_upper_side = bn & (1 << d) ;

         if ( touches_upper_side ) {
            corner_box.lower()(d) = upper_intersection(d);
            if ( corner_box.lower()(d) - box.lower()(d) < d_min_size(d) ) {
               corner_box.lower()(d) = box.lower()(d);
            }
            expansion_rate(d) = -d_cut_factor(d);
         }
         else {
            corner_box.upper()(d) = lower_intersection(d) - 1;
            if ( box.upper()(d) - corner_box.upper()(d) < d_min_size(d) ) {
               corner_box.upper()(d) = box.upper()(d);
            }
            expansion_rate(d) = d_cut_factor(d);
         }

      }

      corner_box_size = corner_box.numberCells();
      corner_box_load = corner_box.size();

      if ( d_print_break_steps ) {
         tbox::plog << "Initial corner box " << bn << " is " << corner_box << corner_box.numberCells() << '|' << corner_box.size() << std::endl;
      }


      int break_acceptance_flags[3] = {0,0,0};

      if ( evaluateBreak( break_acceptance_flags, best_breakoff_load, corner_box_load,
                          ideal_brk_load, low_load, high_load) ) {
         best_breakoff_box = corner_box;
         best_breakoff_size = corner_box_size;
         best_breakoff_load = corner_box_load;
         if ( best_breakoff_load >= low_load && brk_load <= high_load ) {
            break;
         }
      }



      /*
       * stop_growing: whether corner_box_size is already
       * big engough so that it cannot not be grown without breaking
       * off too much.
       */
      hier::IntVector growable(d_dim, 1);
      for (int d = 0; d < d_dim.getValue(); ++d) {
         growable[d] = corner_box_size[d] < box_dims[d];
      }

      while ( ( best_breakoff_load < low_load || best_breakoff_load > high_load ) &&
              corner_box_load < ideal_brk_load ) {
         /*
          * The while loop gradually increases corner_box to bring
          * its size closer to ideal_brk_load.  Stop loop when
          * its size is in the acceptable range or if increasing
          * it only takes it farther from ideal_brk_load.
          *
          * Select inc_dir, the direction to expand corner_box.  Use the
          * smallest direction that is still allowed to grow.
          */

         int inc_dir = -1;
         for (int d = 0; d < d_dim.getValue(); ++d) {
            if ( growable(d) &&
                 (inc_dir == -1 || corner_box_size(d) < corner_box_size(inc_dir)) )
               inc_dir = d;
         }
         if (inc_dir == -1) break; // No growable direction.

         TBOX_ASSERT(corner_box_size(inc_dir) < box_dims(inc_dir));


         /*
          * Grow corner_box, but keep within boundary of box and
          * prevent remainder from violating min size.  Update growability.
          */
         if ( expansion_rate(inc_dir) > 0 ) {
            corner_box.upper()(inc_dir) = tbox::MathUtilities<int>::Min(
               corner_box.upper()(inc_dir) + expansion_rate(inc_dir),
               box.upper()(inc_dir) );
            if ( box.upper()(inc_dir) - corner_box.upper()(inc_dir) < d_min_size(inc_dir) ) {
               corner_box.upper()(inc_dir) = box.upper()(inc_dir);
            }
            growable(inc_dir) = corner_box.upper()(inc_dir) < box.upper()(inc_dir);
         }
         else {
            corner_box.lower()(inc_dir) = tbox::MathUtilities<int>::Max(
               corner_box.lower()(inc_dir) + expansion_rate(inc_dir),
               box.lower()(inc_dir) );
            if ( corner_box.lower()(inc_dir) - box.lower()(inc_dir) < d_min_size(inc_dir) ) {
               corner_box.lower()(inc_dir) = box.lower()(inc_dir);
            }
            growable(inc_dir) = corner_box.lower()(inc_dir) > box.lower()(inc_dir);
         }
         corner_box_size = corner_box.numberCells();
         corner_box_load = corner_box.size();

         if ( d_print_break_steps ) {
            tbox::plog << "Expand corner box " << bn << " to " << corner_box << corner_box.numberCells() << '|' << corner_box.size() << std::endl;
         }


         const bool accept_break = evaluateBreak(
            break_acceptance_flags, best_breakoff_load, corner_box_load,
            ideal_brk_load, low_load, high_load );

         if ( accept_break ) {
            best_breakoff_box = corner_box;
            best_breakoff_size = corner_box_size;
            best_breakoff_load = corner_box_load;
            if ( best_breakoff_load >= low_load && best_breakoff_load <= high_load ) {
               break;
            }
         }

      } // while loop

   } // bn loop


   if ( !best_breakoff_box.empty() ) {
      brk_load = best_breakoff_load;
      breakoff.push_back(best_breakoff_box);

      burstBox(
         leftover,
         box,
         best_breakoff_box );
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
               << "breakoff box " << b << ", with size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "orig box " << box << "\n"
               << "break box " << b << "\n"
               << "break box size " << b.size() << "\n"
               << "ideal brk load " << ideal_brk_load);
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
               << "orig box " << box << "\n"
               << "break box " << b << "\n"
               << "break box size " << b.size() << "\n"
               << "ideal brk load " << ideal_brk_load);
         }
      }
   }
#endif

   return !breakoff.empty();
}






/*
**************************************************************************
* Move Boxes in balance_box_level from ranks outside of
* rank_group to ranks inside rank_group.  Modify the given connectors
* to make them correct following this moving of boxes.
**************************************************************************
*/

void
TreeLoadBalancer::prebalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_anchor,
   const tbox::RankGroup& rank_group) const
{

   if (balance_to_anchor) {
      TBOX_ASSERT(balance_to_anchor->hasTranspose());
      TBOX_ASSERT(balance_to_anchor->getTranspose().checkTransposeCorrectness(*balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor->checkTransposeCorrectness(balance_to_anchor->getTranspose()) == 0);
   }

   /*
    * tmp_box_level will contain the same boxes as
    * balance_box_level, but all will live on the processors
    * specified in rank_group.
    */
   hier::BoxLevel tmp_box_level(balance_box_level.getRefinementRatio(),
      balance_box_level.getGridGeometry(),
      balance_box_level.getMPI());

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
   tbox::AsyncCommPeer<int>* box_send = 0;
   tbox::AsyncCommPeer<int>* box_recv = 0;
   tbox::AsyncCommPeer<int>* id_send = 0;
   tbox::AsyncCommPeer<int>* id_recv = 0;

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
      std::list<int> recv_ranks;
      for (int i = 0; i < d_mpi.getSize(); i++) {
         if (!rank_group.isMember(i) &&
             rank_group.getMappedRank(i % output_nproc) == d_mpi.getRank()) {
            recv_ranks.push_back(i);
         }
      }
      num_recvs = static_cast<int>(recv_ranks.size());
      if (num_recvs > 0) {
         box_recv = new tbox::AsyncCommPeer<int>[num_recvs];
         id_send = new tbox::AsyncCommPeer<int>[num_recvs];
         int recv_count = 0;
         for (std::list<int>::const_iterator ri(recv_ranks.begin());
              ri != recv_ranks.end(); ri++) {
            const int rank = *ri;
            box_recv[recv_count].initialize(&comm_stage);
            box_recv[recv_count].setPeerRank(rank);
            box_recv[recv_count].setMPI(d_mpi);
            box_recv[recv_count].setMPITag(TreeLoadBalancer_PREBALANCE0 + 2 * rank,
               TreeLoadBalancer_PREBALANCE1 + 2 * rank);

            id_send[recv_count].initialize(&comm_stage);
            id_send[recv_count].setPeerRank(rank);
            id_send[recv_count].setMPI(d_mpi);
            id_send[recv_count].setMPITag(TreeLoadBalancer_PREBALANCE0 + 2 * rank,
               TreeLoadBalancer_PREBALANCE1 + 2 * rank);

            recv_count++;
         }
         TBOX_ASSERT(num_recvs == recv_count);
      }
   }

   /*
    * Construct the mapping Connectors which describe the mapping from the box
    * configuration of the given balance_box_level, to the new
    * configuration stored in tmp_box_level.  These mapping Connectors
    * are necessary to modify the two Connectors given in the argument list,
    * so that on return from this method, they will be correct for the new
    * balance_box_level.
    */
   hier::MappingConnector balance_to_tmp(balance_box_level,
         tmp_box_level,
         hier::IntVector::getZero(d_dim));

   hier::MappingConnector tmp_to_balance(tmp_box_level,
         balance_box_level,
         hier::IntVector::getZero(d_dim));

   balance_to_tmp.setTranspose(&tmp_to_balance, false);

   /*
    * Where Boxes already exist on ranks in rank_group,
    * move them directly to tmp_box_level.
    */
   if (!is_sending_rank) {
      const hier::BoxContainer& unchanged_boxes =
         balance_box_level.getBoxes();

      for (hier::BoxContainer::const_iterator ni = unchanged_boxes.begin();
           ni != unchanged_boxes.end(); ++ni) {

         const hier::Box& box = *ni;
         tmp_box_level.addBox(box);
      }
   }

   const int buf_size = hier::Box::commBufferSize(d_dim);

   /*
    * On sending ranks, pack the Boxes into buffers and send.
    */
   if (is_sending_rank) {
      const hier::BoxContainer& sending_boxes =
         balance_box_level.getBoxes();
      const int num_sending_boxes =
         static_cast<int>(sending_boxes.size());

      int* buffer = new int[buf_size * num_sending_boxes];
      int box_count = 0;
      for (hier::BoxContainer::const_iterator ni = sending_boxes.begin();
           ni != sending_boxes.end(); ++ni) {

         const hier::Box& box = *ni;

         box.putToIntBuffer(&buffer[box_count * buf_size]);
         box_count++;
      }
      box_send->beginSend(buffer, buf_size * num_sending_boxes);

      delete[] buffer;
   }

   /*
    * On receiving ranks, complete the receives, add the boxes to local
    * tmp_box_level, insert boxes into tmp_to_balance, and then
    * send the new LocalIdes back to the sending processors.
    */
   if (!is_sending_rank && num_recvs > 0) {
      for (int i = 0; i < num_recvs; i++) {
         box_recv[i].beginRecv();
      }
      int num_completed_recvs = 0;
      std::vector<bool> completed(num_recvs, false);
      while (num_completed_recvs < num_recvs) {
         for (int i = 0; i < num_recvs; i++) {
            if (!completed[i] && box_recv[i].checkRecv()) {
               num_completed_recvs++;
               completed[i] = true;
               const int num_boxes = box_recv[i].getRecvSize() / buf_size;
               const int* buffer = box_recv[i].getRecvData();
               int* id_buffer = new int[num_boxes];

               for (int b = 0; b < num_boxes; b++) {
                  hier::Box box(d_dim);

                  box.getFromIntBuffer(&buffer[b * buf_size]);

                  hier::BoxContainer::const_iterator tmp_iter =
                     tmp_box_level.addBox(box,
                        box.getBlockId());

                  hier::BoxId tmp_box_id = tmp_iter->getBoxId();

                  tmp_to_balance.insertLocalNeighbor(box, tmp_box_id);

                  id_buffer[b] = tmp_box_id.getLocalId().getValue();
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
    * to balance_to_tmp.
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

      const hier::BoxContainer& sending_boxes =
         balance_box_level.getBoxes();
      TBOX_ASSERT(static_cast<int>(id_recv->getRecvSize()) == sending_boxes.size());

      int box_count = 0;
      for (hier::BoxContainer::const_iterator ni = sending_boxes.begin();
           ni != sending_boxes.end(); ++ni) {

         hier::Box new_box(
            *ni,
            (hier::LocalId)buffer[box_count],
            rank_group.getMappedRank(d_mpi.getRank() % output_nproc));

         balance_to_tmp.insertLocalNeighbor(new_box, (*ni).getBoxId());
         box_count++;
      }
   }

   if (balance_to_anchor && balance_to_anchor->hasTranspose()) {
      /*
       * This modify operation copies tmp_box_level to
       * balance_box_level, and changes balance_to_anchor and
       * its transpose such that they are correct for the new state
       * of balance_box_level.
       */
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix(d_object_name);
      mca.modify(balance_to_anchor->getTranspose(),
         balance_to_tmp,
         &balance_box_level,
         &tmp_box_level);

      TBOX_ASSERT(balance_to_anchor->getTranspose().checkTransposeCorrectness(*balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor->checkTransposeCorrectness(balance_to_anchor->getTranspose()) == 0);
   } else {
      hier::BoxLevel::swap(balance_box_level, tmp_box_level);
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

std::ostream&
operator << (
   std::ostream& co,
   const TreeLoadBalancer::BoxInTransit& r)
{
   co << r.d_box
   << r.d_box.numberCells() << '|'
   << r.d_box.numberCells().getProduct();
   co << '-' << r.d_orig_box
   << r.d_box.numberCells() << '|'
   << r.d_box.numberCells().getProduct();
   return co;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
TreeLoadBalancer::setTimers()
{
   /*
    * The first constructor gets timers from the TimerManager.
    * and sets up their deallocation.
    */
   if (!t_load_balance_box_level) {
      t_load_balance_box_level = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::loadBalanceBoxLevel()");
      t_get_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_map");
      t_use_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::use_map");
      t_constrain_size = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constrain_size");
      t_map_big_boxes = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::mapOversizedBoxes()");
      t_load_distribution = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::load_distribution");
t_post_load_distribution_barrier = tbox::TimerManager::getManager()->
   getTimer(d_object_name + "::post_load_distribution_barrier");
      t_compute_local_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_local_load");
      t_compute_global_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_global_load");
      t_compute_tree_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load");

      const int max_cycles_to_time = 4;
      t_compute_tree_load_for_cycle.resize(
         max_cycles_to_time,
         boost::shared_ptr<tbox::Timer>() );
      for ( int i=0; i<max_cycles_to_time; ++i ) {
         t_compute_tree_load_for_cycle[i] = tbox::TimerManager::getManager()->
            getTimer(d_object_name + "::compute_tree_load_for_cycle["
                     + tbox::Utilities::intToString(i) + "]");
      }

      t_break_off_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::breakOffLoad()");
      t_find_bad_cuts = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::find_bad_cuts");
      t_adjust_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::adjustLoad()");
      t_adjust_load_by_swapping = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::adjustLoadBySwapping()");
      t_shift_loads_by_breaking = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::adjustLoadByBreaking()");
      t_find_swap_pair = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::swapLoadPair()");
      t_send_load_to_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_children");
      t_send_load_to_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_parent");
      t_get_load_from_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_children");
      t_get_load_from_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_parent");
      t_construct_semilocal = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constructSemilocalUnbalancedToBalanced()");
      t_construct_semilocal_comm_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constructSemilocalUnbalancedToBalanced()_comm_wait");
      t_report_loads = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::report_loads");
      t_finish_sends = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::finish_sends");
      t_local_balancing = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::local_balancing");
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
      t_child_send_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::child_send_wait");
      t_child_recv_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::child_recv_wait");
      t_parent_send_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_send_wait");
      t_parent_recv_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_recv_wait");
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const tbox::Dimension &dim) :
   d_box(dim),
   d_orig_box(dim)
{
}



/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const hier::Box& origin):
   d_box(origin),
   d_orig_box(origin),
   d_boxload(origin.size())
{
}



/*
 *************************************************************************
 * Construct a new BoxInTransit with the history of an existing box.
 *************************************************************************
 */
TreeLoadBalancer::BoxInTransit::BoxInTransit(
   const BoxInTransit& other,
   const hier::Box& box,
   int rank,
   hier::LocalId local_id):
   d_box(box, local_id, rank),
   d_orig_box(other.d_orig_box),
   d_boxload(d_box.size())
{
}

void
TreeLoadBalancer::BoxInTransit::putToMessageStream(
   tbox::MessageStream &mstream ) const
{
   d_box.putToMessageStream(mstream);
   d_orig_box.putToMessageStream(mstream);
   mstream << d_boxload;
   return;
}

void
TreeLoadBalancer::BoxInTransit::getFromMessageStream(
   tbox::MessageStream &mstream )
{
   d_box.getFromMessageStream(mstream);
   d_orig_box.getFromMessageStream(mstream);
   mstream >> d_boxload;
   return;
}

TreeLoadBalancer::SubtreeData::SubtreeData():
   d_subtree_rank(-1),
   d_num_procs(0),
   d_subtree_load_current(0),
   d_subtree_load_ideal(-1),
   d_subtree_load_upperlimit(-1),
   d_eff_num_procs(0),
   d_eff_load_current(0),
   d_eff_load_ideal(-1),
   d_eff_load_upperlimit(-1),
   d_work_traded(),
   d_wants_work_from_parent(false)
{
}

void
TreeLoadBalancer::SubtreeData::addChild(
   const SubtreeData &child )
{
   /*
    * Sum children load into my_subtree to get data for the whole
    * subtree.
    */

   d_num_procs += child.d_num_procs;
   d_subtree_load_current += child.d_subtree_load_current;
   d_subtree_load_upperlimit += child.d_subtree_load_upperlimit;
   d_subtree_load_ideal += child.d_subtree_load_ideal;

   if ( child.d_wants_work_from_parent ) {
      d_eff_num_procs += child.d_eff_num_procs;
      d_eff_load_current += child.d_eff_load_current;
      d_eff_load_upperlimit += child.d_eff_load_upperlimit;
      d_eff_load_ideal += child.d_eff_load_ideal;
   }

   d_subtree_load_current += child.d_work_traded.getSumLoad();
   d_eff_load_current += child.d_work_traded.getSumLoad();
}

void
TreeLoadBalancer::SubtreeData::printClassData(
   const std::string &border,
   std::ostream &os ) const
{
   os.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
   os.precision(6);
   os << border
      << "Full nproc = " << d_num_procs
      << "   current = " << d_subtree_load_current
      << "   ideal = " << d_subtree_load_ideal
      << "   ratio = " << (d_subtree_load_current/d_subtree_load_ideal)
      << "   avg = " << (d_subtree_load_current / d_num_procs)
      << "   upperlimit = " << d_subtree_load_upperlimit
      << "   surplus = " << surplus()
      << "   excess =  " << excess()
      << '\n' << border
      << "Effective nproc = " << d_eff_num_procs
      << "   current = " << d_eff_load_current
      << "   ideal = " << d_eff_load_ideal
      << "   ratio = " << (d_eff_load_current/d_eff_load_ideal)
      << "   avg = " << (d_eff_load_current / d_eff_num_procs)
      << "   upperlimit = " << d_eff_load_upperlimit
      << "   surplus = " << effSurplus()
      << "   excess =  " << effExcess()
      << '\n' << border
      << "load traded =  " << d_work_traded.getSumLoad()
      << "   wants work from parent = " << d_wants_work_from_parent
      << '\n';
   return;
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
