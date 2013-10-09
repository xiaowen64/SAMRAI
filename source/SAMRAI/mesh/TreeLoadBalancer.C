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

#include "SAMRAI/mesh/BalanceBoxBreaker.h"

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
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
   d_allow_box_breaking(true),
   d_rank_tree(rank_tree ? rank_tree : boost::shared_ptr<tbox::RankTreeStrategy>(new tbox::CenteredRankTree) ),
   d_comm_graph_writer(),
   d_master_workload_data_id(d_default_data_id),
   d_flexible_load_tol(0.0),
   d_load_comparison_tol(1.0e-5),
   d_balance_penalty_wt(1.0),
   d_surface_penalty_wt(1.0),
   d_slender_penalty_wt(1.0),
   d_slender_penalty_threshold(3.0),
   d_precut_penalty_wt(1.0),
   // Output control.
   d_report_load_balance(false),
   d_summarize_map(false),
   // Performance evaluation.
   d_barrier_before(false),
   d_barrier_after(false),
   d_print_steps(false),
   d_print_pop_steps(false),
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

   d_pparams = boost::make_shared<PartitioningParams>(
      *balance_box_level.getGridGeometry(),
      balance_box_level.getRefinementRatio(),
      min_size, max_size, bad_interval, cut_factor);

   /*
    * We expect the domain box_level to be in globalized state.
    */
   TBOX_ASSERT(
      domain_box_level.getParallelState() ==
      hier::BoxLevel::GLOBALIZED);

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
    * If max_size is given (positive), constrain boxes to the given
    * max_size.  If not given, skip the enforcement step to save some
    * communications.
    */

   hier::IntVector max_intvector(d_dim, tbox::MathUtilities<int>::getMax());
   if (max_size != max_intvector) {

      t_constrain_size->barrierAndStart();
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

   d_pparams.reset();

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
            "Errors in load balance mapping found.\n"
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

      if (box_size <= d_pparams->getMaxBoxSize()) {

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
            d_pparams->getMaxBoxSize(),
            d_pparams->getMinBoxSize(),
            d_pparams->getCutFactor(),
            d_pparams->getBadInterval(),
            d_pparams->getDomainBoxes(box.getBlockId()));
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
                 << "  Avg is " << group_avg_load/d_pparams->getMinBoxSize().getProduct()
                 << " times min size of " << d_pparams->getMinBoxSize()
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
    * values can never be the same.  Moreover, imported boxes
    * always take indices from the set associated with the importer,
    * independent their arrival order.
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

   std::vector<hier::SequentialLocalIdGenerator> id_generator(2 + d_rank_tree->getDegree());
   id_generator[0].setLastValue( next_available_index[0] );
   id_generator[0].setIncrement( hier::LocalId(id_generator.size()) );
   for (unsigned int c = 1; c < d_rank_tree->getDegree() + 2; ++c) {
      id_generator[c].setLastValue( next_available_index[0] + c );
      id_generator[c].setIncrement( hier::LocalId(id_generator.size()) );
   }


   /*
    * Data for storing and transfering subtree info.
    */
   SubtreeData my_subtree;
   my_subtree.setPartitioningParams(*d_pparams);
   std::vector<SubtreeData> child_subtrees(num_children);
   for ( size_t i=0; i<child_subtrees.size(); ++i ) {
      child_subtrees[i].setPartitioningParams(*d_pparams);
   }


   /*
    * unassigned is a container of BoxTransitSet::BoxInTransit that has been released by
    * a process and has not yet been assigned to another.  First, put
    * all initial local work in unassigned.  Imported
    * BoxTransitSet::BoxInTransits are placed here before determining whether to keep
    * them or send them to another part of the tree.
    */
   BoxTransitSet unassigned;
   unassigned.setPartitioningParams(*d_pparams);
   unassigned.insertAll(balance_box_level.getBoxes());


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
    * Step 2, remote part:
    *
    * Finish getting tree and load data from children.
    * Add imported BoxTransitSet::BoxInTransit to unassigned bin.
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
         id_generator[cindex],
         mstream);

      unassigned.insertAll( child_subtrees[cindex].d_work_traded );

      if (d_print_steps) {
         tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
         tbox::plog.precision(6);
         tbox::plog << "Received from child "
                    << cindex << ':' << d_rank_tree->getChildRank(cindex) << ": "
                    << child_subtrees[cindex].d_work_traded.size() << " boxes ("
                    << child_subtrees[cindex].d_work_traded.getSumLoad() << " units):";
         if ( child_subtrees[cindex].d_work_traded.size() < 10 ) {
            for ( BoxTransitSet::const_iterator ni=child_subtrees[cindex].d_work_traded.begin();
                  ni!=child_subtrees[cindex].d_work_traded.end(); ++ni ) {
               const BoxTransitSet::BoxInTransit& box_in_transit = *ni;
               tbox::plog << "  " << box_in_transit;
            }
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

            const LoadType export_load_actual =
               my_subtree.d_work_traded.adjustLoad(
                  unassigned,
                  next_available_index[d_rank_tree->getDegree()],
                  id_generator[d_rank_tree->getDegree()],
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
               if ( my_subtree.d_work_traded.size() < 10 ) {
                  for (BoxTransitSet::const_iterator ni = my_subtree.d_work_traded.begin();
                       ni!=my_subtree.d_work_traded.end(); ++ni) {
                     tbox::plog << "  " << *ni;
                  }
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
         id_generator[1 + d_rank_tree->getDegree()],
         mstream);

      unassigned.insertAll( my_subtree.d_work_traded );
      my_subtree.d_subtree_load_current += my_subtree.d_work_traded.getSumLoad();
      my_subtree.d_eff_load_current += my_subtree.d_work_traded.getSumLoad();

      if ( unassigned_highwater < unassigned.size() ) {
         unassigned_highwater = unassigned.size();
      }

      if (d_print_steps) {
         tbox::plog << "Received from parent " << d_rank_tree->getParentRank() << ":"
                    << my_subtree.d_work_traded.size() << " boxes ("
                    << my_subtree.d_work_traded.getSumLoad() << " units):";
         if ( my_subtree.d_work_traded.size() < 10 ) {
            for ( BoxTransitSet::const_iterator ni=my_subtree.d_work_traded.begin();
                  ni!=my_subtree.d_work_traded.end(); ++ni ) {
               const BoxTransitSet::BoxInTransit& box_in_transit = *ni;
               tbox::plog << "  " << box_in_transit;
            }
         }
         tbox::plog << std::endl;
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

            const LoadType export_load_actual =
               recip_subtree.d_work_traded.adjustLoad(
                  unassigned,
                  next_available_index[d_rank_tree->getDegree()],
                  id_generator[d_rank_tree->getDegree()],
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
               if ( recip_subtree.d_work_traded.size() < 10 ) {
                  for (BoxTransitSet::const_iterator ni = recip_subtree.d_work_traded.begin();
                       ni != recip_subtree.d_work_traded.end(); ++ni) {
                     tbox::plog << "  " << *ni;
                  }
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
      if ( unassigned.size() < 10 ) {
         for ( BoxTransitSet::const_iterator bi=unassigned.begin();
               bi!=unassigned.end(); ++bi ) {
            tbox::plog << "    " << *bi << '\n';
         }
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

   t_local_balancing->start();


#if 1
   assignBoxesToLocalProcess(
      balanced_box_level,
      balanced_to_unbalanced,
      unbalanced_to_balanced,
      unassigned );

   removeLocallyOriginatedBoxesFromBoxTransitSet( unassigned, d_mpi.getRank() );
#else
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
   for (BoxTransitSet::iterator
        ni = unassigned.begin();
        ni != unassigned.end(); /* incremented in loop */) {

      const BoxTransitSet::BoxInTransit& box_in_transit = *ni;
      balanced_box_level.addBox(box_in_transit.d_box);

      if (box_in_transit.d_box.isIdEqual(box_in_transit.d_orig_box)) {
         // Unchanged box implies assigned back to local process.  It requires no mapping.
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
#endif

   t_local_balancing->stop();

   t_load_distribution->stop();


   /*
    * Finish messages before starting edge info exchange.
    * We have only sends to complete, so it should not take
    * long to advance them all to completion.
    */
   if ( d_print_steps ) {
      tbox::plog << "TreeLoadBalancer::loadBalanceWithinRankGroup waiting for sends to complete.\n";
   }
   t_finish_sends->start();
   child_send_stage.advanceAll();
   parent_send_stage.advanceAll();
   t_finish_sends->stop();
   if ( d_print_steps ) {
      tbox::plog << "TreeLoadBalancer::loadBalanceWithinRankGroup completed sends.\n";
   }
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
   if ( d_print_steps ) {
      tbox::plog << "TreeLoadBalancer::loadBalanceWithinRankGroup for completed sends.\n";
   }


   if ( d_print_steps ) {
      tbox::plog << "TreeLoadBalancer::loadBalanceWithinRankGroup constructing unbalanced->balanced.\n";
   }
   constructSemilocalUnbalancedToBalanced(
      unbalanced_to_balanced,
      unassigned );
   if ( d_print_steps ) {
      tbox::plog << "TreeLoadBalancer::loadBalanceWithinRankGroup finished constructing unbalanced->balanced.\n";
   }

   t_get_map->stop();


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


   if (balance_to_anchor && balance_to_anchor->hasTranspose()) {
      t_use_map->barrierAndStart();
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix(d_object_name);
      mca.modify(
         balance_to_anchor->getTranspose(),
         unbalanced_to_balanced,
         &balance_box_level,
         &balanced_box_level);
      t_use_map->barrierAndStop();
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
      d_comm_graph_writer->addRecord( d_mpi, size_t(10), size_t(7) );

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
         "origins up",
         double(my_subtree.d_wants_work_from_parent ? 0 : my_subtree.d_work_traded.getNumberOfOriginatingProcesses()),
         tbox::CommGraphWriter::TO,
         prank );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(3),
         "load down",
         double(my_subtree.d_wants_work_from_parent ? my_subtree.d_work_traded.getSumLoad() : 0),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(4),
         "boxes down",
         double(my_subtree.d_wants_work_from_parent ? my_subtree.d_work_traded.size() : 0),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(5),
         "origins down",
         double(my_subtree.d_wants_work_from_parent ? my_subtree.d_work_traded.getNumberOfOriginatingProcesses() : 0),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(6),
         "bytes down",
         double(my_subtree.d_wants_work_from_parent ? parent_recv->getRecvSize() : int(0)),
         tbox::CommGraphWriter::FROM,
         (my_subtree.d_wants_work_from_parent ? prank : -1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(7),
         "child wait",
         t_child_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         d_rank_tree->getChildRank(0) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(8),
         "child wait",
         t_child_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         d_rank_tree->getChildRank(1) );

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(9),
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
 * Assign boxes to local process (put them in the balanced_box_level
 * and put edges in balanced<==>unbalanced Connector.
 *
 * We can generate balanced--->unbalanced edges for all unassigned
 * boxes because we have their origin info.  If the unassigned box
 * originated locally, we can generate the unbalanced--->balanced
 * edge for them as well.  However, we can't generate these edges
 * for boxes originating remotely, so these edges will be missing.
 */
void
TreeLoadBalancer::assignBoxesToLocalProcess(
   hier::BoxLevel& balanced_box_level,
   hier::Connector &balanced_to_unbalanced,
   hier::Connector &unbalanced_to_balanced,
   const BoxTransitSet& unassigned ) const
{
   /*
    * All unassigned boxes should go into balanced_box_level.  Put
    * them there and generate relationships in balanced<==>unbalanced
    * mapping Connectors where required.
    */

   for (BoxTransitSet::iterator ni = unassigned.begin();
        ni != unassigned.end(); ++ni ) {

      const BoxTransitSet::BoxInTransit& box_in_transit = *ni;
      balanced_box_level.addBox(box_in_transit.d_box);

      if (!box_in_transit.d_box.isIdEqual(box_in_transit.d_orig_box)) {

         balanced_to_unbalanced.insertLocalNeighbor(
            box_in_transit.d_orig_box,
            box_in_transit.d_box.getBoxId());

         if (box_in_transit.d_orig_box.getOwnerRank() == d_mpi.getRank()) {
            unbalanced_to_balanced.insertLocalNeighbor(
               box_in_transit.d_box,
               box_in_transit.d_orig_box.getBoxId());
         }
      }

   }

}




/*
 *************************************************************************
 * Remove local boxes from a BoxTransitSet.
 *************************************************************************
 */
void
TreeLoadBalancer::removeLocallyOriginatedBoxesFromBoxTransitSet(
   BoxTransitSet& transit_set,
   int local_rank ) const
{
   for (BoxTransitSet::iterator ni = transit_set.begin();
        ni != transit_set.end(); /* incremented in loop */) {
      if (ni->d_orig_box.getOwnerRank() == local_rank) {
         TBOX_ASSERT(ni->d_orig_box.getOwnerRank() == ni->d_box.getOwnerRank() );
         transit_set.erase(ni++);
      }
      else {
         ++ni;
      }
   }
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
   const BoxTransitSet& for_export = subtree_data.d_work_traded;
   msg << static_cast<int>(for_export.size());
   for (BoxTransitSet::const_iterator
        ni = for_export.begin(); ni != for_export.end(); ++ni) {
      const BoxTransitSet::BoxInTransit& box_in_transit = *ni;
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
   hier::SequentialLocalIdGenerator& id_generator,
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
    * As we pull each BoxTransitSet::BoxInTransit out, give it a new id that reflects
    * its new owner.
    */
   BoxTransitSet::BoxInTransit received_box(d_dim);
   for (int i = 0; i < num_boxes; ++i) {
      received_box.getFromMessageStream(msg);
      BoxTransitSet::BoxInTransit renamed_box(received_box,
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
   const BoxTransitSet& for_export = subtree_data.d_work_traded;
   msg << static_cast<int>(for_export.size());
   for (BoxTransitSet::const_iterator
        ni = for_export.begin(); ni != for_export.end(); ++ni) {
      const BoxTransitSet::BoxInTransit& box_in_transit = *ni;
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
   hier::SequentialLocalIdGenerator& id_generator,
   tbox::MessageStream &msg ) const
{
   t_unpack_load->start();
   int num_boxes = 0;
   msg >> num_boxes;
   /*
    * As we pull each BoxTransitSet::BoxInTransit out, give it a new id that reflects
    * its new owner.
    */
   BoxTransitSet::BoxInTransit received_box(d_dim);
   for (int i = 0; i < num_boxes; ++i) {
      received_box.getFromMessageStream(msg);
      BoxTransitSet::BoxInTransit renamed_box(received_box,
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
   const BoxTransitSet &unassigned,
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
 * BoxTransitSet::BoxInTransit back to the owners of the unbalanced Boxes that
 * originated them.  We don't know what ranks will send back the
 * balanced boxes, so we keep receiving messages from any rank until
 * we have accounted for all the cells in the unbalanced BoxLevel.
 *************************************************************************
 */
void
TreeLoadBalancer::constructSemilocalUnbalancedToBalanced(
   hier::MappingConnector &unbalanced_to_balanced,
   const BoxTransitSet &kept_imports ) const
{
   t_construct_semilocal->start();

   // Stuff the imported BoxTransitSet::BoxInTransits into buffers by their original owners.
   t_pack_edge->start();
   std::map<int,boost::shared_ptr<tbox::MessageStream> > outgoing_messages;
   for ( BoxTransitSet::const_iterator bi=kept_imports.begin();
         bi!=kept_imports.end(); ++bi ) {
      const BoxTransitSet::BoxInTransit &bit = *bi;
      boost::shared_ptr<tbox::MessageStream> &mstream =
         outgoing_messages[bit.d_orig_box.getOwnerRank()];
      if ( !mstream ) {
         mstream.reset(new tbox::MessageStream);
      }
      bit.putToMessageStream(*mstream);
   }
   t_pack_edge->stop();


   /*
    * The incoming unbalanced boxes need a mapping to describe their
    * change, but we don't know what they will become, so create empty
    * maps for now.  Should any not change, we'll erase their
    * neighborhood later.
    */
   for ( hier::BoxContainer::const_iterator bi=unbalanced_to_balanced.getBase().getBoxes().begin();
         bi!=unbalanced_to_balanced.getBase().getBoxes().end(); ++bi ) {
      if ( unbalanced_to_balanced.getHead().getBoxes().find(*bi) ==
           unbalanced_to_balanced.getHead().getBoxes().end() ) {
         unbalanced_to_balanced.makeEmptyLocalNeighborhood(bi->getBoxId());
      }
   }


   /*
    * Send outgoing_messages.  Optimization for mitigating contention:
    * Start by sending to the first recipient with a rank higher than
    * the local rank.
    */

   t_misc1->start();
   std::map<int,boost::shared_ptr<tbox::MessageStream> >::iterator recip_itr =
      outgoing_messages.upper_bound(d_mpi.getRank());
   if ( recip_itr == outgoing_messages.end() ) {
      recip_itr = outgoing_messages.begin();
   }

   int outgoing_messages_size = static_cast<int>(outgoing_messages.size());
   std::vector<tbox::SAMRAI_MPI::Request>
      send_requests( outgoing_messages_size, MPI_REQUEST_NULL );
   t_misc1->stop();

   if ( d_print_edge_steps ) {
      tbox::plog << "TreeLoadBalancer::constructSemilocalUnbalancedToBalanced: starting post-distribution barrier.\n";
   }
   t_post_load_distribution_barrier->start();
   d_mpi.Barrier(); // This barrier seems to speed up the load balancing, maybe by allowing one communication phase to finish before beginning another.
   t_post_load_distribution_barrier->stop();
   if ( d_print_edge_steps ) {
      tbox::plog << "TreeLoadBalancer::constructSemilocalUnbalancedToBalanced: finished post-distribution barrier.\n";
   }

   t_construct_semilocal_send_edges->start();
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
   t_construct_semilocal_send_edges->stop();


   /*
    * Determine number of cells in unbalanced that are not yet accounted
    * for in balanced.
    */
   t_construct_semilocal_local_accounting->start();
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
         // unbalanced_box has changed.  Parts of it may still be local.
         for ( hier::Connector::ConstNeighborIterator ni=unbalanced_to_balanced.begin(neighborhood_itr);
               ni!=unbalanced_to_balanced.end(neighborhood_itr); ++ni ) {
            TBOX_ASSERT( ni->getOwnerRank() == d_mpi.getRank() );
            num_unaccounted_cells -= ni->size();
         }

      }
      else {
         // unbalanced_box has not changed.  All of it is still local.
         num_unaccounted_cells -= unbalanced_box.size();
      }

   }
   if ( d_print_edge_steps ) {
      tbox::plog << num_unaccounted_cells << " unaccounted cells\n";
   }
   t_construct_semilocal_local_accounting->stop();


   /*
    * Receive info about exported cells from processes that now own
    * those cells.  Receive until all cells are accounted for.
    */

   t_misc2->start();
   std::vector<char> incoming_message; // Keep outside loop to avoid reconstructions.
   BoxTransitSet::BoxInTransit balanced_box_in_transit(d_dim);
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

      incoming_message.clear();
   }
   TBOX_ASSERT( num_unaccounted_cells == 0 );
   t_misc2->stop();


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
            sizeof(BoxTransitSet::BoxInTransit)*TreeLoadBalancer_FIRSTDATALEN);
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
         sizeof(BoxTransitSet::BoxInTransit)*TreeLoadBalancer_FIRSTDATALEN);

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
   const BoxTransitSet& a,
   const BoxTransitSet& b) const
{
   double surface_penalty = 0;
   for (BoxTransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->d_box);
   }
   for (BoxTransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
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
   const BoxTransitSet& a,
   const BoxTransitSet& b) const
{
   double slender_penalty = 0;
   for (BoxTransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->d_box);
   }
   for (BoxTransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
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
      d_print_pop_steps =
         input_db->getBoolWithDefault("DEV_print_pop_steps", d_print_pop_steps);
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

      d_allow_box_breaking =
         input_db->getBoolWithDefault("DEV_allow_box_breaking",
                                      d_allow_box_breaking);

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
      t_construct_semilocal_send_edges = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constructSemilocalUnbalancedToBalanced()_send_edges");
      t_construct_semilocal_local_accounting = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constructSemilocalUnbalancedToBalanced()_local_accounting");
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
      t_misc1 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::misc1");
      t_misc2 = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::misc2");
   }
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
