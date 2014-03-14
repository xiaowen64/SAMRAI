/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner_C
#define included_mesh_CascadePartitioner_C

#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/mesh/VoucherTransitLoad.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/hier/BoxContainer.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/tbox/PIO.h"
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

const int CascadePartitioner::CascadePartitioner_LOADTAG0;
const int CascadePartitioner::CascadePartitioner_LOADTAG1;
const int CascadePartitioner::CascadePartitioner_FIRSTDATALEN;

const int CascadePartitioner::s_default_data_id = -1;


/*
 *************************************************************************
 * CascadePartitioner constructor.
 *************************************************************************
 */

CascadePartitioner::CascadePartitioner(
   const tbox::Dimension& dim,
   const std::string& name,
   const boost::shared_ptr<tbox::Database>& input_db,
   const boost::shared_ptr<tbox::RankTreeStrategy> &rank_tree):
   d_dim(dim),
   d_object_name(name),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mpi_is_dupe(false),
   d_master_workload_data_id(s_default_data_id),
   d_flexible_load_tol(0.05),
   d_mca(),
   // Shared data.
   d_comm_stage(),
   // Performance evaluation.
   d_barrier_before(false),
   d_barrier_after(false),
   d_report_load_balance(false),
   d_summarize_map(false),
   d_print_steps(false),
   d_check_connectivity(false),
   d_check_map(false),
   d_num_ag_cycles(0)
{
   for ( int i=0; i<4; ++i ) d_comm_peer[i].initialize(&d_comm_stage);

   TBOX_ASSERT(!name.empty());
   getFromInput(input_db);
   setTimers();
   d_mca.setTimerPrefix(d_object_name);
}



/*
 *************************************************************************
 * CascadePartitioner constructor.
 *************************************************************************
 */

CascadePartitioner::~CascadePartitioner()
{
   freeMPICommunicator();
}



/*
 *************************************************************************
 * Accessory functions to get/set load balancing parameters.
 *************************************************************************
 */

bool
CascadePartitioner::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}



/*
**************************************************************************
**************************************************************************
*/
void
CascadePartitioner::setWorkloadPatchDataIndex(
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
 * This method implements the abstract LoadBalanceStrategy interface.
 *************************************************************************
 */
void
CascadePartitioner::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_reference,
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
   NULL_USE(domain_box_level);
   TBOX_ASSERT(!balance_to_reference || balance_to_reference->hasTranspose());
   TBOX_ASSERT(!balance_to_reference || balance_to_reference->isTransposeOf(balance_to_reference->getTranspose()));
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
            TBOX_ERROR("CascadePartitioner::loadBalanceBoxLevel:\n"
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

   if (d_print_steps) {
      tbox::plog << "CascadePartitioner::loadBalanceBoxLevel called with:"
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
    * reference<==>balance.
    */
   balance_box_level.removePeriodicImageBoxes();
   if (balance_to_reference) {
      balance_to_reference->getTranspose().removePeriodicRelationships();
      balance_to_reference->getTranspose().setHead(balance_box_level, true);
      balance_to_reference->removePeriodicRelationships();
      balance_to_reference->setBase(balance_box_level, true);
   }

   t_load_balance_box_level->start();

   d_pparams = boost::make_shared<PartitioningParams>(
      *balance_box_level.getGridGeometry(),
      balance_box_level.getRefinementRatio(),
      min_size, max_size, bad_interval, cut_factor,
      d_flexible_load_tol);

   LoadType local_load = computeLocalLoad(balance_box_level);

   LoadType global_sum_load = local_load;


   /*
    * Determine the total load and number of processes that has any
    * initial load.
    */
   global_sum_load = local_load;
   if (d_mpi.getSize() > 1) {
      d_mpi.AllReduce( &global_sum_load, 1, MPI_SUM );
   }

   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << "CascadePartitioner::loadBalanceBoxLevel"
                 << " global_sum_load=" << global_sum_load
                 << " (initially born across " << d_mpi.getSize()
                 << " procs, averaging " << global_sum_load / d_mpi.getSize()
                 << " or " << pow(global_sum_load / d_mpi.getSize(), 1.0 / d_dim.getValue())
                 << "^" << d_dim << " per proc." << std::endl;
   }


   d_global_load_avg = global_sum_load / rank_group.size();

   // Run the partitioning algorithm.
   partitionByCascade(
      balance_box_level,
      balance_to_reference,
      global_sum_load );


   /*
    * Finished load balancing.  Clean up and wrap up.
    */

   d_pparams.reset();

   t_load_balance_box_level->stop();

   local_load = computeLocalLoad(balance_box_level);
   d_load_stat.push_back(local_load);
   d_box_count_stat.push_back(
      static_cast<int>(balance_box_level.getBoxes().size()));

   if (d_print_steps) {
      tbox::plog << "Post balanced:\n" << balance_box_level.format("", 2);
   }

   if (d_report_load_balance) {
      tbox::plog
         << "CascadePartitioner::loadBalanceBoxLevel results  ";
      BalanceUtilities::gatherAndReportLoadBalance(local_load,
         balance_box_level.getMPI());
   }

   if (d_check_connectivity && balance_to_reference) {
      hier::Connector& reference_to_balance = balance_to_reference->getTranspose();
      tbox::plog << "CascadePartitioner checking balance-reference connectivity."
                 << std::endl;
      int errs = 0;
      if (reference_to_balance.checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in reference_to_balance!\n";
      }
      if (balance_to_reference->checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in balance_to_reference!\n";
      }
      if (reference_to_balance.checkTransposeCorrectness(*balance_to_reference)) {
         ++errs;
         tbox::perr << "Error found in balance-reference transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found.\n"
            << "reference_box_level:\n" << reference_to_balance.getBase().format("", 2)
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "reference_to_balance:\n" << reference_to_balance.format("", 2)
            << "balance_to_reference:\n" << balance_to_reference->format("", 2));
      }
      tbox::plog << "CascadePartitioner checked balance-reference connectivity."
                 << std::endl;
   }

}



/*
 *************************************************************************
 * This method implements the cascade partitioner algorithm.
 *
 * It calls distributeLoad to do distribute the load then assigns the
 * distributed the loads to their new owners and update Connectors.
 *************************************************************************
 */
void
CascadePartitioner::partitionByCascade(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_reference,
   double global_sum_load
   ) const
{
   if ( d_print_steps ) {
      tbox::plog << "CascadePartitioner::partitionByCascade: entered" << std::endl;
   }

   BoxTransitSet local_work(*d_pparams), shipment(*d_pparams);
   local_work.insertAll( balance_box_level.getBoxes() );
   d_local_load = &local_work;
   d_shipment = &shipment;

   /*
    * Initialize empty balanced_box_level and mappings so they are
    * ready to be populated.
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


   double group_load = local_work.getSumLoad();

   const int lg_size = lgInt(d_mpi.getSize());

   const int ag_group_size = 1 << d_num_ag_cycles;

   t_get_map->start();

#if 0
   if ( d_mpi.getSize() >= ag_group_size && ag_group_size > 1 ) {
      agglomerate( ag_group_size, minMessageBytes );
   }
#endif


tbox::plog << "lg_size="<<lg_size << "  ag_group_size="<<ag_group_size << std::endl;
   CascadePartitionerTree groups(*this);
   groups.balanceAll();


#if 0
   if ( d_mpi.getSize() >= ag_group_size && ag_group_size > 1 ) {
      deagglomerate( ag_group_size, minMessageBytes );
   }
#endif

   if ( d_print_steps ) {
      tbox::plog << "CascadePartitioner: local_work after load distribution:\n";
      local_work.recursivePrint(tbox::plog, "LL->\t", 2);
      tbox::plog
         << "CascadePartitioner::partitionByCascade constructing unbalanced<==>balanced.\n";
   }

   t_assign_to_local_and_populate_maps->start();
   local_work.assignToLocalAndPopulateMaps(
      balanced_box_level,
      balanced_to_unbalanced,
      unbalanced_to_balanced,
      d_flexible_load_tol );
   t_assign_to_local_and_populate_maps->stop();

   t_get_map->stop();

   if ( d_summarize_map ) {
      tbox::plog << "CascadePartitioner::partitionByCascade unbalanced--->balanced map:\n"
                 << unbalanced_to_balanced.format("\t",0)
                 << "Map statistics:\n" << unbalanced_to_balanced.formatStatistics("\t")
                 << "CascadePartitioner::partitionByCascade balanced--->unbalanced map:\n"
                 << balanced_to_unbalanced.format("\t",0)
                 << "Map statistics:\n" << balanced_to_unbalanced.formatStatistics("\t")
                 << '\n';
   }

   if (d_check_map) {
      if (unbalanced_to_balanced.findMappingErrors() != 0) {
         TBOX_ERROR(
            "CascadePartitioner::partitionByCascade Mapping errors found in unbalanced_to_balanced!");
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         TBOX_ERROR(
            "CascadePartitioner::partitionByCascade Transpose errors found!");
      }
   }


   if ( d_summarize_map ) {
      tbox::plog << "CascadePartitioner::partitionByCascade: unbalanced--->balanced map:\n"
                 << unbalanced_to_balanced.format("\t",0)
                 << "Map statistics:\n" << unbalanced_to_balanced.formatStatistics("\t")
                 << "CascadePartitioner::partitionByCascade: balanced--->unbalanced map:\n"
                 << balanced_to_unbalanced.format("\t",0)
                 << "Map statistics:\n" << balanced_to_unbalanced.formatStatistics("\t")
                 << '\n';
   }


   if (balance_to_reference && balance_to_reference->hasTranspose()) {
      t_use_map->barrierAndStart();
      d_mca.modify(
         balance_to_reference->getTranspose(),
         unbalanced_to_balanced,
         &balance_box_level,
         &balanced_box_level);
      t_use_map->barrierAndStop();
   } else {
      hier::BoxLevel::swap(balance_box_level, balanced_box_level);
   }

   if ( d_print_steps ) {
      tbox::plog
         << "CascadePartitioner::partitionByCascade finished constructing unbalanced<==>balanced.\n";
   }

   d_local_load = 0;

   if ( d_print_steps ) {
      tbox::plog << "CascadePartitioner::partitionByCascade: leaving" << std::endl;
   }
}



/*
 *************************************************************************
 * Compute log-base-2 of integer, rounded up.
 *************************************************************************
 */
int CascadePartitioner::lgInt(int s) {
   int lg_s = 0;
   while ( (1<<lg_s) < s ) {
      ++lg_s;
   }
   return lg_s;
}



/*
 *************************************************************************
 * Set the MPI commuicator.  If there's a private communicator, free
 * it first.  It's safe to free the private communicator because no
 * other code have access to it.
 *************************************************************************
 */
void
CascadePartitioner::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& samrai_mpi)
{
   if (samrai_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      TBOX_ERROR("CascadePartitioner::setSAMRAI_MPI error: Given\n"
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
CascadePartitioner::freeMPICommunicator()
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
 *************************************************************************
 */
CascadePartitioner::LoadType
CascadePartitioner::computeLocalLoad(
   const hier::BoxLevel& box_level) const
{
   double load = 0.0;
   const hier::BoxContainer& boxes = box_level.getBoxes();
   for (hier::BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end();
        ++ni) {
      double box_load = double(ni->size());
      load += box_load;
   }
   return static_cast<LoadType>(load);
}



/*
 *************************************************************************
 *
 * Read values (described in the class header) from input database.
 *
 *************************************************************************
 */

void
CascadePartitioner::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{

   if (input_db) {

      d_print_steps = input_db->getBoolWithDefault("DEV_print_steps", false);
      d_check_connectivity =
         input_db->getBoolWithDefault("DEV_check_connectivity", d_check_connectivity);
      d_check_map =
         input_db->getBoolWithDefault("DEV_check_map", d_check_map);

      d_summarize_map = input_db->getBoolWithDefault("DEV_summarize_map",
         d_summarize_map);

      d_report_load_balance = input_db->getBoolWithDefault(
         "DEV_report_load_balance", d_report_load_balance);
      d_barrier_before = input_db->getBoolWithDefault("DEV_barrier_before",
         d_barrier_before);
      d_barrier_after = input_db->getBoolWithDefault("DEV_barrier_after",
         d_barrier_after);

      d_flexible_load_tol =
         input_db->getDoubleWithDefault("flexible_load_tolerance",
            d_flexible_load_tol);

   }
}



/*
 ***************************************************************************
 *
 ***************************************************************************
 */
void
CascadePartitioner::assertNoMessageForPrivateCommunicator() const
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
            << "CascadePartitioner detected before or\n"
            << "after using a private communicator that there\n"
            << "is a message yet to be received.  This is\n"
            << "an error because all messages using the\n"
            << "private communicator should have been\n"
            << "accounted for.  Message status:\n"
            << "source " << mpi_status.MPI_SOURCE << '\n'
            << "tag " << mpi_status.MPI_TAG << '\n'
            << "count " << count << " (assuming integers)\n"
            << "current tags: "
            << ' ' << CascadePartitioner_LOADTAG0 << ' '
            << CascadePartitioner_LOADTAG1
            );
      }
   }
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
CascadePartitioner::setTimers()
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

      t_assign_to_local_and_populate_maps = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::assign_to_local_and_populate_maps");

   }
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
CascadePartitioner::printStatistics(
   std::ostream& output_stream) const
{
   BalanceUtilities::gatherAndReportLoadBalance(
      d_load_stat,
      tbox::SAMRAI_MPI::getSAMRAIWorld(),
      output_stream);
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
