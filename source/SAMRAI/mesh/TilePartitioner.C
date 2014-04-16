/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TilePartitioner_C
#define included_mesh_TilePartitioner_C

#include "SAMRAI/mesh/TilePartitioner.h"

#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {


/*
 *************************************************************************
 * TilePartitioner constructor.
 *************************************************************************
 */

TilePartitioner::TilePartitioner(
   const tbox::Dimension& dim,
   const std::string& name,
   const boost::shared_ptr<tbox::Database>& database):
   d_dim(dim),
   d_object_name(name),
   d_cap(dim, name + ":ChopAndPackLoadBalancer",
         database ? database->getDatabaseWithDefault("ChopAndPackLoadBalancer",
                                                     boost::shared_ptr<tbox::Database>())
         : boost::shared_ptr<tbox::Database>()),
   d_tlb(dim, name + ":TreeLoadBalancer",
         database ? database->getDatabaseWithDefault("TreeLoadBalancer",
                                                     boost::shared_ptr<tbox::Database>())
         : boost::shared_ptr<tbox::Database>()),
   d_cp(dim, name + ":CascadePartitioner",
         database ? database->getDatabaseWithDefault("CascadePartitioner",
                                                     boost::shared_ptr<tbox::Database>())
         : boost::shared_ptr<tbox::Database>()),
   d_graphlb(dim, name + ":GraphLoadBalancer",
         database ? database->getDatabaseWithDefault("GraphLoadBalancer",
                                                     boost::shared_ptr<tbox::Database>())
         : boost::shared_ptr<tbox::Database>()),
   d_internal_load_balancer('t'),
   d_box_size(dim),
   // Output control.
   d_report_load_balance(false),
   // Performance evaluation.
   d_barrier_and_time(false),
   d_print_steps(false)
{
   TBOX_ASSERT(!name.empty());
   getFromInput(database);

   t_load_balance_box_level = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::loadBalanceBoxLevel()");
}



/*
 *************************************************************************
 * TilePartitioner constructor.
 *************************************************************************
 */

TilePartitioner::~TilePartitioner()
{
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
TilePartitioner::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector *balance_to_anchor,
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


   if (d_report_load_balance) {
      tbox::plog
         << "TilePartitioner::loadBalanceBoxLevel loads before partitioning:"
         << std::endl;
      BalanceUtilities::gatherAndReportLoadBalance(
         static_cast<double>(balance_box_level.getLocalNumberOfCells()),
         balance_box_level.getMPI());
   }


   if (d_print_steps) {
      tbox::plog << "TilePartitioner::loadBalanceBoxLevel called with:"
                 << "\n  min_size = " << min_size
                 << "\n  max_size = " << max_size
                 << "\n  bad_interval = " << bad_interval
                 << "\n  cut_factor = " << cut_factor
                 << "\n  prebalance:\n"
                 << balance_box_level.format("  ", 2);
   }

   // cut_factor must be <= the box size.
   if ( !(cut_factor <= d_box_size) ) {
      TBOX_ERROR("TilePartitioner::loadBalanceBoxLevel:\n"
                 << " cut_factor " << cut_factor
                 << " must be <= box_size " << d_box_size << '\n');
   }

   // d_box_size must be divisible by cut_factor.
   if ( d_box_size/cut_factor*cut_factor != d_box_size ) {
      TBOX_ERROR("TilePartitioner::loadBalanceBoxLevel:\n"
                 << " box_size " << d_box_size
                 << " must be divisible by cut_factor " << cut_factor << '\n');
   }

   if (d_barrier_and_time) {
      t_load_balance_box_level->barrierAndStart();
   }

   hier::IntVector tile_cut_factor = d_box_size;

   switch (d_internal_load_balancer) {
   case 'd':
      d_cp.loadBalanceBoxLevel(
         balance_box_level,
         balance_to_anchor,
         hierarchy,
         level_number,
         min_size,
         max_size,
         domain_box_level,
         bad_interval,
         tile_cut_factor,
         rank_group);
      break;
   case 'g':
      d_graphlb.loadBalanceBoxLevel(
         balance_box_level,
         balance_to_anchor,
         hierarchy,
         level_number,
         min_size,
         max_size,
         domain_box_level,
         bad_interval,
         tile_cut_factor,
         rank_group);
      break;
   case 't':
      d_tlb.loadBalanceBoxLevel(
         balance_box_level,
         balance_to_anchor,
         hierarchy,
         level_number,
         min_size,
         max_size,
         domain_box_level,
         bad_interval,
         tile_cut_factor,
         rank_group);
      break;
   case 'c':
      /*
       * TODO: There may be a bug in the ChopAndPackLoadBalancer.
       * Its cuts can fall off the tile boundaries.
       */
      TBOX_WARNING("TilePartitioner: Warning: using the \"ChopAndPackLoadBalancer\"\n"
                   <<"internal load balancer may produce cuts that fall off of\n"
                   <<"the tile boundaries.  This is a known bug.");
      d_cap.loadBalanceBoxLevel(
         balance_box_level,
         balance_to_anchor,
         hierarchy,
         level_number,
         min_size,
         max_size,
         domain_box_level,
         bad_interval,
         tile_cut_factor,
         rank_group);
      break;
   default:
      TBOX_ERROR("TilePartitioner error: We should not be here!");
   }

   if (d_barrier_and_time) {
      t_load_balance_box_level->stop();
   }


   if (d_report_load_balance) {
      tbox::plog
         << "TilePartitioner::loadBalanceBoxLevel loads after partitioning:"
         << std::endl;
      BalanceUtilities::gatherAndReportLoadBalance(
         static_cast<double>(balance_box_level.getLocalNumberOfCells()),
         balance_box_level.getMPI());
   }


   // Save results for analysis of quality.
   d_load_stat.push_back(static_cast<double>(balance_box_level.getLocalNumberOfCells()));

}



/*
 *************************************************************************
 * Read values (described in the class header) from input database.
 *************************************************************************
 */

void
TilePartitioner::getFromInput(
   const boost::shared_ptr<tbox::Database>& database)
{

   if (database) {

      if (database->isInteger("box_size")) {
         database->getIntegerArray("box_size",
            &d_box_size[0],
            d_dim.getValue());
      }

      if ( database->isString("internal_load_balancer") ) {

         std::string internal_load_balancer =
            database->getString("internal_load_balancer");

         if ( internal_load_balancer == "ChopAndPackLoadBalancer" ) {
            d_internal_load_balancer = 'c';
         }
         else if ( internal_load_balancer == "CascadePartitioner" ) {
            d_internal_load_balancer = 'd';
         }
         else if ( internal_load_balancer == "TreeLoadBalancer" ) {
            d_internal_load_balancer = 't';
         }
         else if ( internal_load_balancer == "GraphLoadBalancer" ) {
            d_internal_load_balancer = 'g';
         }
         else {
            TBOX_ERROR("TilePartitioner::getFromInput error:\n"
                       <<"internal_load_balancer must be set to\n"
                       <<"\"ChopAndPackLoadBalancer\" or \"TreeLoadBalancer\" or \"CascadePartitioner\".\n"
                       <<"You specified \"" << internal_load_balancer << "\".");
         }

      }

      d_print_steps = database->getBoolWithDefault("DEV_print_steps", d_print_steps);

      d_report_load_balance = database->getBoolWithDefault(
         "DEV_report_load_balance",
         d_report_load_balance);

   }
}



/*
 *************************************************************************
 * Write out statistics recorded for the most recent load
 * balancing result.
 *************************************************************************
 */
void
TilePartitioner::printStatistics(
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
