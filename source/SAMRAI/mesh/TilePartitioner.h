/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer for fixed-size, tiled patches.
 *
 ************************************************************************/

#ifndef included_mesh_TilePartitioner
#define included_mesh_TilePartitioner

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RankGroup.h"

#include "boost/shared_ptr.hpp"
#include <iostream>
#include <vector>

namespace SAMRAI {
namespace mesh {





/*!
 * @brief Load balance by cutting boxes along designated mesh lines
 * defining a tile box grid.  Implemementing the LoadBalancerStrategy.
 *
 * This class delegates its partitioning to the TreeLoadBalancer
 * but overrides the cut_factor parameter in loadBalanceBoxLevel()
 * to control the patch size.
 *
 * @verbatim
 *
 * User inputs (default):
 *
 * - IntVector @b box_size Box size in the index space of the tag level.
 *
 * - report_load_balance = TRUE // Write out load balance report in log.
 *
 * // Internal load balancer to use: TreeLoadBalancer, ChopAndPackLoadBalancer
 * // TODO: There may be a bug in the ChopAndPackLoadBalancer.
 * internal_load_balancer = "TreeLoadBalancer"
 * TreeLoadBalancer { ... } // Input for internal TreeLoadBalancer.
 * ChopAndPackLoadBalancer { ... } // Input for internal ChopAndPackLoadBalancer.
 *
 * @endverbatim
 *
 * @see mesh::LoadBalanceStrategy
 */

class TilePartitioner:
   public LoadBalanceStrategy
{
public:
   /*!
    * @brief Initializing constructor sets object state to default or,
    * if database provided, to parameters in database.
    *
    * @param[in] dim
    *
    * @param[in] name User-defined std::string identifier used for error
    * reporting and timer names.  If omitted, "TilePartitioner"
    * is used.
    *
    * @param[in] rank_tree How to arange a contiguous range of MPI ranks
    * into a tree.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    */
   TilePartitioner(
      const tbox::Dimension& dim,
      const std::string& name,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~TilePartitioner();

   /*!
    * @brief Return true if load balancing procedure for given level
    * depends on patch data on mesh; otherwise return false.
    *
    * @param[in] level_number  Integer patch level number.
    */
   bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const { return false; }

   /*!
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
    *
    * Note: This implementation does not yet support non-uniform load
    * balancing.
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector *balance_to_anchor,
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const;

   /*!
    * @brief Write out statistics recorded for the most recent load
    * balancing result.
    *
    * @param[in] output_stream
    */
   void
   printStatistics(
      std::ostream& output_stream = tbox::plog) const
   {
      BalanceUtilities::gatherAndReportLoadBalance(d_load_stat,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         output_stream);
   }



private:

   typedef double LoadType;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   TilePartitioner(
      const TilePartitioner&);

   void
   operator = (
      const TilePartitioner&);


   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& input_db);

   /*
    * Object dimension.
    */
   const tbox::Dimension d_dim;

   /*
    * String identifier for load balancer object.
    */
   std::string d_object_name;

   ChopAndPackLoadBalancer d_cap;
   TreeLoadBalancer d_tlb;

   /*
    * @brief Which internal load balancer to use.
    *
    * 'c' = ChopAndPackLoadBalancer
    * 't' = TreeLoadBalancer
    */
   char d_internal_load_balancer;

   hier::IntVector d_box_size;

   /*!
    * @brief Whether to immediately report the results of the load balancing cycles
    * in the log files.
    */
   bool d_report_load_balance;

   //@{
   //! @name Used for evaluating peformance.
   bool d_barrier_and_time;
   boost::shared_ptr<tbox::Timer> t_load_balance_box_level;
   //@}

   /*
    * Statistics on number of cells and patches generated.
    */
   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;


   // Extra debug independent of optimization/debug flag.
   char d_print_steps;

};

}
}

#endif
