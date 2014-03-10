/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner
#define included_mesh_CascadePartitioner

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/CascadePartitionerGroup.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/RankTreeStrategy.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/shared_ptr.hpp"
#include <iostream>
#include <vector>

namespace SAMRAI {
namespace mesh {



/*!
 * @brief Provides load balancing routines for AMR hierarchy by
 * implemementing the LoadBalancerStrategy.
 *
 * Currently, only uniform load balancing is supported.  Eventually,
 * non-uniform load balancing should be supported.  (Non-uniform load
 * balancing is supported by the CutAndPackLoadBalancer class.)
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *
 *   - \b flexible_load_tolerance
 *   Fraction of ideal load a process can
 *   take on in order to reduce box cutting and load movement.  Higher
 *   values often reduce partitioning time and box count but produce
 *   less balanced work loads.  Surplus work greater than this
 *   tolerance can still result due to other constraints, such as
 *   minimum box size.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>flexible_load_tolerance</td>
 *     <td>double</td>
 *     <td>0.05</td>
 *     <td>0-1</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @internal The following are developer inputs.  Defaults listed
 * in parenthesis:
 *
 * @internal DEV_voucher_mode (false)
 * bool
 * Whether to use experimental voucher mode.
 *
 * @internal DEV_allow_box_breaking (true)
 * bool
 * Whether to allow box-breaking.  Set to false when boxes have
 * been pre-cut.
 *
 * @see mesh::LoadBalanceStrategy
 */

class CascadePartitioner:
   public LoadBalanceStrategy
{
public:
   /*!
    * @brief Initializing constructor sets object state to default or,
    * if database provided, to parameters in database.
    *
    * @param[in] dim
    *
    * @param[in] name User-defined identifier used for error reporting
    * and timer names.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    *
    * @param[in] rank_tree How to arange a contiguous range of MPI ranks
    * into a tree.  If omitted, we use a tbox::CenteredRankTree.
    *
    * @pre !name.empty()
    */
   CascadePartitioner(
      const tbox::Dimension& dim,
      const std::string& name,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>(),
      const boost::shared_ptr<tbox::RankTreeStrategy> &rank_tree =
         boost::shared_ptr<tbox::RankTreeStrategy>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~CascadePartitioner();

   /*!
    * @brief Set the internal SAMRAI_MPI to a duplicate of the given
    * SAMRAI_MPI.
    *
    * The given SAMRAI_MPI must have a valid communicator.
    *
    * The given SAMRAI_MPI is duplicated for private use.  This
    * requires a global communication, so all processes in the
    * communicator must call it.  The advantage of a duplicate
    * communicator is that it ensures the communications for the
    * object won't accidentally interact with unrelated
    * communications.
    *
    * If the duplicate SAMRAI_MPI it is set, the CascadePartitioner will
    * only balance BoxLevels with congruent SAMRAI_MPI objects and
    * will use the duplicate SAMRAI_MPI for communications.
    * Otherwise, the SAMRAI_MPI of the BoxLevel will be used.  The
    * duplicate MPI communicator is freed when the object is
    * destructed, or freeMPICommunicator() is called.
    *
    * @pre samrai_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull
    */
   void
   setSAMRAI_MPI(
      const tbox::SAMRAI_MPI& samrai_mpi);

   /*!
    * @brief Free the internal MPI communicator, if any has been set.
    *
    * This is automatically done by the destructor, if needed.
    *
    * @see setSAMRAI_MPI().
    */
   void
   freeMPICommunicator();

   /*!
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
    *
    * Note: This implementation does not yet support non-uniform load
    * balancing.
    *
    * @pre !balance_to_reference || balance_to_reference->hasTranspose()
    * @pre !balance_to_reference || balance_to_reference->isTransposeOf(balance_to_reference->getTranspose())
    * @pre (d_dim == balance_box_level.getDim()) &&
    *      (d_dim == min_size.getDim()) && (d_dim == max_size.getDim()) &&
    *      (d_dim == domain_box_level.getDim()) &&
    *      (d_dim == bad_interval.getDim()) && (d_dim == cut_factor.getDim())
    * @pre !hierarchy || (d_dim == hierarchy->getDim())
    * @pre !d_mpi_is_dupe || (d_mpi.getSize() == balance_box_level.getMPI().getSize())
    * @pre !d_mpi_is_dupe || (d_mpi.getSize() == balance_box_level.getMPI().getRank())
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_reference,
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
      std::ostream& output_stream = tbox::plog) const;


   /*!
    * @brief Get the name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

   /*!
    * @brief Configure the load balancer to use the data stored
    * in the hierarchy at the specified descriptor index
    * for estimating the workload on each cell.
    *
    * Note: This method currently does not affect the results because
    * this class does not yet support uniform load balancing.
    *
    * @param data_id
    * Integer value of patch data identifier for workload
    * estimate on each cell.  An invalid value (i.e., < 0)
    * indicates that a spatially-uniform work estimate
    * will be used.  The default value is -1 (undefined)
    * implying the uniform work estimate.
    *
    * @param level_number
    * Optional integer number for level on which data id
    * is used.  If no value is given, the data will be
    * used for all levels.
    *
    * @pre hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataFactory(data_id) is actually a  boost::shared_ptr<pdat::CellDataFactory<double> >
    */
   void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1);

   /*!
    * @brief Return true if load balancing procedure for given level
    * depends on patch data on mesh; otherwise return false.
    *
    * @param[in] level_number  Integer patch level number.
    */
   bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const;

private:

   typedef double LoadType;


   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int CascadePartitioner_LOADTAG0 = 1;
   static const int CascadePartitioner_LOADTAG1 = 2;
   static const int CascadePartitioner_FIRSTDATALEN = 500;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   CascadePartitioner(
      const CascadePartitioner&);

   void
   operator = (
      const CascadePartitioner&);


   /*
    * @brief Check if there is any pending messages for the private
    * communication and throw an error if there is.
    */
   void
   assertNoMessageForPrivateCommunicator() const;

   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& input_db);

   /*
    * Utility functions to determine parameter values for level.
    */
   int
   getWorkloadDataId(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return (level_number < static_cast<int>(d_workload_data_id.size()) ?
         d_workload_data_id[level_number] :
         d_master_workload_data_id);
   }

   /*
    * Count the local workload.
    */
   LoadType
   computeLocalLoad(
      const hier::BoxLevel& box_level) const;

   /*!
    *@brief Implements the cascade partitioner algorithm.
    */
   void
   partitionByCascade(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_reference,
      double global_sum_load ) const;

   /*!
    * @brief Executes a single cascade cycle.
    */
   void
   cascadeCycle(
      double &group_surplus,
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_reference,
      int outerCycle,
      int innerCycle ) const;

   //! @brief Compute log-base-2 of integer, rounded up.
   int lgInt(int s) const;

   /*!
    * @brief Set up timers for the object.
    */
   void
   setTimers();

   /*
    * CascadePartitioner and CascadePartitionerGroup are tightly
    * coupled.  CascadePartitioner has the common parts of the data
    * and algorithm.  CascadePartitionerGroup has the group-specific
    * parts.  CascadePartitionerGroup can be made a private subclass
    * of CascadePartitioner, but that would make a big file.
    */
   friend CascadePartitionerGroup;

   /*
    * Object dimension.
    */
   const tbox::Dimension d_dim;

   /*
    * String identifier for load balancer object.
    */
   std::string d_object_name;

   //! @brief Duplicated communicator object.  See setSAMRAI_MPI().
   mutable tbox::SAMRAI_MPI d_mpi;

   //! @brief Whether d_mpi is an internal duplicate.  See setSAMRAI_MPI().
   bool d_mpi_is_dupe;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   std::vector<int> d_workload_data_id;

   int d_master_workload_data_id;

   /*!
    * @brief Fraction of ideal load a process can accept over and above
    * the ideal.
    *
    * See input parameter "flexible_load_tolerance".
    */
   double d_flexible_load_tol;

   /*!
    * @brief Metadata operations with timers set according to this object.
    */
   hier::MappingConnectorAlgorithm d_mca;

   //@{
   //! @name Data shared with private methods during balancing.
   mutable boost::shared_ptr<PartitioningParams> d_pparams;
   mutable LoadType d_global_load_avg;
   mutable LoadType d_min_load;
   //@}

   static const int s_default_data_id;

   int d_num_ag_cycles;

   //@{
   //! @name Used for evaluating peformance.

   bool d_barrier_before;
   bool d_barrier_after;

   /*!
    * @brief Whether to immediately report the results of the load
    * balancing cycles in the log files.
    */
   bool d_report_load_balance;

   /*!
    * @brief See "summarize_map" input parameter.
    */
   char d_summarize_map;

   /*
    * Performance timers.
    */
   boost::shared_ptr<tbox::Timer> t_load_balance_box_level;
   boost::shared_ptr<tbox::Timer> t_assign_to_local_and_populate_maps;
   boost::shared_ptr<tbox::Timer> t_get_map;
   boost::shared_ptr<tbox::Timer> t_use_map;

   //@}

   // Extra checks independent of optimization/debug.
   char d_print_steps;
   char d_check_connectivity;
   char d_check_map;

   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;

};

}
}

#endif
