/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer
#define included_mesh_TreeLoadBalancer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/CommGraphWriter.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/RankTreeStrategy.h"
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/shared_ptr.hpp"
#include <iostream>
#include <vector>
#include <set>

namespace SAMRAI {
namespace mesh {





/*!
 * @brief Provides load balancing routines for AMR hierarchy by
 * implemementing the LoadBalancerStrategy.
 *
 * This class implements a tree-based load balancer.  The MPI
 * processes are arranged in a tree.  Work load is transmitted from
 * process to process along the edges of the tree.
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
 *   Fraction of ideal load a process can take on in order to avoid excessive
 *   box cutting and load movement.  This is not a hard limit and some
 *   processes can still exceed this amount.  Higher values help the load
 *   balancer run faster but produces less balanced work loads.
 *
 *   - \b max_cycle_spread_ratio
 *   This parameter limits how many processes may receive the load of one
 *   process in a load fan-out cycle.  If a process has too much initial load,
 *   this limit causes the load to fan out the load over multiple cycles.  It
 *   alleviates the bottle-neck of one process having to work with too many
 *   other processes in any cycle.
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
 *     <td>0.0</td>
 *     <td>0-1</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>max_cycle_spread_ratio</td>
 *     <td>int</td>
 *     <td>1000000</td>
 *     <td> > 1</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @internal The following are developer inputs.  Defaults listed
 * in parenthesis:
 *
 * @internal DEV_allow_box_breaking (true)
 * bool
 * Whether to allow box-breaking.  Set to false when boxes have
 * been pre-cut.
 *
 * @see mesh::LoadBalanceStrategy
 */

class TreeLoadBalancer:
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
    * reporting and timer names.  If omitted, "TreeLoadBalancer"
    * is used.
    *
    * @param[in] rank_tree How to arange a contiguous range of MPI ranks
    * into a tree.  If omitted, we use a tbox::CenteredRankTree.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    *
    * @pre !name.empty()
    */
   TreeLoadBalancer(
      const tbox::Dimension& dim,
      const std::string& name,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>(),
      const boost::shared_ptr<tbox::RankTreeStrategy> &rank_tree =
         boost::shared_ptr<tbox::RankTreeStrategy>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~TreeLoadBalancer();

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
    * object won't accidentally interact with other communications.
    *
    * If the duplicate SAMRAI_MPI it is set, the TreeLoadBalancer will
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

   /*!
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
    *
    * Note: This implementation does not yet support non-uniform load
    * balancing.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre !balance_to_anchor || balance_to_anchor->isTransposeOf(balance_to_anchor->getTranspose())
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
      hier::Connector* balance_to_anchor,
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const;

   /*!
    * @brief Print out all members of the class instance to given
    * output stream.
    *
    * @param[in] output_stream
    */
   virtual void
   printClassData(
      std::ostream& output_stream) const;

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


   /*!
    * @brief Enable or disable saving of tree data for diagnostics.
    *
    * @param [in] comm_graph_writer
    * External CommGraphWriter to save tree data to.
    * Use NULL to disable saving.
    */
   void
   setCommGraphWriter(
      const boost::shared_ptr<tbox::CommGraphWriter> &comm_graph_writer )
   {
      d_comm_graph_writer = comm_graph_writer;
   }


   /*!
    * @brief Get the name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:

   typedef double LoadType;


   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int TreeLoadBalancer_LOADTAG0 = 1;
   static const int TreeLoadBalancer_LOADTAG1 = 2;
   static const int TreeLoadBalancer_EDGETAG0 = 3;
   static const int TreeLoadBalancer_EDGETAG1 = 4;
   static const int TreeLoadBalancer_PREBALANCE0 = 5;
   static const int TreeLoadBalancer_PREBALANCE1 = 6;
   static const int TreeLoadBalancer_FIRSTDATALEN = 500;

   static const int TreeLoadBalancer_MIN_NPROC_FOR_AUTOMATIC_MULTICYCLE = 65;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   TreeLoadBalancer(
      const TreeLoadBalancer&);

   void
   operator = (
      const TreeLoadBalancer&);


   /*!
    * @brief Data to save for each sending/receiving process and the
    * subtree at that process.
    *
    * Terminology: The "pruned" parts of the tree are branches that
    * are not open to receiving work from above.  These parts are not
    * counted in the "effective" tree for the purpose of sending down
    * work.
    */
   class SubtreeData {

   public:
      //! @brief Constructor.
      SubtreeData();

      void setPartitioningParams( const PartitioningParams &pparams )
         {
            d_pparams = &pparams;
            d_work_traded.setPartitioningParams(pparams);
         }

      /*!
       * @brief Set the ideal, current and upper limit of the load for
       * the local process.
       */
      void setStartingLoad(
         LoadType ideal,
         LoadType current,
         LoadType upperlimit );

      //! @brief Number of processes in subtree.
      int numProcs() const { return d_num_procs; }
      //! @brief Number of processes in effective subtree.
      int numProcsEffective() const { return d_eff_num_procs; }

      //@{
      //! @name Amount of work in subtree, compared to various references.
      // surplus and deficit are current load compared to ideal.
      LoadType surplus() const { return d_subtree_load_current - d_subtree_load_ideal; }
      LoadType deficit() const { return d_subtree_load_ideal - d_subtree_load_current; }
      LoadType effSurplus() const { return d_eff_load_current - d_eff_load_ideal; }
      LoadType effDeficit() const { return d_eff_load_ideal - d_eff_load_current; }
      // excess and margin are current load compared to upper limit.
      LoadType excess() const { return d_subtree_load_current - d_subtree_load_upperlimit; }
      LoadType margin() const { return d_subtree_load_upperlimit - d_subtree_load_current; }
      LoadType effExcess() const { return d_eff_load_current - d_eff_load_upperlimit; }
      LoadType effMargin() const { return d_eff_load_upperlimit - d_eff_load_current; }
      //@}

      //! @brief Set whether this subtree want work from its parents.
      void setWantsWorkFromParent(bool wants) { d_wants_work_from_parent = wants; }

      //! @brief Get whether this subtree want work from its parents.
      bool getWantsWorkFromParent() const { return d_wants_work_from_parent; }

      //@{
      //! @name Information on work exchanged
      //! @brief Get amount of work exchanged.
      LoadType getExchangeLoad() const { return d_work_traded.getSumLoad(); }
      //! @brief Get count of work exchanged.
      size_t getExchangePackageCount() const { return d_work_traded.size(); }
      //! @brief Get count of originators of the work exchanged.
      size_t getExchangeOriginatorCount() const { return d_work_traded.getNumberOfOriginatingProcesses(); }
      //@}


      //@{
      //! @name Methods supporting load import/export.
      /*!
       * @brief Adjust load to be sent away by taking work from or
       * dumping work into a reserve container.
       */
      LoadType
      adjustOutboundLoad(
         BoxTransitSet& reserve,
         hier::SequentialLocalIdGenerator &id_generator,
         LoadType ideal_load,
         LoadType low_load,
         LoadType high_load );

      //! @brief Move inbound load to the given reserve container.
      void moveInboundLoadToReserve( BoxTransitSet &reserve );

      /*!
       * @brief Incorporate child subtree into this subtree and add
       * child's excess to reserve container.
       */
      void incorporateChild( BoxTransitSet &reserve,
                             const SubtreeData &child );
      //@}


      //@{
      //! @name Packing/unpacking for communication up and down the tree.

      //! @brief Pack load/boxes for sending up to parent.
      void
      packDataToParent(
         tbox::MessageStream &msg) const;

      //! @brief Unpack load/boxes received from child.
      void
      unpackDataFromChild(
         hier::SequentialLocalIdGenerator &id_generator,
         int mpi_rank,
         tbox::MessageStream &msg );

      //! @brief Pack load/boxes for sending down to child.
      void
      packDataToChild(
         tbox::MessageStream &msg) const;

      //! @brief Unpack load/boxes received from parent.
      void
      unpackDataFromParent(
         hier::SequentialLocalIdGenerator &id_generator,
         int mpi_rank,
         tbox::MessageStream &msg );

      //@}

      //! @brief Diagnostic printing.
      void printClassData( const std::string &border, std::ostream &os ) const;

      //! @brief Setup names of timers.
      void setTimerPrefix(const std::string& timer_prefix);

      //! @brief Whether to print steps for debugging.
      void setPrintSteps(bool print_steps) { d_print_steps = print_steps; }

   private:

      /*!
       * @brief Number of processes in subtree
       */
      int d_num_procs;

      /*!
       * @brief Current amount of work in the subtree, including local unassigned
       */
      LoadType d_subtree_load_current;

      /*!
       * @brief Ideal amount of work for the subtree
       */
      LoadType d_subtree_load_ideal;

      /*!
       * @brief Amount of work the subtree is willing to have, based
       * on the load tolerance and upper limits of children.
       */
      LoadType d_subtree_load_upperlimit;

      /*!
       * @brief Number of processes in subtree after pruning independent descendants
       */
      int d_eff_num_procs;

      /*!
       * @brief Current amount of work in the pruned subtree, including local unassigned
       */
      LoadType d_eff_load_current;

      /*!
       * @brief Ideal amount of work for the pruned subtree
       */
      LoadType d_eff_load_ideal;

      /*!
       * @brief Amount of work the pruned subtree is willing to have, based
       * on the load tolerance and upper limit of dependent children.
       */
      LoadType d_eff_load_upperlimit;

      /*!
       * @brief Work to traded (or to be traded).
       *
       * If the object is for the local process, work_traded means
       * traded with the process's *parent*.
       */
      BoxTransitSet d_work_traded;

      /*!
       * @brief Whether subtree expects its parent to send work down.
       */
      bool d_wants_work_from_parent;

      //! @brief Common partitioning parameters.
      const PartitioningParams *d_pparams;

      //@{
      //! @name Debugging and diagnostic data.
      boost::shared_ptr<tbox::Timer> t_pack_load;
      boost::shared_ptr<tbox::Timer> t_unpack_load;
      bool d_print_steps;
      //@}

   }; // SubtreeData declaration.



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

   /*!
    * Move Boxes in balance_box_level from ranks outside of
    * rank_group to ranks inside rank_group.  Modify the given connectors
    * to make them correct following this moving of boxes.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre !balance_to_anchor || (balance_to_anchor->getTranspose().checkTransposeCorrectness(*balance_to_anchor) == 0)
    * @pre !balance_to_anchor || (balance_to_anchor->checkTransposeCorrectness(balance_to_anchor->getTranspose()) == 0)
    */
   void
   prebalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_anchor,
      const tbox::RankGroup& rank_group) const;

   /*!
    * @brief Assign unassigned boxes to local process and set
    * mapping edges where possible without communicating.
    */
   void
   assignUnassignedToLocalProcess(
      hier::BoxLevel& balanced_box_level,
      hier::Connector &balanced_to_unbalanced,
      hier::Connector &unbalanced_to_balanced,
      BoxTransitSet& unassigned ) const;

   void
      removeLocallyOriginatedBoxesFromBoxTransitSet(
      BoxTransitSet& transit_set,
      int local_rank ) const;

   /*!
    * @brief Construct semilocal relationships in
    * unbalanced--->balanced Connector.
    *
    * Constructing semilocal unbalanced--->balanced relationships
    * require communication to determine where exported work ended up.
    * This methods does the necessary communication and constructs
    * these relationship in the given Connector.
    *
    * @param [out] unbalanced_to_balanced Connector to store
    * relationships in.
    *
    * @param [in] kept_imports Work that was imported and locally kept.
    */
   void
   constructSemilocalUnbalancedToBalanced(
      hier::MappingConnector &unbalanced_to_balanced,
      const BoxTransitSet &kept_imports ) const;

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

   /*!
    * @brief Compute the load for a Box.
    */
   double
   computeLoad(
      const hier::Box& box) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.  For non-uniform loads, this method
       * needs the patch data index for the load.  It would summ up
       * the individual cell loads in the cell.
       */
      return double(box.size());
   }

   /*!
    * @brief Compute the load for the Box, restricted to where it
    * intersects a given box.
    */
   double
   computeLoad(
      const hier::Box& box,
      const hier::Box& restriction) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.  For non-uniform loads, this method
       * needs the patch data index for the load.  It would summ up
       * the individual cell loads in the overlap region.
       */
      return double((box * restriction).size());
   }

   /*
    * Count the local workload.
    */
   LoadType
   computeLocalLoads(
      const hier::BoxLevel& box_level) const;

   /*!
    * @brief Given an "unbalanced" BoxLevel, compute the BoxLevel that
    * is load-balanced within the given rank_group and compute the
    * mapping between the unbalanced and balanced BoxLevels.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre d_dim == balance_box_level.getDim()
    */
   void
   loadBalanceWithinRankGroup(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_anchor,
      const tbox::RankGroup& rank_group,
      const double group_sum_load ) const;

   /*!
    * @brief Constrain maximum box sizes in the given BoxLevel and
    * update given Connectors to the changed BoxLevel.
    *
    * @pre !anchor_to_level || anchor_to_level->hasTranspose()
    * @pre d_dim == box_level.getDim()
    */
   void
   constrainMaxBoxSizes(
      hier::BoxLevel& box_level,
      hier::Connector* anchor_to_level) const;

   /*!
    * @brief Compute surplus load per descendent who is still waiting
    * for load from parents.
    */
   LoadType
   computeSurplusPerEffectiveDescendent(
      const LoadType &unassigned_load,
      const LoadType &group_avg_load,
      const std::vector<SubtreeData> &child_subtrees,
      int first_child ) const;

   /*!
    * @brief Create the cycle-based RankGroups the local process
    * belongs in.
    *
    * The RankGroup size increases exponentially with the cycle
    * number such that for the last cycle the rank group includes
    * all processes in d_mpi.
    *
    * @param [out] rank_group
    * @param [out] num_groups
    * @param [out] group_num
    * @param [in] cycle_number
    * @param [in] number_of_cycles
    */
   void
   createBalanceRankGroupBasedOnCycles(
      tbox::RankGroup &rank_group,
      int &num_groups,
      int &group_num,
      const int cycle_number,
      const int number_of_cycles) const;

   /*!
    * @brief Set up the asynchronous communication objects for the
    * given RankGroup.
    *
    * Based on a conceptual process tree with num_children children,
    * set the AsyncCommPeer objects for communication with children
    * and parent.
    *
    * @param [out] child_stage
    * @param [out] child_comms
    * @param [out] parent_stage
    * @param [out] parent_comm
    * @param [in] rank_group
    */
   void
   setupAsyncCommObjects(
      tbox::AsyncCommStage& child_stage,
      tbox::AsyncCommPeer<char> *& child_comms,
      tbox::AsyncCommStage& parent_stage,
      tbox::AsyncCommPeer<char> *& parent_comm,
      const tbox::RankGroup &rank_group ) const;

   /*
    * @brief Undo the set-up done by setupAsyncCommObjects.
    *
    * @pre (d_mpi.getSize() != 1) || ((child_comms == 0) && (parent_comms == 0))
    */
   void
   destroyAsyncCommObjects(
      tbox::AsyncCommPeer<char> *& child_comms,
      tbox::AsyncCommPeer<char> *& parent_comm) const;

   /*!
    * @brief Set up timers for the object.
    */
   void
   setTimers();

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

   //! @brief Max number of processes the a single process may spread its load onto per root cycle.
   int d_max_cycle_spread_ratio;

   //! @brief Whether to allow box breaking.
   bool d_allow_box_breaking;

   //! @brief How to arange a contiguous range of MPI ranks in a tree.
   const boost::shared_ptr<tbox::RankTreeStrategy> d_rank_tree;

   /*!
    * @brief Utility to save data for communication graph output.
    */
   boost::shared_ptr<tbox::CommGraphWriter> d_comm_graph_writer;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   std::vector<int> d_workload_data_id;

   int d_master_workload_data_id;

   /*!
    * @brief Fraction of ideal load a process can accept over and above
    * the ideal it should have.
    *
    * See input parameter "flexible_load_tol".
    */
   double d_flexible_load_tol;

   /*!
    * @brief Load comparison tolerance factor.
    *
    * When low-level methods check whether one candidate is better
    * than the other, ignore improvements less than
    * d_load_comparison_tol*d_global_avg_load.  This prevents infinite
    * loops when the improvement is very near zero.
    */
   double d_load_comparison_tol;

   //@{
   //! @name Data shared with private methods during balancing.
   mutable boost::shared_ptr<PartitioningParams> d_pparams;
   mutable LoadType d_global_avg_load;
   mutable LoadType d_min_load;
   //@}


   /*!
    * @brief Whether to immediately report the results of the load balancing cycles
    * in the log files.
    */
   bool d_report_load_balance;

   /*!
    * @brief See "summarize_map" input parameter.
    */
   char d_summarize_map;

   //@{
   //! @name Used for evaluating peformance.
   bool d_barrier_before;
   bool d_barrier_after;
   //@}

   static const int d_default_data_id;

   /*
    * Performance timers.
    */
   boost::shared_ptr<tbox::Timer> t_load_balance_box_level;
   boost::shared_ptr<tbox::Timer> t_get_map;
   boost::shared_ptr<tbox::Timer> t_use_map;
   boost::shared_ptr<tbox::Timer> t_constrain_size;
   boost::shared_ptr<tbox::Timer> t_map_big_boxes;
   boost::shared_ptr<tbox::Timer> t_load_distribution;
   boost::shared_ptr<tbox::Timer> t_post_load_distribution_barrier;
   boost::shared_ptr<tbox::Timer> t_compute_local_load;
   boost::shared_ptr<tbox::Timer> t_compute_global_load;
   boost::shared_ptr<tbox::Timer> t_compute_tree_load;
   std::vector<boost::shared_ptr<tbox::Timer> > t_compute_tree_load_for_cycle;
   boost::shared_ptr<tbox::Timer> t_send_load_to_children;
   boost::shared_ptr<tbox::Timer> t_send_load_to_parent;
   boost::shared_ptr<tbox::Timer> t_get_load_from_children;
   boost::shared_ptr<tbox::Timer> t_get_load_from_parent;
   boost::shared_ptr<tbox::Timer> t_construct_semilocal;
   boost::shared_ptr<tbox::Timer> t_construct_semilocal_comm_wait;
   boost::shared_ptr<tbox::Timer> t_construct_semilocal_send_edges;
   boost::shared_ptr<tbox::Timer> t_construct_semilocal_local_accounting;
   boost::shared_ptr<tbox::Timer> t_report_loads;
   boost::shared_ptr<tbox::Timer> t_local_balancing;
   boost::shared_ptr<tbox::Timer> t_finish_sends;
   boost::shared_ptr<tbox::Timer> t_pack_edge;
   boost::shared_ptr<tbox::Timer> t_unpack_edge;
   boost::shared_ptr<tbox::Timer> t_children_load_comm;
   boost::shared_ptr<tbox::Timer> t_parent_load_comm;
   boost::shared_ptr<tbox::Timer> t_children_edge_comm;
   boost::shared_ptr<tbox::Timer> t_parent_edge_comm;
   boost::shared_ptr<tbox::Timer> t_barrier_before;
   boost::shared_ptr<tbox::Timer> t_barrier_after;
   boost::shared_ptr<tbox::Timer> t_child_send_wait;
   boost::shared_ptr<tbox::Timer> t_child_recv_wait;
   boost::shared_ptr<tbox::Timer> t_parent_send_wait;
   boost::shared_ptr<tbox::Timer> t_parent_recv_wait;
   boost::shared_ptr<tbox::Timer> t_misc1;
   boost::shared_ptr<tbox::Timer> t_misc2;

   /*
    * Statistics on number of cells and patches generated.
    */
   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;

   // Extra checks independent of optimization/debug.
   char d_print_steps;
   char d_print_pop_steps;
   char d_print_break_steps;
   char d_print_swap_steps;
   char d_print_edge_steps;
   char d_check_connectivity;
   char d_check_map;

};

}
}

#endif
