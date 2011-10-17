/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer
#define included_mesh_TreeLoadBalancer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/Timer.h"

#include <iostream>
#include <vector>
#include <set>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Data to save for each Box that gets passed along the
 * tree edges.
 *
 * The purpose of the BoxInTransit is to associate extra data
 * with a Box as the it is broken up and passed from processor
 * to processor.  A BoxInTransit is a Box going through
 * these changes.  It has a current work load and an orginating
 * Box.  It is passed from process to process and keeps a
 * history of the processes it passed through.  It is assigned a
 * LocalId on every process it passes through, but the LocalId history
 * is not kept.
 */
struct BoxInTransit {
   typedef hier::Box Box;
   typedef hier::LocalId LocalId;

   /*!
    * @brief Constructor
    *
    * @param[in] dim
    */
   BoxInTransit(
      const tbox::Dimension& dim):
      box(dim),
      orig_box(dim)
   {
   }

   /*!
    * @brief Construct a newly birthed BoxInTransit.
    *
    * @param[in] other
    */
   BoxInTransit(
      const hier::Box& other):
      box(other),
      orig_box(other),
      load(other.size()),
      proc_hist()
   {
   }

   /*!
    * @brief Construct new object using integer data packed by
    * putToIntBuffer().
    *
    * @param[i/o] ptr Pointer to integer data in buffer.  On return,
    * @c ptr will be advanced past the data used by this object.
    *
    * @param[i] dim
    */
   BoxInTransit(
      const int *&ptr,
      const tbox::Dimension &dim ):
      box(dim),
      orig_box(dim),
      load(0),
      proc_hist(0)
   {
      getFromIntBuffer(ptr);
      ptr += commBufferSize();
   }

   /*!
    * @brief Construct new object based on an existing object, taking
    * the exiting object's origin data and process history, but
    * appending a new process in the history and using a new box.
    *
    * @param[in] other
    *
    * @param[in] box
    *
    * @param[in] rank
    *
    * @param[in] local_id
    */
   BoxInTransit(
      const BoxInTransit& other,
      const hier::Box& box,
      int rank,
      LocalId local_id):
      box(box, local_id, rank, other.orig_box.getBlockId()),
      orig_box(other.orig_box),
      load(box.size()),
      proc_hist(other.proc_hist)
   {
      if (rank != other.getOwnerRank()) {
         proc_hist.push_back(other.getOwnerRank());
      }
   }

   /*!
    * @brief Assignment operator
    *
    * @param[in] other
    */
   const BoxInTransit& operator = (
      const BoxInTransit& other) {
      box = other.box;
      orig_box = other.orig_box;
      load = other.load;
      proc_hist = other.proc_hist;
      return *this;
   }

   //! @brief The Box.
   hier::Box box;

   //! @brief Originating Box.
   hier::Box orig_box;

   //! @brief Normalized work load.
   int load;

   /*!
    * @brief History of processors passed in transit.  Each time a
    * processor passes on a Box, it should append its rank to
    * the proc_hist.
    */
   std::vector<int> proc_hist;

   //! @brief Return the owner rank.
   int getOwnerRank() const {
      return box.getOwnerRank();
   }

   //! @brief Return the LocalId.
   LocalId getLocalId() const {
      return box.getLocalId();
   }

   //! @brief Return the Box.
   hier::Box& getBox() {
      return box;
   }

   //! @brief Return the Box.
   const hier::Box& getBox() const {
      return box;
   }

   /*!
    * @brief Return number of ints required for putting a putting the
    * object in message passing buffer.
    */
   int commBufferSize() const {
      const tbox::Dimension& dim(box.getDim());
      return 2 * hier::Box::commBufferSize(dim) + 2 + static_cast<int>(proc_hist.size());
   }

   /*!
    * @brief Put self into a int buffer.
    *
    * This is the opposite of getFromIntBuffer().  Number of ints
    * written is given by commBufferSize().
    *
    * If skip_last_owner is true, and the last owner in proc_hist will
    * be skipped (proc_hist must be non-empty) and the number of integers
    * put in the buffer would be one less than commBufferSize().
    */
   void putToIntBuffer(
      int* buffer,
      bool skip_last_owner=false) const {
      const tbox::Dimension& dim(box.getDim());
      box.putToIntBuffer(buffer);
      buffer += hier::Box::commBufferSize(dim);
      orig_box.putToIntBuffer(buffer);
      buffer += hier::Box::commBufferSize(dim);
      *(buffer++) = load;
      if (skip_last_owner) { TBOX_ASSERT( !proc_hist.empty() ); }
      *(buffer++) = static_cast<int>(proc_hist.size()-skip_last_owner);
      for (unsigned int i = 0; i < proc_hist.size()-skip_last_owner; ++i) {
         buffer[i] = proc_hist[i];
      }
   }

   /*!
    * @brief Set attributes according to data in int buffer.
    *
    * This is the opposite of putToIntBuffer().  Number of ints read
    * is given by what commBufferSize() AFTER this method is called.
    */
   void getFromIntBuffer(
      const int* buffer) {
      const tbox::Dimension& dim(box.getDim());
      box.getFromIntBuffer(buffer);
      buffer += hier::Box::commBufferSize(dim);
      orig_box.getFromIntBuffer(buffer);
      buffer += hier::Box::commBufferSize(dim);
      load = *(buffer++);
      proc_hist.clear();
      proc_hist.insert(proc_hist.end(), *(buffer++), 0);
      for (unsigned int i = 0; i < proc_hist.size(); ++i) {
         proc_hist[i] = buffer[i];
      }
   }

   /*!
    * @brief Stuff into an outgoing message destined for the previous
    * owner.
    *
    * prev_owner must not be empty.  This method pops the next value
    * from prev_owner (changing the object's state!) and packs the
    * object into the message for that owner.
    *
    * @param outgoing_messages Map of outgoing messages, indexed by
    * recipient rank.  On return, the message corresponding to the
    * previous owner will be grown by commBufferSize()-1.
    */
   void packForPreviousOwner(
      std::map<int,std::vector<int> > &outgoing_messages ) const {
      const int prev_owner = proc_hist.back();
      std::vector<int> &msg(outgoing_messages[prev_owner]);
      const int cbs = commBufferSize();
      msg.insert(msg.end(), cbs-1, 0);
      putToIntBuffer(&msg[msg.size() - (cbs-1)], true);
   }

   //! @brief Stream-insert operator.
   friend std::ostream&
   operator << (
      std::ostream& co,
      const BoxInTransit& r);
};

/*!
 * @brief Comparison functor for sorting BoxInTransit
 * from bigger to smaller boxes.
 */
struct BoxInTransitMoreLoad {
   bool operator () (
      const BoxInTransit& a,
      const BoxInTransit& b) const {
      if (a.getBox().size() != b.getBox().size()) {
         return a.load > b.load;
      }
      return a.box.getId() < b.box.getId();
   }
};

/*!
 * @brief Provides load balancing routines for AMR hierarchy by
 * implemementing the LoadBalancerStrategy.
 *
 * Currently, only uniform load balancing is supported.  Eventually,
 * non-uniform load balancing should be supported.  (Non-uniform load
 * balancing is supported by the CutAndPackLoadBalancer class.)
 *
 * Inputs and their default values:
 *
 * No special inputs are required for this class.
 *
 * @verbatim
 * report_load_balance = TRUE // Write out load balance report in log
 * n_root_cycles = -1         // Number of root cycles to use for
 *                            // reaching final partitioning. Nominally 1.
 *                            // Can be set higher to reduce negative
 *                            // performance effects of extremely poor initial
 *                            // load balance.  Set to -1 for "automatic".
 *                            // Set to zero to effectively bypass load balancing.
 * balance_penalty_wt = 1.0   // Relative weight for computing combined box breaking
 *                            // penalty:  How much to penalize imbalance.
 * surface_penalty_wt = 1.0   // Relative weight for computing combined box breaking
 *                            // penalty:  How much to penalize new surfaces.
 * slender_penalty_wt = 1.0   // Relative weight for computing combined box breaking
 *                            // penalty:  How much to penalize slender boxes.
 * @endverbatim
 *
 * @see mesh::LoadBalanceStrategy
 */

class TreeLoadBalancer:
   public LoadBalanceStrategy
{
public:
   /*!
    * @brief Initializing constructor sets object state to default or,
    * if provided, to parameters in database.
    *
    * @param[in] dim
    *
    * @param[in] name User-defined std::string identifier used for error
    * reporting and timer names.  If omitted, "TreeLoadBalancer"
    * is used.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    */
   TreeLoadBalancer(
      const tbox::Dimension& dim,
      const std::string& name = std::string("TreeLoadBalancer"),
      tbox::Pointer<tbox::Database> input_db =
         tbox::Pointer<tbox::Database>(NULL));

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~TreeLoadBalancer();

   /*!
    * @brief Set the internal communicator to a duplicate of the given
    * communicator.
    *
    * The given communicator must be a valid communicator.
    *
    * The given communicator is duplicated for private use.  This
    * requires a global communication, so all processes in the
    * communicator must call it.  The advantage of a duplicate
    * communicator is that it ensures the communications for the
    * object won't accidentally interact with other communications.
    *
    * If the duplicate MPI communicator it is set, the
    * TreeLoadBalancer will only balance BoxLevels with
    * congruent SAMRAI_MPI objects and will use the duplicate
    * communicator for communications.  Otherwise, the communicator of
    * the BoxLevel will be used.  The duplicate MPI communicator
    * is freed when the object is destructed, or freeMPICommunicator()
    * is called.
    *
    * TODO: For uniformity, we should pass in a SAMRAI_MPI instead of
    * a communicator.  This class is the only one that still uses the
    * communicator instead of the SAMRAI_MPI.
    */
   void
   setSAMRAI_MPI(
      const tbox::SAMRAI_MPI& samrai_mpi);

   /*!
    * @brief Free the internal MPI communicator, if any has been set.
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
    */
   void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1);

   /*!
    * @brief Configure the load balancer to load balance boxes by
    * assuming all cells on the specified level or all hierarchy
    * levels are weighted equally.
    *
    * @param level_number
    * Optional integer number for level on which uniform
    * workload estimate will be used.  If the level
    * number is not specified, a uniform workload
    * estimate will be used on all levels.
    */
   void
   setUniformWorkload(
      int level_number = -1);

   /*!
    * @brief Return true if load balancing procedure for given level
    * depends on patch data on mesh; otherwise return false.
    *
    * This can be used to determine whether a level needs to be
    * rebalanced although its box configuration is unchanged.  This
    * function is pure virtual in the LoadBalanceStrategy base class.
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
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_mapped_box_level,
      hier::Connector& balance_to_anchor,
      hier::Connector& anchor_to_balance,
      const tbox::Pointer<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const hier::Connector& unbalanced_to_attractor,
      const hier::Connector& attractor_to_unbalanced,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_mapped_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& = tbox::RankGroup()) const;

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
      std::ostream& output_stream = tbox::plog) const;

   /*!
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const;

private:
   /*
    * Static integer constants.
    */
   static const int TreeLoadBalancer_LOADTAG0;
   static const int TreeLoadBalancer_LOADTAG1;
   static const int TreeLoadBalancer_EDGETAG0;
   static const int TreeLoadBalancer_EDGETAG1;
   static const int TreeLoadBalancer_PREBALANCE0;
   static const int TreeLoadBalancer_PREBALANCE1;
   static const int TreeLoadBalancer_FIRSTDATALEN;

   typedef hier::Box Box;

   typedef hier::LocalId LocalId;

   typedef hier::BoxLevel BoxLevel;

   typedef mesh::BoxInTransitMoreLoad BoxInTransitMoreLoad;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   TreeLoadBalancer(
      const TreeLoadBalancer&);

   void
   operator = (
      const TreeLoadBalancer&);

   /*!
    * @brief A set of BoxInTransit, sorted from highest load to lowest load.
    */
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad> TransitSet;

   /*!
    * @brief Data to save for each subtree.
    */
   struct SubtreeLoadData {
      SubtreeLoadData():num_procs(0),
         total_work(0),
         load_exported(0),
         load_imported(0) {
      }
      /*!
       * @brief Number of nodes in child's subtree
       * (unused for parent's data).
       */
      int num_procs;
      /*!
       * @brief Current total work amount in child's subtree
       */
      int total_work;
      /*!
       * @brief Work exported (or scheduled for export) to nonlocal process.
       *
       * Value is often used to determine if more communication is needed.
       *
       * For local data, this refers to the work exported to parent.
       */
      int load_exported;
      /*!
       * @brief Work imported from nonlocal process.
       *
       * Value is often used to determine if more communication is needed.
       *
       * For local data, this refers to the work imported from parent.
       */
      int load_imported;
      /*!
       * @brief Ideal work amount in subtree
       */
      int ideal_work;
      /*!
       * @brief Nodes to export.
       *
       * For local data, this refers to exporting to parent.
       */
      TransitSet for_export;
   };

   void
   assertNoMessageForPrivateCommunicator() const;

   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      tbox::Pointer<tbox::Database> db);

   void
   sortIntVector(
      hier::IntVector& ordered_dims,
      const hier::IntVector& vector) const;

   /*!
    * Move Boxes in balance_mapped_box_level from ranks outside of
    * rank_group to ranks inside rank_group.  Modify the given connectors
    * to make them correct following this moving of boxes.
    */
   void
   prebalanceBoxLevel(
      hier::BoxLevel& balance_mapped_box_level,
      hier::Connector& balance_to_anchor,
      hier::Connector& anchor_to_balance,
      const tbox::RankGroup& rank_group) const;

   /*!
    * @brief Reassign loads from one TransitSet to another.
    *
    * @param ideal_transfer Amount of load to reassign from src to
    * dst.  If negative, reassign the load from dst to src.
    *
    * @param actual_transfer Amount of load transfered.  If positive,
    * transfer load from src to dst (if negative, from dst to src).
    *
    * @param next_available_index Index for guaranteeing new
    * Boxes are uniquely numbered.
    */
   void
   reassignLoads(
      TransitSet& src,
      TransitSet& dst,
      int& actual_transfer,
      hier::LocalId& next_available_index,
      const int ideal_transfer ) const;

   /*
    * @brief Shift load from src to dst by swapping BoxInTransit
    * between them.
    *
    * @param ideal_transfer Amount of load to reassign from src to
    * dst.  If negative, reassign the load from dst to src.
    *
    * @param actual_transfer Amount of load transfered.  If positive,
    * transfer load from src to dst (if negative, from dst to src).
    *
    * @param actual_transfer Amount of load transfered.  If positive,
    * transfer load from src to dst (if negative, from dst to src).
    */
   bool
   shiftLoadsBySwapping(
      TransitSet& src,
      TransitSet& dst,
      int& actual_transfer,
      const int ideal_transfer ) const;

   /*
    * @brief Shift load from src to dst by various box breaking strategies.
    * choosing the break that gives the best overall penalty.
    *
    * @param ideal_transfer Amount of load to reassign from src to
    * dst.  If negative, reassign the load from dst to src.
    *
    * @param actual_transfer Amount of load transfered.  If positive,
    * transfer load from src to dst (if negative, from dst to src).
    *
    * @param next_available_index Index for guaranteeing new
    * Boxes are uniquely numbered.
    */
   bool
   shiftLoadsByBreaking(
      TransitSet& src,
      TransitSet& dst,
      int& actual_transfer,
      hier::LocalId& next_available_index,
      const int ideal_transfer ) const;

   bool
   findLoadSwapPair(
      TransitSet& src,
      TransitSet& dst,
      int& actual_transfer,
      TransitSet::iterator& isrc,
      TransitSet::iterator& idst,
      const int ideal_transfer ) const;

   /*!
    * @brief Pack load/boxes for sending.
    */
   void
   packSubtreeLoadData(
      std::vector<int>& msg,
      const SubtreeLoadData& load_data) const;

   /*!
    * @brief Unpack load/boxes received.
    */
   void
   unpackSubtreeLoadData(
      SubtreeLoadData& proc_data,
      TransitSet& receiving_bin,
      hier::LocalId& next_available_index,
      const int* received_data,
      int received_data_length ) const;

   void
   unpackAndRouteNeighborhoodSets(
      std::map<int,std::vector<int> > &outgoing_messages,
      hier::Connector& unbalanced_to_balanced,
      const int* received_data,
      int received_data_length ) const;

   /*!
    * @brief Break off a given load size from a given Box.
    *
    * @param mapped_box Box to break.
    * @param @ideal_load_to_break Ideal load to break.
    * This is not guaranteed to be met.  The actual
    * amount broken off may be slightly over or under.
    * @param breakoff Boxes broken off (usually just one).
    * @parem leftover Remainder of Box after breakoff is gone.
    */
   bool
   breakOffLoad(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& mapped_box,
      double ideal_load_to_break ) const;

   /*!
    * @brief Computes surface area of a list of boxes.
    */
   double
   computeBoxSurfaceArea(
      const std::vector<hier::Box>& boxes) const;

   /*!
    * @brief Computes the surface area of a box.
    */
   int
   computeBoxSurfaceArea(
      const hier::Box& box) const;

   double
   combinedBreakingPenalty(
      double balance_penalty,
      double surface_penalty,
      double slender_penalty) const;

   double
   computeBalancePenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b,
      double off_balance) const;
   double
   computeBalancePenalty(
      const TransitSet& a,
      const TransitSet& b,
      double off_balance) const;
   double
   computeBalancePenalty(
      const hier::Box& a,
      double off_balance) const;

   double
   computeSurfacePenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b) const;
   double
   computeSurfacePenalty(
      const TransitSet& a,
      const TransitSet& b) const;
   double
   computeSurfacePenalty(
      const hier::Box& a) const;

   double
   computeSlenderPenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b) const;
   double
   computeSlenderPenalty(
      const TransitSet& a,
      const TransitSet& b) const;
   double
   computeSlenderPenalty(
      const hier::Box& a) const;

   bool
   breakOffLoad_planar(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& mapped_box,
      double ideal_load_to_break,
      const tbox::Array<tbox::Array<bool> >& bad_cuts ) const;

   bool
   breakOffLoad_cubic1(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& mapped_box,
      double ideal_load_to_give,
      const tbox::Array<tbox::Array<bool> >& bad_cuts ) const;

   void
   burstBox(
      std::vector<hier::Box>& boxes,
      const hier::Box& bursty,
      const hier::Box& solid ) const;

   /*
    * Utility functions to determine parameter values for level.
    */
   int
   getWorkloadDataId(
      int level_number) const;

   /*!
    * @brief Compute the load for a Box.
    */
   double
   computeLoad(
      const hier::Box& mapped_box) const;

   /*!
    * @brief Compute the load for the Box where it intersects the given box.
    */
   double
   computeLoad(
      const hier::Box& mapped_box,
      const hier::Box& box) const;

   /*
    * Count the local workload.
    */
   double
   computeLocalLoads(
      const hier::BoxLevel& mapped_box_level) const;

   /*!
    * @brief Given an "unbalanced" BoxLevel, compute the BoxLevel that
    * is load-balanced within the given rank_group and compute the
    * mapping between the unbalanced and balanced BoxLevels.
    */
   void
   computeLoadBalancingMapWithinRankGroup(
      hier::BoxLevel& balanced_mapped_box_level,
      hier::Connector& unbalanced_to_balanced,
      hier::Connector& balanced_to_unbalanced,
      const hier::BoxLevel& unbalanced_mapped_box_level,
      const tbox::RankGroup& rank_group,
      const int cycle_number,
      const int number_of_cycles,
      const double local_load,
      const double global_sum_load ) const;

   /*!
    * @brief Compute BoxLevel conforming to max size constraint and
    * the mapping to that BoxLevel.
    *
    * The mapping is entirely local (no transfering of work).
    */
   void
   mapOversizedBoxes(
      hier::BoxLevel& constrained,
      hier::Connector& unconstrained_to_constrained,
      const hier::BoxLevel& unconstrained ) const;

   /*!
    * @brief Create the cycle-based RankGroups the local process
    * belongs in.
    *
    * The RankGroup size increases exponentially with the cycle
    * number such that for the last cycle the rank group includes
    * all processes in d_mpi.
    *
    * @param [o] rank_group
    * @param [o] num_groups
    * @param [o] group_num
    * @param [i] cycle_number
    * @param [i] number_of_cycles
    */
   void createBalanceRankGroupBasedOnCycles(
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
    * @param [o] num_children
    * @param [o] child_comms
    * @param [o] parent_send
    * @param [o] parent_recv
    * @param [i/o] parent_recv
    * @param [i] rank_group
    */
   void setupAsyncCommObjects(
      int& num_children,
      tbox::AsyncCommPeer<int> *& child_comms,
      tbox::AsyncCommPeer<int> *& parent_send,
      tbox::AsyncCommPeer<int> *& parent_recv,
      tbox::AsyncCommStage& comm_stage,
      const tbox::RankGroup &rank_group ) const;

   /*
    * @brief Undo the set-up done by setupAsyncCommObjects.
    */
   void
   destroyAsyncCommObjects(
      tbox::AsyncCommPeer<int> *& child_comms,
      tbox::AsyncCommPeer<int> *& parent_send,
      tbox::AsyncCommPeer<int> *& parent_recv) const;

   void
   setShadowData(
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_mapped_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const hier::IntVector& refinement_ratio) const;
   void
   unsetShadowData() const;

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
   tbox::SAMRAI_MPI d_mpi_dup;

   int d_n_root_cycles;

   const int d_degree;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   tbox::Array<int> d_workload_data_id;

   int d_master_workload_data_id;

   /*!
    * @brief Weighting factor for penalizing imbalance.
    *
    * @see combinedBreakingPenalty().
    */
   double d_balance_penalty_wt;

   /*!
    * @brief Weighting factor for penalizing new suraces.
    *
    * @see combinedBreakingPenalty().
    */
   double d_surface_penalty_wt;

   /*!
    * @brief Weighting factor for penalizing slenderness.
    *
    * @see combinedBreakingPenalty().
    */
   double d_slender_penalty_wt;

   /*!
    * @brief How high a slenderness ratio we can tolerate before penalizing.
    */
   double d_slender_penalty_threshold;

   /*!
    * @brief Extra penalty weighting applied before cutting.
    *
    * Set to range [1,ininity).
    * Higher value forces more agressive cutting but can produce more slivers.
    */
   double d_precut_penalty_wt;

   //@{
   //! @name Data shared with private methods during balancing.
   mutable tbox::SAMRAI_MPI d_mpi;
   mutable hier::IntVector d_min_size;
   mutable hier::IntVector d_max_size;
   mutable hier::BoxList d_domain_boxes;
   mutable hier::IntVector d_bad_interval;
   mutable hier::IntVector d_cut_factor;
   mutable double d_global_avg_load;
   //@}

   mutable tbox::Array<int> d_output_procs;
   bool d_using_all_procs;

   /*!
    * @brief Whether to immediately report the results of the load balancing cycles
    * in the log files.
    */
   bool d_report_load_balance;

   //@{
   //! @name Used for evaluating peformance.
   bool d_barrier_before;
   bool d_barrier_after;
   //@}

   static const int d_default_data_id;

   /*
    * Performance timers.
    */
   tbox::Pointer<tbox::Timer> t_load_balance_mapped_box_level;
   tbox::Pointer<tbox::Timer> t_get_map;
   tbox::Pointer<tbox::Timer> t_use_map;
   tbox::Pointer<tbox::Timer> t_constrain_size;
   tbox::Pointer<tbox::Timer> t_map_big_boxes;
   tbox::Pointer<tbox::Timer> t_compute_local_load;
   tbox::Pointer<tbox::Timer> t_compute_global_load;
   tbox::Pointer<tbox::Timer> t_compute_tree_load;
   tbox::Pointer<tbox::Timer> t_compute_tree_load0;
   tbox::Pointer<tbox::Timer> t_compute_tree_load1;
   tbox::Pointer<tbox::Timer> t_compute_tree_load2;
   tbox::Pointer<tbox::Timer> t_reassign_loads;
   tbox::Pointer<tbox::Timer> t_shift_loads_by_swapping;
   tbox::Pointer<tbox::Timer> t_shift_loads_by_breaking;
   tbox::Pointer<tbox::Timer> t_find_swap_pair;
   tbox::Pointer<tbox::Timer> t_break_off_load;
   tbox::Pointer<tbox::Timer> t_find_bad_cuts;
   tbox::Pointer<tbox::Timer> t_send_load_to_children;
   tbox::Pointer<tbox::Timer> t_send_load_to_parent;
   tbox::Pointer<tbox::Timer> t_get_load_from_children;
   tbox::Pointer<tbox::Timer> t_get_load_from_parent;
   tbox::Pointer<tbox::Timer> t_send_edge_to_children;
   tbox::Pointer<tbox::Timer> t_send_edge_to_parent;
   tbox::Pointer<tbox::Timer> t_get_edge_from_children;
   tbox::Pointer<tbox::Timer> t_get_edge_from_parent;
   tbox::Pointer<tbox::Timer> t_report_loads;
   tbox::Pointer<tbox::Timer> t_local_balancing;
   tbox::Pointer<tbox::Timer> t_local_edges;
   tbox::Pointer<tbox::Timer> t_finish_comms;
   tbox::Pointer<tbox::Timer> t_misc1;
   tbox::Pointer<tbox::Timer> t_misc2;
   tbox::Pointer<tbox::Timer> t_misc3;
   tbox::Pointer<tbox::Timer> t_pack_load;
   tbox::Pointer<tbox::Timer> t_unpack_load;
   tbox::Pointer<tbox::Timer> t_unpack_edge;
   tbox::Pointer<tbox::Timer> t_children_load_comm;
   tbox::Pointer<tbox::Timer> t_parent_load_comm;
   tbox::Pointer<tbox::Timer> t_children_edge_comm;
   tbox::Pointer<tbox::Timer> t_parent_edge_comm;
   tbox::Pointer<tbox::Timer> t_barrier_before;
   tbox::Pointer<tbox::Timer> t_barrier_after;
   tbox::Pointer<tbox::Timer> t_MPI_wait;

   /*
    * Statistics on number of cells and patches generated.
    */
   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;

   // Extra checks independent of optimization/debug.
   char d_print_steps;
   char d_print_break_steps;
   char d_print_swap_steps;
   char d_print_edge_steps;
   char d_check_connectivity;
   char d_check_map;

};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/mesh/TreeLoadBalancer.I"
#endif
#endif
