/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Asynchronous Berger-Rigoutsos algorithm wrapper
 *
 ************************************************************************/
#ifndef included_mesh_BergerRigoutsos_C
#define included_mesh_BergerRigoutsos_C

#include <stdlib.h>

#include "SAMRAI/mesh/BergerRigoutsos.h"

#include "SAMRAI/mesh/BergerRigoutsosNode.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

namespace SAMRAI {
namespace mesh {


/*
 ************************************************************************
 * Constructor stores parameters of options for ussing
 * the asynchronous Berger-Rigoutsos implementation.
 ************************************************************************
 */
BergerRigoutsos::BergerRigoutsos(
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),

   d_object_timers(0),
   d_relaunch_queue(),
   d_comm_stage(),
   d_algo_advance_mode(ADVANCE_SOME),
   d_new_box_level(),
   d_tag_to_new(),
   d_root_boxes(),
   // Parameters not from clustering algorithm interface ...
   d_max_inflection_cut_from_center(1.0),
   d_inflection_cut_threshold_ar(0.0),
   d_max_box_size(hier::IntVector(dim, tbox::MathUtilities<int>::getMax())),
   d_min_box_size_from_cutting(dim, 0),
   // Parameters from clustering algorithm interface ...
   d_tag_data_index(-1),
   d_tag_val(1),
   d_min_box(dim),
   d_efficiency_tol(0.80),
   d_combine_tol(0.80),
   // Implementation flags and data...
   d_compute_relationships(2),
   d_relationship_senders(),
   d_relationship_messages(),
   d_max_gcw(dim, 1),
   d_owner_mode(MOST_OVERLAP),
   // Communication parameters ...
   d_mpi_object(MPI_COMM_NULL),
   d_tag_upper_bound(-1),
   d_available_mpi_tag(-1),
   // Analysis support ...
   d_log_node_history(false),
   d_num_tags_in_all_nodes(0),
   d_max_tags_owned(0),
   d_num_nodes_allocated(0),
   d_max_nodes_allocated(0),
   d_num_nodes_active(0),
   d_max_nodes_active(0),
   d_num_nodes_owned(0),
   d_max_nodes_owned(0),
   d_num_nodes_commwait(0),
   d_max_nodes_commwait(0),
   d_num_nodes_completed(0),
   d_max_generation(0),
   d_num_boxes_generated(0),
   d_num_conts_to_complete(0),
   d_max_conts_to_complete(0),
   d_num_nodes_existing(0),

   // d_max_box_size(hier::IntVector(d_dim, tbox::MathUtilities<int>::getMax())),
   // d_max_inflection_cut_from_center(1.0),
   // d_inflection_cut_threshold_ar(0.0),
   // d_log_node_history(false),
   d_log_cluster_summary(false),
   d_log_cluster(false),
   // d_owner_mode("MOST_OVERLAP"),
   // d_algo_advance_mode("ADVANCE_SOME"),
   d_sort_output_nodes(false),
   d_check_min_box_size('w'),
   // d_min_box_size_from_cutting(d_dim, 0),
   d_barrier_before(false),
   d_barrier_after(false)
{

   /*
    * Set database-dependent parameters or cache them for use
    * when we construct a dendogram root.
    */
   getFromInput(input_db);

   setObjectTimers(s_default_timer_prefix);
   // Set the timer for the communication stage's MPI waiting.
   d_comm_stage.setCommunicationWaitTimer(d_object_timers->t_MPI_wait);
}

BergerRigoutsos::~BergerRigoutsos()
{
   if (d_mpi_object.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      // Free the private communicator (if SAMRAI_MPI has not been finalized).
      int flag;
      tbox::SAMRAI_MPI::Finalized(&flag);
      if (!flag) {
         d_mpi_object.freeCommunicator();
      }
   }
}

void
BergerRigoutsos::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{
   if (input_db) {

      if (input_db->isInteger("max_box_size")) {
         input_db->getIntegerArray("max_box_size",
            &d_max_box_size[0],
            d_dim.getValue());
         for (int i = 0; i < d_dim.getValue(); ++i) {
            if (!(d_max_box_size[i] > 0)) {
               INPUT_RANGE_ERROR("max_box_size");
            }
         }
      }

      d_max_inflection_cut_from_center =
         input_db->getDoubleWithDefault("DEV_max_inflection_cut_from_center",
                                        d_max_inflection_cut_from_center);
      d_inflection_cut_threshold_ar =
         input_db->getDoubleWithDefault("DEV_inflection_cut_threshold_ar",
                                        d_inflection_cut_threshold_ar);
      if (input_db->isInteger("DEV_min_box_size_from_cutting")) {
         input_db->getIntegerArray("DEV_min_box_size_from_cutting",
            &d_min_box_size_from_cutting[0],
            d_dim.getValue());
      }
      if (input_db->isInteger("DEV_min_box_size_from_cutting")) {
         input_db->getIntegerArray("DEV_min_box_size_from_cutting",
            &d_min_box_size_from_cutting[0],
            d_dim.getValue());
      }
      d_log_node_history =
         input_db->getBoolWithDefault("DEV_log_node_history", false);
      d_log_cluster_summary =
         input_db->getBoolWithDefault("DEV_log_cluster_summary", false);
      d_log_cluster =
         input_db->getBoolWithDefault("DEV_log_cluster", false);

      std::string algo_advance_mode =
         input_db->getStringWithDefault("DEV_algo_advance_mode", "ADVANCE_SOME");
      if (!(algo_advance_mode == "ADVANCE_SOME" ||
            algo_advance_mode == "ADVANCE_ANY" ||
            algo_advance_mode == "SYNCHRONOUS")) {
         INPUT_VALUE_ERROR("DEV_algo_advance_mode");
      }
      setAlgorithmAdvanceMode(algo_advance_mode);

      std::string owner_mode =
         input_db->getStringWithDefault("DEV_owner_mode", "MOST_OVERLAP");
      if (!(owner_mode == "SINGLE_OWNER" ||
            owner_mode == "MOST_OVERLAP" ||
            owner_mode == "FEWEST_OWNED" ||
            owner_mode == "LEAST_ACTIVE")) {
         INPUT_VALUE_ERROR("DEV_owner_mode");
      }
      setOwnerMode(owner_mode);

      d_sort_output_nodes =
         input_db->getBoolWithDefault("sort_output_nodes", false);

      std::string tmp_str;

      tmp_str =
         input_db->getStringWithDefault("check_min_box_size", std::string("WARN"));
      if (!(tmp_str == "IGNORE" || tmp_str == "WARN" || tmp_str == "ERROR")) {
         INPUT_VALUE_ERROR("check_min_box_size");
      }
      d_check_min_box_size = char(tolower(*tmp_str.c_str()));

      d_barrier_before =
         input_db->getBoolWithDefault("DEV_barrier_before", false);
      d_barrier_after =
         input_db->getBoolWithDefault("DEV_barrier_after", false);
   }
}


/*
 ************************************************************************
 *
 * Implement the BoxGeneratorStrategy interface method using
 * the asynchronous Berger-Rigoutsos implementation.
 *
 * Create objects for using the ABR recursion tree, set options for
 * using the ABR implementation, then run it.
 *
 * The output boxes from the dendogram root is in the form of a
 * BoxLevel.  This method postprocess that data to
 * convert the output to the box list form required by the
 * box clustering strategy interface.
 *
 ************************************************************************
 */
void
BergerRigoutsos::findBoxesContainingTags(
   boost::shared_ptr<hier::BoxLevel>& new_box_level,
   boost::shared_ptr<hier::Connector>& tag_to_new,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_data_index,
   const int tag_val,
   const hier::BoxContainer& bound_boxes,
   const hier::IntVector& min_box,
   const double efficiency_tol,
   const double combine_tol,
   const hier::IntVector& max_gcw)
{
   TBOX_ASSERT(!bound_boxes.isEmpty());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*tag_level,
      *(bound_boxes.begin()),
      min_box,
      max_gcw);

   tbox::SAMRAI_MPI mpi(tag_level->getBoxLevel()->getMPI());

   for (hier::BoxContainer::const_iterator bb_itr = bound_boxes.begin();
        bb_itr != bound_boxes.end(); ++bb_itr) {
      if (!(bb_itr->numberCells() >= min_box)) {
         if (d_check_min_box_size == 'e') {
            TBOX_ERROR("BergerRigoutsos::findBoxesContainingTags input error:\n"
               << "Input box " << *bb_itr << " has size " << bb_itr->numberCells()
               << "\nwhich is already smaller than the minimum box size\n"
               << min_box << "\n\n"
               << "To ignore or just issue a warning, see the input parameter\n"
               << "check_min_box_size.\n");
         } else if (d_check_min_box_size == 'w') {
            TBOX_WARNING("BergerRigoutsos::findBoxesContainingTags input warning:\n"
               << "Input box " << *bb_itr << " has size " << bb_itr->numberCells()
               << "\nwhich is already smaller than the minimum box size\n"
               << min_box << "\n\n"
               << "To ignore or issue error, see the input parameter\n"
               << "check_min_box_size.\n");
         }
         if (bb_itr->empty()) {
            TBOX_ERROR("BergerRigoutsos: empty bounding box not allowed.");
         }
      }
   }


   if (d_barrier_before) {
      d_object_timers->t_barrier_before->start();
      mpi.Barrier();
      d_object_timers->t_barrier_before->stop();
   }

   for ( hier::BoxContainer::const_iterator bi=bound_boxes.begin();
         bi!=bound_boxes.end(); ++bi ) {
      if (bi->empty()) {
         TBOX_ERROR("BergerRigoutsos: empty bounding box not allowed.");
      }
   }

   d_object_timers->t_find_boxes_containing_tags->start();


   /*
    * Set up some internal data then call
    * clusterAndComputeRelationships run the clustering algorithm.
    */

   resetCounters();


   /*
    * Set parameters received from findBoxesContainingTags() virtual
    * interface.
    */

   d_tag_data_index = tag_data_index;
   d_tag_val = tag_val;
   d_min_box = min_box;
   d_efficiency_tol = efficiency_tol;
   d_combine_tol = combine_tol;

   setComputeRelationships("BIDIRECTIONAL", max_gcw);

   d_tag_level = tag_level;
   d_root_boxes = bound_boxes;

   /*
    * If d_mpi_object has not been set, then user wants to do use the
    * MPI in tag_level (nothing special).  If it has been set, it is a
    * duplicate MPI, so don't change it.
    */
   if ( d_mpi_object.getCommunicator() == MPI_COMM_NULL ) {
      d_mpi_object = d_tag_level->getBoxLevel()->getMPI();
      setupMPIDependentData();
   }
#if defined(DEBUG_CHECK_ASSERTIONS)
   else {
      if ( !checkMPICongruency() ) {
         TBOX_ERROR("BergerRigoutsosNode::clusterAndComputeRelationships:\n"
                    << "The communicator of the input tag BoxLevel ("
                    << d_tag_level->getBoxLevel()->getMPI().getCommunicator()
                    << " is not congruent with the MPI communicator ("
                    << d_mpi_object.getCommunicator()
                    << " duplicated in the call to useDuplicateMPI().\n"
                    << "If you call useDuplicateMPI(), you are restricted\n"
                    << "to using SAMRAI_MPI objects that are congruent with\n"
                    << "the duplicated object."
                    << std::endl);
      }
   }
#endif

   clusterAndComputeRelationships();


   if (d_sort_output_nodes == true) {
      /*
       * Sorting the node indices is not required.
       * This optional step makes the results order
       * deterministic, which makes the results repeatable.
       * (The natural order of the output of the asynchronous
       * clustering algorithm is non-deterministic because
       * it depends on the order of asynchronous messages.)
       */
      sortOutputBoxes();
   }

   /*
    * Get some global parameters.  Do it before logging to prevent
    * the logging flag from having an undue side effect on performance.
    */
   d_object_timers->t_global_reductions->start();
   d_new_box_level->getGlobalNumberOfBoxes();
   d_new_box_level->getGlobalNumberOfCells();
   for (hier::BoxContainer::const_iterator bi = bound_boxes.begin();
        bi != bound_boxes.end(); ++bi) {
      d_new_box_level->getGlobalBoundingBox(bi->getBlockId().getBlockValue());
   }
   d_object_timers->t_global_reductions->stop();

   if (d_log_cluster) {
      d_object_timers->t_logging->start();
      tbox::plog << "BergerRigoutsos cluster log:\n"
                 << "\tNew box_level clustered by BergerRigoutsos:\n" << d_new_box_level->format("",
                                                                                               2)
                 << "\tBergerRigoutsos tag_to_new:\n" << d_tag_to_new->format("", 2)
                 << "\tBergerRigoutsos new_to_tag:\n" << d_tag_to_new->getTranspose().format("", 2);
      d_object_timers->t_logging->stop();
   }
   if (d_log_cluster_summary) {
      /*
       * Log summary of clustering and dendogram.
       */
      d_object_timers->t_logging->start();
      tbox::plog << "BergerRigoutsos summary:\n"
                 << "\tAsync BR on proc " << mpi.getRank()
                 << " owned "
                 << getMaxOwnership() << " participating in "
                 << getMaxNodes() << " nodes ("
                 << (double)getMaxOwnership() / getMaxNodes()
                 << ") in " << getMaxGeneration() << " generations,"
                 << "   " << getNumBoxesGenerated()
                 << " boxes generated.\n\t"
                 << getMaxTagsOwned() << " locally owned tags on new BoxLevel.\n\t";

      for (hier::BoxContainer::const_iterator bi = bound_boxes.begin();
           bi != bound_boxes.end(); ++bi) {
         const int bn = bi->getBlockId().getBlockValue();
         tbox::plog << "Block " << bn
                    << " initial bounding box = " << *bi << ", "
                    << bi->size() << " cells, "
                    << "final global bounding box = "
                    << d_new_box_level->getGlobalBoundingBox(bn)
                    << ", "
                    << d_new_box_level->getGlobalBoundingBox(bn).size()
                    << " cells.\n\t";
      }

      tbox::plog << "Final output has " << getNumTags()
                 << " tags in "
                 << d_new_box_level->getGlobalNumberOfCells()
                 << " global cells [" << d_new_box_level->getMinNumberOfCells()
                 << "-" << d_new_box_level->getMaxNumberOfCells() << "], "
                 << "over-refinement " << double(d_new_box_level->getGlobalNumberOfCells())/getNumTags()-1 << ", "
                 << d_new_box_level->getGlobalNumberOfBoxes()
                 << " global boxes [" << d_new_box_level->getMinNumberOfBoxes()
                 << "-" << d_new_box_level->getMaxNumberOfBoxes() << "]\n\t"
                 << "Number of continuations: avg = "
                 << getAvgNumberOfCont()
                 << "   max = " << getMaxNumberOfCont() << '\n'
                 << "\tBergerRigoutsos new_level summary:\n" << d_new_box_level->format("\t\t",0)
                 << "\tBergerRigoutsos new_level statistics:\n" << d_new_box_level->formatStatistics("\t\t")
                 << "\tBergerRigoutsos new_to_tag summary:\n" << d_tag_to_new->getTranspose().format("\t\t",0)
                 << "\tBergerRigoutsos new_to_tag statistics:\n" << d_tag_to_new->getTranspose().formatStatistics("\t\t")
                 << "\tBergerRigoutsos tag_to_new summary:\n" << d_tag_to_new->format("\t\t",0)
                 << "\tBergerRigoutsos tag_to_new statistics:\n" << d_tag_to_new->formatStatistics("\t\t")
                 << "\n";
      d_object_timers->t_logging->stop();
   }

   /*
    * Set outputs.  Clear temporary parameters that are only used
    * during active clustering.
    */
   new_box_level = d_new_box_level;
   tag_to_new = d_tag_to_new;
   d_new_box_level.reset();
   d_tag_to_new.reset();
   d_tag_level.reset();


   if (d_barrier_after) {
      d_object_timers->t_barrier_after->start();
      mpi.Barrier();
      d_object_timers->t_barrier_after->stop();
   }

   d_object_timers->t_find_boxes_containing_tags->stop();
}

/*
 ********************************************************************
 ********************************************************************
 */
void
BergerRigoutsos::clusterAndComputeRelationships()
{
   d_object_timers->t_cluster_and_compute_relationships->start();


   /*
    * During the algorithm, we kept the results in primitive
    * containers to avoid the overhead of fine-grain changes to the
    * output objects.  Now initialize the outputs using those
    * primitive containers.
    */

   d_new_box_level.reset(new hier::BoxLevel(
      d_tag_level->getRatioToLevelZero(),
      d_tag_level->getGridGeometry(),
      d_tag_level->getBoxLevel()->getMPI()));

   if (d_compute_relationships >= 1) {
      d_tag_to_new.reset(new hier::Connector(*d_tag_level->getBoxLevel(),
         *d_new_box_level,
         d_max_gcw));
   }
   if (d_compute_relationships >= 2) {
      hier::Connector* new_to_tag =
         new hier::Connector(*d_new_box_level,
                             *d_tag_level->getBoxLevel(),
                             d_max_gcw);
      d_tag_to_new->setTranspose(new_to_tag, true);
   }

   d_object_timers->t_cluster->start();

   /*
    * If compute_relationships == 1:
    *   - Compute relationships from tagged level to new levels.
    *     These relationships are organized around the tagged nodes.
    *     They do not need to be shared with the owners of the
    *     new nodes.
    *
    * If compute_relationships == 2:
    *   - Compute relationships as in compute_relationships == 1 case.
    *   - Owners of new relationships send new relationship data to owners
    *     of new nodes.  This creates the neighbor data
    *     organized around the new nodes.
    */

   if (d_compute_relationships > 0) {

      /*
       * Create empty neighbor lists for nodes on tagged box_level.
       * As new nodes are finalized, they will be added to
       * these lists.
       */
      const hier::BoxContainer& tag_boxes = d_tag_level->getBoxLevel()->getBoxes();
      for (hier::RealBoxConstIterator ni(tag_boxes.realBegin());
           ni != tag_boxes.realEnd(); ++ni) {
         d_tag_to_new->makeEmptyLocalNeighborhood(ni->getBoxId());
      }
      TBOX_ASSERT(
         static_cast<int>(d_tag_level->getBoxLevel()->getLocalNumberOfBoxes()) ==
         d_tag_to_new->getLocalNumberOfNeighborSets());

   }

   TBOX_ASSERT(d_algo_advance_mode == ADVANCE_SOME ||
      d_algo_advance_mode == ADVANCE_ANY ||
      d_algo_advance_mode == SYNCHRONOUS);            // No other supported currently.
   {

      /*
       * Create a BergerRigoutsosNode for each incoming root box and
       * push into the relaunch queue for execution.
       *
       * We use extent data from each incoming root box but do not
       * assume its id is properly set, because the interface did not
       * guarantee they would be.
       */
      hier::LocalId root_box_local_id(0);
      std::list< boost::shared_ptr<BergerRigoutsosNode> > block_nodes_to_delete;
      for (hier::BoxContainer::const_iterator rb = d_root_boxes.begin();
           rb != d_root_boxes.end(); ++rb) {

         const hier::Box block_box(*rb, root_box_local_id, 0);

         BergerRigoutsosNode *block_node(
            new BergerRigoutsosNode( this, block_box ) );

         d_relaunch_queue.push_back(block_node);

         block_nodes_to_delete.push_back(boost::shared_ptr<BergerRigoutsosNode>(block_node));
      }


      int n_comm_group_completed = 0;
      do {

         d_object_timers->t_compute->start();
         while (!d_relaunch_queue.empty()) {
            BergerRigoutsosNode* node_for_relaunch = d_relaunch_queue.front();
            d_relaunch_queue.pop_front();
            if (0) {
               tbox::plog << "Continuing from queue ";
               node_for_relaunch->printNodeState(tbox::plog);
               tbox::plog << std::endl;
            }
            node_for_relaunch->continueAlgorithm();
            if (0) {
               tbox::plog << "Exiting continueAlgorithm ";
               node_for_relaunch->printNodeState(tbox::plog);
               tbox::plog << std::endl;
            }
         }
         d_object_timers->t_compute->stop();

         d_object_timers->t_comm_wait->start();
         n_comm_group_completed =
            static_cast<int>(d_comm_stage.advanceSome());
         d_object_timers->t_comm_wait->stop();

         d_object_timers->t_compute->start();
         while ( d_comm_stage.numberOfCompletedMembers() > 0 ) {
            BergerRigoutsosNode* node_for_relaunch =
               (BergerRigoutsosNode *)(d_comm_stage.popCompletionQueue()->getHandler());
            if (0) {
               tbox::plog << "Continuing from stage ";
               node_for_relaunch->printNodeState(tbox::plog);
               tbox::plog << std::endl;
            }
            node_for_relaunch->continueAlgorithm();
            if (0) {
               tbox::plog << "Exiting continueAlgorithm ";
               node_for_relaunch->printNodeState(tbox::plog);
               tbox::plog << std::endl;
            }
         }
         d_object_timers->t_compute->stop();

         if (0) {
            tbox::plog << "relaunch_queue size "
                       << d_relaunch_queue.size()
                       << "   groups completed: " << n_comm_group_completed
                       << std::endl;
            tbox::plog << "Stage has " << d_comm_stage.numberOfMembers()
                       << " members, "
                       << d_comm_stage.numberOfPendingMembers()
                       << " pending members, "
                       << d_comm_stage.numberOfPendingRequests()
                       << " pending requests." << std::endl;
         }
      } while ( !d_relaunch_queue.empty() || d_comm_stage.hasPendingRequests() );

   }

   TBOX_ASSERT( d_relaunch_queue.empty() );
   TBOX_ASSERT( !d_comm_stage.hasPendingRequests() );

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_compute_relationships > 2) {
      // Each new node should have its own neighbor list.
      TBOX_ASSERT(d_new_box_level->getBoxes().size() ==
                  d_tag_to_new->getTranspose().getLocalNumberOfNeighborSets());
   }
#endif

   // Barrier to separate clustering cost from relationship sharing cost.
   d_mpi_object.Barrier();

   d_object_timers->t_cluster->stop();

   /*
    * Share relationships with owners, if requested.
    * This is a one-time operation that is not considered a part
    * of continueAlgorithm(), so it lies outside that timimg.
    */
   if (d_compute_relationships > 1) {
      shareNewNeighborhoodSetsWithOwners();
   }

   // Clear out communication data or it will mess up next clustering run.
   d_relationship_senders.clear();
   d_relationship_messages.clear();

   d_object_timers->t_cluster_and_compute_relationships->stop();

   d_new_box_level->finalize();

   TBOX_ASSERT(d_tag_to_new->checkConsistencyWithBase() == 0);
   TBOX_ASSERT(d_tag_to_new->checkConsistencyWithHead() == 0);
   TBOX_ASSERT(d_tag_to_new->getTranspose().checkConsistencyWithBase() == 0);
   TBOX_ASSERT(d_tag_to_new->getTranspose().checkConsistencyWithHead() == 0);

#ifdef DEBUG_CHECK_ASSERTIONS
   assertNoMessageForPrivateCommunicator();
#endif

   if ( d_mpi_object == d_tag_level->getBoxLevel()->getMPI() ) {
      /*
       * We have been using an external SAMRAI_MPI.
       * Reset it to avoid mistaking it for an internal one.
       */
      d_mpi_object.setCommunicator(MPI_COMM_NULL);
   }

}

/*
 **********************************************************************
 *
 * Send new relationships found by local process to owners of the new nodes
 * associated with those relationships.  Receive similar data from other
 * processes.
 *
 * Messages to be sent out were placed in d_relationship_messages by
 * computeNewNeighborhoodSets().  This method sends out these messages
 * and receives anticipated messages from processes listed in
 * d_relationship_senders.  Received messages are unpacked to get
 * data on new relationships.
 *
 **********************************************************************
 */
void
BergerRigoutsos::shareNewNeighborhoodSetsWithOwners()
{
   tbox::SAMRAI_MPI mpi(d_mpi_object);
   if (mpi.getSize() == 1) {
      return;
   }

   d_object_timers->t_share_new_relationships->start();

   IntSet relationship_senders = d_relationship_senders;
   std::map<int, VectorOfInts>& relationship_messages = d_relationship_messages;

   const int ints_per_node = hier::Box::commBufferSize(getDim());

   int ierr;
   tbox::SAMRAI_MPI::Status mpi_status;

   // Nonblocking send of relationship data.
   d_object_timers->t_share_new_relationships_send->start();
   tbox::Array<tbox::SAMRAI_MPI::Request> mpi_request(
      static_cast<int>(relationship_messages.size()));
   std::map<int, VectorOfInts>::iterator send_i;
   int nsend = 0;
   for (send_i = relationship_messages.begin(), nsend = 0;
        send_i != relationship_messages.end();
        ++send_i, ++nsend) {
      const int& owner = (*send_i).first;
      VectorOfInts& msg = (*send_i).second;
      ierr = mpi.Isend(&msg[0],
            static_cast<int>(msg.size()),
            MPI_INT,
            owner,
            d_tag_upper_bound,
            &mpi_request[nsend]);
#ifndef DEBUG_CHECK_ASSERTIONS
      NULL_USE(ierr);
#endif
      TBOX_ASSERT(ierr == MPI_SUCCESS);
   }
   d_object_timers->t_share_new_relationships_send->stop();

   {
      /*
       * The rest of this method assumes current process is NOT
       * in relationship_senders, so remove it.  For efficiency, method
       * computeNewNeighborhoodSets() (which created the relationship senders)
       * did not remove it.
       */
      IntSet::iterator local = relationship_senders.find(d_mpi_object.getRank());
      if (local != relationship_senders.end()) {
         relationship_senders.erase(local);
      }
   }

   /*
    * Create set recved_from which is to contain ranks of
    * processes from which we've received the expected relationship data.
    * The while loop goes until all expected messages have
    * been received from relationship_senders.
    *
    * In the while loop:
    *    - Probe for an incomming message.
    *    - Determine its size allocate memory for receiving the message.
    *    - Receive the message.
    *    - Get relationship data from the message.
    */
   IntSet recved_from;
   while (recved_from.size() < relationship_senders.size()) {

      d_object_timers->t_share_new_relationships_recv->start();
      ierr = mpi.Probe(MPI_ANY_SOURCE,
            d_tag_upper_bound,
            &mpi_status);
      TBOX_ASSERT(ierr == MPI_SUCCESS);

      const int sender = mpi_status.MPI_SOURCE;
      int mesg_size = -1;
      mpi.Get_count(&mpi_status, MPI_INT, &mesg_size);
      TBOX_ASSERT(relationship_senders.find(sender) != relationship_senders.end());
      TBOX_ASSERT(recved_from.find(sender) == recved_from.end());
      TBOX_ASSERT(mesg_size >= 0);

      tbox::Array<int> buf(mesg_size);
      int* ptr = buf.getPointer();
      ierr = mpi.Recv(ptr,
            mesg_size,
            MPI_INT,
            sender,
            d_tag_upper_bound,
            &mpi_status);
      TBOX_ASSERT(ierr == MPI_SUCCESS);
      d_object_timers->t_share_new_relationships_recv->stop();

      d_object_timers->t_share_new_relationships_unpack->start();
      int consumed = 0;
      while (ptr < buf.getPointer() + buf.size()) {
         const hier::LocalId new_local_id(*(ptr++));
         hier::BoxId box_id(new_local_id, d_mpi_object.getRank());
         int n_new_relationships = *(ptr++);
         TBOX_ASSERT(d_tag_to_new->getTranspose().hasNeighborSet(box_id));
         if (n_new_relationships > 0) {
            hier::Connector::NeighborhoodIterator base_box_itr =
               d_tag_to_new->getTranspose().makeEmptyLocalNeighborhood(box_id);
            for (int n = 0; n < n_new_relationships; ++n) {
               hier::Box node(getDim());
               node.getFromIntBuffer(ptr);
               ptr += ints_per_node;
               d_tag_to_new->getTranspose().insertLocalNeighbor(node, base_box_itr);
            }
         }
         consumed += 2 + n_new_relationships * ints_per_node;
      }
      recved_from.insert(sender);
      d_object_timers->t_share_new_relationships_unpack->stop();
   }

   if (nsend > 0) {
      // Make sure all nonblocking sends completed.
      d_object_timers->t_share_new_relationships_send->start();
      tbox::Array<tbox::SAMRAI_MPI::Status> mpi_statuses(
         static_cast<int>(relationship_messages.size()));
      ierr = mpi.Waitall(static_cast<int>(relationship_messages.size()),
            mpi_request.getPointer(),
            mpi_statuses.getPointer());
      TBOX_ASSERT(ierr == MPI_SUCCESS);
      d_object_timers->t_share_new_relationships_send->stop();
   }

   d_object_timers->t_share_new_relationships->stop();

}



/*
 **********************************************************************
 *
 * Methods for setting algorithm parameters before running.
 *
 **********************************************************************
 */

void
BergerRigoutsos::setAlgorithmAdvanceMode(
   const std::string& mode)
{
   if (mode == "ADVANCE_ANY") {
      d_algo_advance_mode = ADVANCE_ANY;
   } else if (mode == "ADVANCE_SOME") {
      d_algo_advance_mode = ADVANCE_SOME;
   } else if (mode == "SYNCHRONOUS") {
      d_algo_advance_mode = SYNCHRONOUS;
   } else {
      TBOX_ERROR("No such algorithm choice: " << mode << "\n");
   }
}

void
BergerRigoutsos::setOwnerMode(
   const std::string& mode)
{
   if (mode == "SINGLE_OWNER") {
      d_owner_mode = SINGLE_OWNER;
   } else if (mode == "MOST_OVERLAP") {
      d_owner_mode = MOST_OVERLAP;
   } else if (mode == "FEWEST_OWNED") {
      d_owner_mode = FEWEST_OWNED;
   } else if (mode == "LEAST_ACTIVE") {
      d_owner_mode = LEAST_ACTIVE;
   } else {
      TBOX_ERROR("BergerRigoutsos: Unrecognized owner mode request: "
         << mode << std::endl);
   }
}

void
BergerRigoutsos::setComputeRelationships(
   const std::string mode,
   const hier::IntVector& ghost_cell_width)
{
   if (mode == "NONE") {
      d_compute_relationships = 0;
   } else if (mode == "TAG_TO_NEW") {
      d_compute_relationships = 1;
   } else if (mode == "BIDIRECTIONAL") {
      d_compute_relationships = 2;
   } else {
      TBOX_ERROR("BergerRigoutsos::setComputeRelationships error:\n"
         << "bad mode '" << mode << "' specified.\n"
         << "Should be one of NONE, TAG_TO_NEW, BIDIRECTIONAL" << std::endl);
   }
   TBOX_ASSERT(ghost_cell_width >= hier::IntVector::getZero(ghost_cell_width.getDim()));
   d_max_gcw = ghost_cell_width;
}


/*
 **************************************************************************
 **************************************************************************
 */
void
BergerRigoutsos::resetCounters()
{
   d_num_tags_in_all_nodes = 0;
   d_max_tags_owned = 0;
   d_num_nodes_allocated = 0;
   d_max_nodes_allocated = 0;
   d_num_nodes_active = 0;
   d_max_nodes_active = 0;
   d_num_nodes_owned = 0;
   d_max_nodes_owned = 0;
   d_num_nodes_commwait = 0;
   d_max_nodes_commwait = 0;
   d_num_nodes_completed = 0;
   d_max_generation = 0;
   d_num_boxes_generated = 0;
   d_num_conts_to_complete = 0;
   d_max_conts_to_complete = 0;
   d_num_nodes_existing = 0;
}

/*
 **************************************************************************
 **************************************************************************
 */
void
BergerRigoutsos::useDuplicateMPI(
   const tbox::SAMRAI_MPI& mpi_object)
{
   TBOX_ASSERT( !d_tag_level ); // Setting MPI during clustering makes a mess.

   // If needed, free current private communicator.
   if ( d_mpi_object.getCommunicator() != MPI_COMM_NULL ) {
      d_mpi_object.freeCommunicator();
      TBOX_ASSERT( d_mpi_object.getCommunicator() == MPI_COMM_NULL );
   }

   if (mpi_object.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      d_mpi_object.dupCommunicator(mpi_object);
   }

   setupMPIDependentData();
}

/*
 **************************************************************************
 * Check the congruency between d_mpi and d_tag_level's MPI.
 * Writes out warning in log if not congruent.
 * Returns whether the two are congruent.
 *
 * Note: Sequential runs (no MPI or MPI with 1 process) always means
 * congruency.  d_mpi with a Null communicator indicates that we will
 * use the d_tag_level's MPI, so that is also automatically congruent.
 **************************************************************************
 */
bool
BergerRigoutsos::checkMPICongruency() const
{

   if ( !tbox::SAMRAI_MPI::usingMPI() ||
        ( d_mpi_object.getCommunicator() == MPI_COMM_NULL ) ||
        ( d_mpi_object.getSize() == 1 &&
          d_tag_level->getBoxLevel()->getMPI().getSize() == 1 ) ) {
      return true;
   }

   /*
    * If a valid MPI communicator is given, use it instead of the
    * tag BoxLevel's communicator.  It must be congruent with
    * the tag BoxLevel's.
    */

   bool is_congruent = true;
   /*
    * Make sure mpi_object is compatible with the BoxLevel
    * involved.
    */
   tbox::SAMRAI_MPI mpi1(d_mpi_object);
   tbox::SAMRAI_MPI mpi2(d_tag_level->getBoxLevel()->getMPI());
   TBOX_ASSERT(mpi1.getSize() == mpi2.getSize());
   TBOX_ASSERT(mpi1.getRank() == mpi2.getRank());
   if (mpi1.getSize() > 1) {
      int compare_result;
      tbox::SAMRAI_MPI::Comm_compare(
         d_mpi_object.getCommunicator(),
         d_tag_level->getBoxLevel()->getMPI().getCommunicator(),
         &compare_result);
      is_congruent =
         (compare_result == MPI_CONGRUENT) ||
         (compare_result == MPI_IDENT);
   }

   return is_congruent;
}

/*
 **************************************************************************
 **************************************************************************
 */
void
BergerRigoutsos::setupMPIDependentData()
{
   /*
    * Reserve the tag upper bound for the relationship-sharing phase.
    * Divide the rest into tag pools divided among all processes.
    */
   if (tbox::SAMRAI_MPI::usingMPI()) {
      /*
       * For some MPI implementations, I cannot get the attribute for
       * any communicator except for MPI_COMM_WORLD.  Assuming the tag
       * upper bound is the same for all communicators, I will try
       * some other communicators to get it.
       */
      int* tag_upper_bound_ptr, flag;
      d_mpi_object.Attr_get(
         MPI_TAG_UB,
         &tag_upper_bound_ptr,
         &flag);
      if (tag_upper_bound_ptr == 0) {
         tbox::SAMRAI_MPI::getSAMRAIWorld().Attr_get(
            MPI_TAG_UB,
            &tag_upper_bound_ptr,
            &flag);
      }
      if (tag_upper_bound_ptr == 0) {
         tbox::SAMRAI_MPI mpi1(tbox::SAMRAI_MPI::commWorld);
         mpi1.Attr_get(
            MPI_TAG_UB,
            &tag_upper_bound_ptr,
            &flag);
      }
      TBOX_ASSERT(tag_upper_bound_ptr != 0);
      d_tag_upper_bound = *tag_upper_bound_ptr;

   } else {
      // MPI not used, so choose a sufficiently big tag upper bound.
      d_tag_upper_bound = 1000000;
   }

   // Divide the rest into tag pools divided among all processes.
   d_available_mpi_tag =
      d_tag_upper_bound / d_mpi_object.getSize() * d_mpi_object.getRank();

}

/*
***********************************************************************
***********************************************************************
*/
void
BergerRigoutsos::sortOutputBoxes()
{

   d_object_timers->t_sort_output_nodes->start();

   TBOX_ASSERT(d_tag_to_new->hasTranspose());
   hier::Connector& new_to_tag = d_tag_to_new->getTranspose();

   if (0) {
      // Check inputs.
      int errs = 0;
      if (d_tag_to_new->checkOverlapCorrectness(false, true)) {
         ++errs;
         tbox::perr << "Error found in tag_to_new!\n";
      }
      if (new_to_tag.checkOverlapCorrectness(false, true)) {
         ++errs;
         tbox::perr << "Error found in new_to_tag!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors found before sorting nodes."
            << "new_box_level:\n" << d_new_box_level->format("", 2)
            << "tag box_level:\n" << d_tag_to_new->getBase().format("", 2)
            << "tag_to_new:\n" << d_tag_to_new->format("", 2)
            << "new_to_tag:\n" << new_to_tag.format("", 2) << std::endl);
      }
   }

   /*
    * Sort local indices by corners to make the output deterministic.
    */
   boost::shared_ptr<hier::MappingConnector> sorting_map;
   boost::shared_ptr<hier::BoxLevel> sorted_box_level;
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      sorted_box_level,
      sorting_map,
      *d_new_box_level,
      true /* sort nodes by corners */,
      false /* don't sequentialize indices globally */);
   if (0) {
      tbox::plog
         << "tag box_level:\n" << d_tag_to_new->getBase().format("", 2)
         << "tag_to_new:\n" << d_tag_to_new->format("", 2)
         << "new_to_tag:\n" << new_to_tag.format("", 2)
         << "Sorting map:\n" << sorting_map->format("", 2);
   }
   if (0) {
      // Check sorting_map before using it.
      int errs = 0;
      if (sorting_map->checkOverlapCorrectness(false, true)) {
         ++errs;
         tbox::perr << "Error found in sorting_map!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "presorted box_level:\n" << d_new_box_level->format("", 2)
            << "sorted box_level:\n" << sorted_box_level->format("", 2)
            << "sorting_map:\n" << sorting_map->format("", 2) << std::endl);
      }
   }
   hier::MappingConnectorAlgorithm mca;
   mca.modify(*d_tag_to_new,
              *sorting_map,
              d_new_box_level.get());
   if (0) {
      // Check result of mapping.
      int errs = 0;
      if (d_tag_to_new->checkOverlapCorrectness(false, true)) {
         ++errs;
         tbox::perr << "Error found in tag_to_new!\n";
      }
      if (new_to_tag.checkOverlapCorrectness(false, true)) {
         ++errs;
         tbox::perr << "Error found in new_to_tag!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors found after sorting nodes."
            << "new_box_level:\n" << d_new_box_level->format("", 2)
            << "tag box_level:\n" << d_tag_to_new->getBase().format("", 2)
            << "tag_to_new:\n" << d_tag_to_new->format("", 2)
            << "new_to_tag:\n" << new_to_tag.format("", 2) << std::endl);
      }
   }

   d_object_timers->t_sort_output_nodes->stop();
}

/*
***************************************************************************
*
***************************************************************************
*/
void
BergerRigoutsos::assertNoMessageForPrivateCommunicator() const
{
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if (d_mpi_object.getCommunicator() != tbox::SAMRAI_MPI::commNull &&
       d_mpi_object != d_tag_level->getBoxLevel()->getMPI() ) {
      int flag;
      tbox::SAMRAI_MPI::Status mpi_status;
      int mpi_err = d_mpi_object.Iprobe(MPI_ANY_SOURCE,
                                        MPI_ANY_TAG,
                                        &flag,
                                        &mpi_status);
      if (mpi_err != MPI_SUCCESS) {
         TBOX_ERROR("Error probing for possible lost messages." << std::endl);
      }
      if (flag == true) {
         int count = -1;
         mpi_err = tbox::SAMRAI_MPI::Get_count(&mpi_status, MPI_INT, &count);
         TBOX_ERROR("Library error!\n"
                    << "BergerRigoutsos detected before or after\n"
                    << "the clustering algorithm that there\n"
                    << "is a message yet to be received.  This is\n"
                    << "an error because all messages using the\n"
                    << "private communicator should have been\n"
                    << "accounted for.  Message status:\n"
                    << "source " << mpi_status.MPI_SOURCE << '\n'
                    << "tag " << mpi_status.MPI_TAG << '\n'
                    << "count " << count << " (assuming integers)\n");
      }
   }
}


/*
 ***********************************************************************
 ***********************************************************************
 */
void
BergerRigoutsos::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix));
   if (ti != s_static_timers.end()) {
      d_object_timers = &(ti->second);
   } else {
      setObjectTimers(timer_prefix);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BergerRigoutsos::setObjectTimers(
   const std::string& timer_prefix)
{
   d_object_timers = &s_static_timers[timer_prefix];

   tbox::TimerManager *tm = tbox::TimerManager::getManager();

   d_object_timers->t_find_boxes_containing_tags = tm->
      getTimer("mesh::BergerRigoutsos::findBoxesContainingTags()");

   d_object_timers->t_cluster = tm->
      getTimer(timer_prefix + "::cluster");
   d_object_timers->t_cluster_and_compute_relationships = tm->
      getTimer(timer_prefix + "::clusterAndComputeRelationships()");
   d_object_timers->t_continue_algorithm = tm->
      getTimer(timer_prefix + "::continueAlgorithm()");

   d_object_timers->t_compute = tm->
      getTimer(timer_prefix + "::compute");
   d_object_timers->t_comm_wait = tm->
      getTimer(timer_prefix + "::Comm_wait");
   d_object_timers->t_MPI_wait = tm->
      getTimer(timer_prefix + "::MPI_wait");

   d_object_timers->t_compute_new_graph_relationships = tm->
      getTimer(timer_prefix + "::computeNewNeighborhoodSets()");
   d_object_timers->t_share_new_relationships = tm->
      getTimer(timer_prefix + "::shareNewNeighborhoodSetsWithOwners()");
   d_object_timers->t_share_new_relationships_send = tm->
      getTimer(timer_prefix + "::shareNewNeighborhoodSetsWithOwners()_send");
   d_object_timers->t_share_new_relationships_recv = tm->
      getTimer(timer_prefix + "::shareNewNeighborhoodSetsWithOwners()_recv");
   d_object_timers->t_share_new_relationships_unpack = tm->
      getTimer(timer_prefix + "::shareNewNeighborhoodSetsWithOwners()_unpack");

   d_object_timers->t_local_histogram = tm->
      getTimer(timer_prefix + "::makeLocalTagHistogram()");
   d_object_timers->t_local_tasks = tm->
      getTimer(timer_prefix + "::continueAlgorithm()_local_tasks");

   // Multi-stage timers
   d_object_timers->t_reduce_histogram = tm->
      getTimer(timer_prefix + "::reduce_histogram");
   d_object_timers->t_bcast_acceptability = tm->
      getTimer(timer_prefix + "::bcast_acceptability");
   d_object_timers->t_gather_grouping_criteria = tm->
      getTimer(timer_prefix + "::gather_grouping_criteria");
   d_object_timers->t_bcast_child_groups = tm->
      getTimer(timer_prefix + "::bcast_child_groups");
   d_object_timers->t_bcast_to_dropouts = tm->
      getTimer(timer_prefix + "::bcast_to_dropouts");

   // Pre- and post-processing timers.
   d_object_timers->t_barrier_before = tm->
      getTimer("mesh::BergerRigoutsos::barrier_before");
   d_object_timers->t_barrier_after = tm->
      getTimer("mesh::BergerRigoutsos::barrier_after");
   d_object_timers->t_global_reductions = tm->
      getTimer("mesh::BergerRigoutsos::global_reductions");
   d_object_timers->t_sort_output_nodes = tm->
      getTimer("mesh::BergerRigoutsos::sort_output_nodes");
   d_object_timers->t_logging = tm->
      getTimer("mesh::BergerRigoutsos::logging");
}


}
}
#endif
