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

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/BergerRigoutsosNode.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

namespace SAMRAI {
namespace mesh {

boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_barrier_before;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_barrier_after;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_find_boxes_with_tags;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_cluster_and_compute_relationships;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_run_abr;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_global_reductions;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_logging;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_sort_output_nodes;

tbox::StartupShutdownManager::Handler
BergerRigoutsos::s_initialize_finalize_handler(
   BergerRigoutsos::initializeCallback,
   0,
   0,
   BergerRigoutsos::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

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

#if 0
/*
 *************************************************************************
 * Set the MPI communicator.
 *************************************************************************
 */
void
BergerRigoutsos::useDuplicateMPI(
   const tbox::SAMRAI_MPI& mpi)
{
   d_common.useDuplicateMPI(mpi);
}
#endif

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
      t_barrier_before->start();
      mpi.Barrier();
      t_barrier_before->stop();
   }

   for ( hier::BoxContainer::const_iterator bi=bound_boxes.begin();
         bi!=bound_boxes.end(); ++bi ) {
      if (bi->empty()) {
         TBOX_ERROR("BergerRigoutsos: empty bounding box not allowed.");
      }
   }

   const hier::BoxLevel& tag_box_level = *tag_level->getBoxLevel();

   t_find_boxes_with_tags->start();

   setParameters(
      tag_data_index,
      tag_val,
      min_box,
      efficiency_tol,
      combine_tol,
      d_max_box_size,
      d_max_inflection_cut_from_center,
      d_inflection_cut_threshold_ar);

   // Set the parallel algorithm.
   // setAlgorithmAdvanceMode(d_algo_advance_mode);
   // setOwnerMode(d_owner_mode);
   setComputeRelationships("BIDIRECTIONAL", max_gcw);
   setMinBoxSizeFromCutting(d_min_box_size_from_cutting);

   // Set debugging/verbosity parameters.
   setLogNodeHistory(d_log_node_history);

   t_cluster_and_compute_relationships->start();
   clusterAndComputeRelationships(new_box_level,
                                           tag_to_new,
                                           tag_level,
                                           bound_boxes);
   t_cluster_and_compute_relationships->stop();

   if (d_sort_output_nodes == true) {
      /*
       * Sorting the node indices is not required.
       * This optional step makes the results order
       * deterministic, which makes the results repeatable.
       * (The natural order of the output of the asynchronous
       * clustering algorithm is non-deterministic because
       * it depends on the order of asynchronous messages.)
       */
      sortOutputBoxes(*new_box_level,
                      *tag_to_new);
   }

   /*
    * Get some global parameters.  Do it before logging to prevent
    * the logging flag from having an undue side effect on performance.
    */
   t_global_reductions->start();
   new_box_level->getGlobalNumberOfBoxes();
   new_box_level->getGlobalNumberOfCells();
   for (hier::BoxContainer::const_iterator bi = bound_boxes.begin();
        bi != bound_boxes.end(); ++bi) {
      new_box_level->getGlobalBoundingBox(bi->getBlockId().getBlockValue());
   }
   t_global_reductions->stop();

   if (d_log_cluster) {
      t_logging->start();
      tbox::plog << "BergerRigoutsos cluster log:\n"
                 << "\tNew box_level clustered by BergerRigoutsos:\n" << new_box_level->format("",
                                                                                               2)
                 << "\tBergerRigoutsos tag_to_new:\n" << tag_to_new->format("", 2)
                 << "\tBergerRigoutsos new_to_tag:\n" << tag_to_new->getTranspose().format("", 2);
      t_logging->stop();
   }
   if (d_log_cluster_summary) {
      /*
       * Log summary of clustering and dendogram.
       */
      t_logging->start();
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
                    << new_box_level->getGlobalBoundingBox(bn)
                    << ", "
                    << new_box_level->getGlobalBoundingBox(bn).size()
                    << " cells.\n\t";
      }

      tbox::plog << "Final output has " << getNumTags()
                 << " tags in "
                 << new_box_level->getGlobalNumberOfCells()
                 << " global cells [" << new_box_level->getMinNumberOfCells()
                 << "-" << new_box_level->getMaxNumberOfCells() << "], "
                 << "over-refinement " << double(new_box_level->getGlobalNumberOfCells())/getNumTags()-1 << ", "
                 << new_box_level->getGlobalNumberOfBoxes()
                 << " global boxes [" << new_box_level->getMinNumberOfBoxes()
                 << "-" << new_box_level->getMaxNumberOfBoxes() << "]\n\t"
                 << "Number of continuations: avg = "
                 << getAvgNumberOfCont()
                 << "   max = " << getMaxNumberOfCont() << '\n'
                 << "\tBergerRigoutsos new_level summary:\n" << new_box_level->format("\t\t",0)
                 << "\tBergerRigoutsos new_level statistics:\n" << new_box_level->formatStatistics("\t\t")
                 << "\tBergerRigoutsos new_to_tag summary:\n" << tag_to_new->getTranspose().format("\t\t",0)
                 << "\tBergerRigoutsos new_to_tag statistics:\n" << tag_to_new->getTranspose().formatStatistics("\t\t")
                 << "\tBergerRigoutsos tag_to_new summary:\n" << tag_to_new->format("\t\t",0)
                 << "\tBergerRigoutsos tag_to_new statistics:\n" << tag_to_new->formatStatistics("\t\t")
                 << "\n";
      t_logging->stop();
   }

   if (d_barrier_after) {
      t_barrier_after->start();
      mpi.Barrier();
      t_barrier_after->stop();
   }

   t_find_boxes_with_tags->stop();
}

/*
***********************************************************************
***********************************************************************
*/
void
BergerRigoutsos::sortOutputBoxes(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new) const
{

   t_sort_output_nodes->start();

   TBOX_ASSERT(tag_to_new.hasTranspose());
   hier::Connector& new_to_tag = tag_to_new.getTranspose();

   if (0) {
      // Check inputs.
      int errs = 0;
      if (tag_to_new.checkOverlapCorrectness(false, true)) {
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
            << "new_box_level:\n" << new_box_level.format("", 2)
            << "tag box_level:\n" << tag_to_new.getBase().format("", 2)
            << "tag_to_new:\n" << tag_to_new.format("", 2)
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
      new_box_level,
      true /* sort nodes by corners */,
      false /* don't sequentialize indices globally */);
   if (0) {
      tbox::plog
         << "tag box_level:\n" << tag_to_new.getBase().format("", 2)
         << "tag_to_new:\n" << tag_to_new.format("", 2)
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
            << "presorted box_level:\n" << new_box_level.format("", 2)
            << "sorted box_level:\n" << sorted_box_level->format("", 2)
            << "sorting_map:\n" << sorting_map->format("", 2) << std::endl);
      }
   }
   hier::MappingConnectorAlgorithm mca;
   mca.modify(tag_to_new,
              *sorting_map,
              &new_box_level);
   if (0) {
      // Check result of mapping.
      int errs = 0;
      if (tag_to_new.checkOverlapCorrectness(false, true)) {
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
            << "new_box_level:\n" << new_box_level.format("", 2)
            << "tag box_level:\n" << tag_to_new.getBase().format("", 2)
            << "tag_to_new:\n" << tag_to_new.format("", 2)
            << "new_to_tag:\n" << new_to_tag.format("", 2) << std::endl);
      }
   }

   t_sort_output_nodes->stop();
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
                    << "running BergerRigoutsosNode that there\n"
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
BergerRigoutsos::initializeCallback()
{
   TBOX_ASSERT(!t_global_reductions);
   t_run_abr = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::run_abr");
   t_find_boxes_with_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::find_boxes_with_tags");
   t_cluster_and_compute_relationships = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::cluster_and_compute_relationships");
   t_global_reductions = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::global_reductions");
   t_logging = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::logging");
   t_sort_output_nodes = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::sort_output_nodes");
   t_barrier_before = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::barrier_before");
   t_barrier_after = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::barrier_after");
}

/*
***************************************************************************
*
* Release static timers.  To be called by shutdown registry to make sure
* memory for timers does not leak.
*
***************************************************************************
*/
void
BergerRigoutsos::finalizeCallback()
{
   t_barrier_before.reset();
   t_barrier_after.reset();
   t_find_boxes_with_tags.reset();
   t_run_abr.reset();
   t_global_reductions.reset();
   t_sort_output_nodes.reset();
}

}
}
#endif
