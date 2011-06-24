/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines. 
 *
 ************************************************************************/

#ifndef included_mesh_MultiblockGriddingAlgorithm_C
#define included_mesh_MultiblockGriddingAlgorithm_C

#include "SAMRAI/mesh/MultiblockGriddingAlgorithm.h"

#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/math/PatchCellDataOpsInteger.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/MultiblockRefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/MappedBoxContainerUtils.h"
#include "SAMRAI/hier/MappedBoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MappedBoxTree.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/GriddingAlgorithmConnectorWidthRequestor.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineScheduleConnectorWidthRequestor.h"

#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

#ifndef SAMRAI_INLINE
#include "SAMRAI/mesh/MultiblockGriddingAlgorithm.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int MultiblockGriddingAlgorithm::ALGS_GRIDDING_ALGORITHM_VERSION = 3;

char MultiblockGriddingAlgorithm::s_check_overflow_nesting = 0;
char MultiblockGriddingAlgorithm::s_check_proper_nesting = 0;
char MultiblockGriddingAlgorithm::s_check_connectors = 0;
char MultiblockGriddingAlgorithm::s_print_mapped_box_level_hierarchy = 0;
char MultiblockGriddingAlgorithm::s_print_steps = 0;



/*
 *************************************************************************
 *
 * Initialize the static data members.
 *
 *************************************************************************
 */
tbox::Array<int> * MultiblockGriddingAlgorithm::s_tag_indx = 0;
tbox::Array<int> * MultiblockGriddingAlgorithm::s_buf_tag_indx = 0;

tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_find_domain_complement;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_load_balance;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_load_balance0;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_load_balance_setup;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bdry_fill_tags_create;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_coarsest;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_finer;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_finer_setup;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_finer_tagging;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_finer_create;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_regrid_all_finer;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_regrid_finer_create;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bridge_links;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_fill_tags_from_mapped_box_level;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_tag_cells_for_refinement;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_buffer_tags;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bdry_fill_tags_comm;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_second_finer_tagging;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_find_refinement;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bridge_new_to_new;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_find_new_to_new;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bridge_new_to_coarser;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bridge_new_to_finer;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_bridge_new_to_old;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_find_boxes_containing_tags;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_enforce_nesting;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_nesting_map;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_nesting_map_compute;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_nesting_map_convert;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_use_nesting_map;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_overflow_map;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_overflow_map_compute;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_overflow_map_convert;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_use_overflow_map;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_compute_external_parts;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_compute_nesting_violator;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_extend_to_domain_boundary;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_extend_within_domain;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_grow_boxes_within_domain;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_sort_nodes;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_modify_connector;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_misc1;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_misc2;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_misc3;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_misc4;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_misc5;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_domain;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_get_balance;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_use_balance;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_make_new;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_process_error;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_limit_overflow;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_reset_hier;
tbox::Pointer<tbox::Timer> MultiblockGriddingAlgorithm::t_box_massage;

tbox::StartupShutdownManager::Handler
MultiblockGriddingAlgorithm::s_startup_shutdown_handler(
   0,
   MultiblockGriddingAlgorithm::startupCallback,
   MultiblockGriddingAlgorithm::shutdownCallback,
   0,
   tbox::StartupShutdownManager::priorityListElements);

tbox::StartupShutdownManager::Handler
MultiblockGriddingAlgorithm::s_initialize_handler(
   MultiblockGriddingAlgorithm::initializeCallback,
   0,
   0,
   MultiblockGriddingAlgorithm::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *
 * Constructor and destructor for MultiblockGriddingAlgorithm.
 *
 *************************************************************************
 */
MultiblockGriddingAlgorithm::MultiblockGriddingAlgorithm(
   const tbox::Pointer<hier::PatchHierarchy> &mb_hierarchy,
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer<TagAndInitializeStrategy> tag_init_strategy,
   tbox::Pointer<BoxGeneratorStrategy> generator,
   tbox::Pointer<LoadBalanceStrategy> balancer,
   tbox::Pointer<LoadBalanceStrategy> balancer0,
   MultiblockGriddingTagger* mb_tagger_strategy,
   bool register_for_restart):
   BaseGriddingAlgorithm(mb_hierarchy),
   d_mb_hierarchy(mb_hierarchy),
   d_dim(mb_hierarchy->getDim()),
   d_domain_complement_tree(0, hier::MappedBoxTree(d_dim)),
   d_check_nonrefined_tags('w'),
   d_check_overlapping_patches('i'),
   d_sequentialize_patch_indices(false),
   d_enforce_proper_nesting(true),
   d_extend_to_domain_boundary(true),
   d_load_balance(true)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!tag_init_strategy.isNull());
   TBOX_ASSERT(!generator.isNull());
   TBOX_ASSERT(!balancer.isNull());

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(d_object_name, this);
   }

   d_tag_init_strategy = tag_init_strategy;
   d_box_generator = generator;
   d_load_balancer = balancer;
   d_load_balancer0 = balancer0;
   if (d_load_balancer0.isNull()) {
      d_load_balancer0 = d_load_balancer;
   }
   d_mb_tagger_strategy = mb_tagger_strategy;

   if (mb_tagger_strategy) {
      d_mb_tagger_strategy = mb_tagger_strategy;
      d_internal_tagger_strategy = false;
   } else {
      d_mb_tagger_strategy = new MultiblockGriddingTagger(d_dim);
      d_internal_tagger_strategy = true;
   }

   /*
    * Construct integer tag variables and add to variable database.
    * Note that variables and patch data indices are shared among
    * all MultiblockGriddingAlgorithm instances.  The VariableDatabase
    * holds the variables, once contructed and registered via the
    * VariableDatabase::makeInternalSAMRAIWorkVariablePatchDataIndex()
    * function call.  Note that variables are registered and patch data
    * indices are made only for the first time through the constructor.
    */
   hier::VariableDatabase* var_db =
      hier::VariableDatabase::getDatabase();

   static std::string tag_interior_variable_name(
      "MultiblockGriddingAlgorithm__tag-interior");
   static std::string tag_buffer_variable_name(
      "MultiblockGriddingAlgorithm__tag-buffer");

   d_tag = var_db->getVariable(tag_interior_variable_name);
   if (d_tag.isNull()) {
      d_tag = new pdat::CellVariable<int>(d_dim, tag_interior_variable_name, 1);
   }

   d_buf_tag = var_db->getVariable(tag_buffer_variable_name);
   if (d_buf_tag.isNull()) {
      d_buf_tag = new pdat::CellVariable<int>(d_dim,
                                              tag_buffer_variable_name,
                                              1);
   }

   if ((*s_tag_indx)[d_dim.getValue() - 1] < 0) {
      (*s_tag_indx)[d_dim.getValue() - 1] =
         var_db->registerInternalSAMRAIVariable(d_tag,
            hier::IntVector::getZero(d_dim));
   }
   if ((*s_buf_tag_indx)[d_dim.getValue() - 1] < 0) {
      (*s_buf_tag_indx)[d_dim.getValue() - 1] =
         var_db->registerInternalSAMRAIVariable(d_buf_tag,
            hier::IntVector::getOne(d_dim));
   }

   d_tag_indx = (*s_tag_indx)[d_dim.getValue() - 1];
   d_buf_tag_indx = (*s_buf_tag_indx)[d_dim.getValue() - 1];

   d_mb_tagger_strategy->setScratchTagPatchDataIndex(d_buf_tag_indx);

   /*
    * Tag value for refined cells is one; others are zero.
    */
   d_true_tag = 1;
   d_false_tag = 0;

   /*
    * Initialize communication algorithm for exchanging buffered tag data.
    */
   d_mblk_bdry_fill_tags = new xfer::MultiblockRefineAlgorithm(d_mb_hierarchy,
                                                               d_dim);

   d_mblk_bdry_fill_tags->
   registerRefine(d_buf_tag_indx,
      d_buf_tag_indx,
      d_buf_tag_indx,
      ((tbox::Pointer<hier::RefineOperator>)NULL));


   d_base_ln = -1;

   d_efficiency_tolerance.resizeArray(d_mb_hierarchy->getMaxNumberOfLevels());
   d_combine_efficiency.resizeArray(d_mb_hierarchy->getMaxNumberOfLevels());

   for (int ln = 0; ln < d_mb_hierarchy->getMaxNumberOfLevels(); ln++) {
      d_efficiency_tolerance[ln] = 0.8e0;
      d_combine_efficiency[ln] = 0.8e0;
   }

   //d_proper_nesting_complement.resize(d_mb_hierarchy->getMaxNumberOfLevels(), hier::MappedBoxLevel(d_dim));
   d_proper_nesting_complement.resize(0);
   d_to_nesting_complement.resize(0);
   d_from_nesting_complement.resize(0);

   /*
    * Initialize object with data read from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart && d_registered_for_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

   /*
    * This part must go after getFromInput(), where d_tag_init_strategy
    * has a critical modification.
    */
   if ( d_tag_init_strategy->getErrorCoarsenRatio() > 1 ) {
      tbox::Pointer<StandardTagAndInitialize> std_tag_init(d_tag_init_strategy);
      if ( ! std_tag_init.isNull() ) {
         d_mb_hierarchy->registerConnectorWidthRequestor(
            std_tag_init->getConnectorWidthRequestor() );
      }
   }

   d_bdry_sched_tags.resizeArray(d_mb_hierarchy->getMaxNumberOfLevels());

   d_singleblock_domain_mapped_box_level.resize(d_mb_hierarchy->getGridGeometry()->getNumberBlocks(),
                                                hier::MappedBoxLevel(d_dim));
   for ( int bn=0; bn<d_mb_hierarchy->getGridGeometry()->getNumberBlocks(); ++bn ) {
      makeSingleBlockMappedBoxLevel(
         d_singleblock_domain_mapped_box_level[bn],
         d_mb_hierarchy->getDomainMappedBoxLevel(),
         hier::BlockId(bn) );
   }

#ifdef MGA_RECORD_STATS

   d_boxes_stat.resizeArray(d_mb_hierarchy->getMaxNumberOfLevels());
   d_cells_stat.resizeArray(d_mb_hierarchy->getMaxNumberOfLevels());
   d_timestamp_stat.resizeArray(d_mb_hierarchy->getMaxNumberOfLevels());

   for (int ln = 0; ln < d_mb_hierarchy->getMaxNumberOfLevels(); ++ln) {
      std::string ln_text = tbox::Utilities::intToString(ln, 2);
      const std::string num_boxes_str = std::string("MGA_BoxesL") + ln_text;
      const std::string num_cells_str = std::string("MGA_CellsL") + ln_text;
      const std::string timestamp_str = std::string("MGA_TimeL") + ln_text;
      d_boxes_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(num_boxes_str, "PROC_STAT");
      d_cells_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(num_cells_str, "PROC_STAT");
      d_timestamp_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(timestamp_str, "PROC_STAT");
   }

#endif

}

/*
 *************************************************************************
 *
 * Destructor tells the tbox::RestartManager to remove this object from
 * the list of restart items.
 *
 *************************************************************************
 */
MultiblockGriddingAlgorithm::~MultiblockGriddingAlgorithm()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }

   if (d_internal_tagger_strategy && d_mb_tagger_strategy) {
      delete d_mb_tagger_strategy;
   }

}

/*
 *************************************************************************
 *
 * Construct coarsest level in the patch hierarchy (i.e., level 0).
 * This routine operates in two modes:
 *
 * (1) If level 0 does not exist in the hierarchy, then a new level
 *  will be created based on the domain description provided by
 *  the grid geometry associated with the hierarchy.
 * (2) If level 0 exists in the hierarchy, then a new level is made
 *  by re-applying the load balancing routines to the existing level.
 *  The pre-existing level will be removed from the hierarchy.
 *
 * In either case, the new level is placed in the hierarchy as level 0.
 * If level 0 can be refined (i.e., it will be subject to tagging cells
 * for refinement), the domain boxes are checked against the constraints
 * of regridding process.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::makeCoarsestLevel(
   const double level_time,
   const hier::MappedBoxLevel& override_mapped_box_level)
{

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, *d_mb_hierarchy, override_mapped_box_level);

   if (override_mapped_box_level.isInitialized()) {
      TBOX_ERROR("MultiblockGriddingAlgorithm::makeCoarsestLevel is not\n"
         "supporting override_mapped_box_level yet, due to incomplete coding.");
   }

   t_make_coarsest->barrierAndStart();

   const hier::MappedBoxLevelConnectorUtils dlbg_edge_utils;
   const hier::IntVector& one_vec = hier::IntVector::getOne(d_dim);


   TBOX_ASSERT(d_mb_hierarchy->getMaxNumberOfLevels() > 0);

   const int level_number = 0;
   const unsigned int nblocks = d_mb_hierarchy->getGridGeometry()->getNumberBlocks();

   if (d_proper_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_proper_nesting_complement.size());
      d_proper_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_proper_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_proper_nesting_complement[nb].resize(
               d_mb_hierarchy->getMaxNumberOfLevels(), hier::MappedBoxLevel(d_dim));
         }
      }
   }
   if (d_to_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_to_nesting_complement.size());
      d_to_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_to_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_to_nesting_complement[nb].resize(d_mb_hierarchy->getMaxNumberOfLevels());
         }
      }
   }
   if (d_from_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_from_nesting_complement.size());
      d_from_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_from_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_from_nesting_complement[nb].resize(d_mb_hierarchy->getMaxNumberOfLevels());
         }
      }
   }

   if (static_cast<unsigned int>(d_domain_complement_tree.size()) < nblocks) {
      d_domain_complement_tree.resizeArray(nblocks, hier::MappedBoxTree(d_dim));
   }

   bool level_zero_exists = (d_mb_hierarchy->getNumberOfLevels() > 0);

   hier::Connector* domain_to_domain = new hier::Connector[nblocks];
   hier::Connector* domain_to_new = new hier::Connector[nblocks];
   hier::Connector* new_to_domain = new hier::Connector[nblocks];
   hier::Connector* new_to_new = new hier::Connector[nblocks];

   tbox::Array<hier::MappedBoxLevel> new_mapped_box_level(nblocks,
      hier::MappedBoxLevel(d_dim));
   tbox::Pointer<hier::PatchLevel> old_level;

   for (unsigned int nb = 0; nb < nblocks; nb++) {

      const hier::BlockId block_id(nb); 

      hier::BoxList domain_boxes =
         d_mb_hierarchy->getGridGeometry()->getPhysicalDomain(nb);

      hier::IntVector smallest_patch(d_dim);
      hier::IntVector largest_patch(d_dim);
      hier::IntVector extend_ghosts(d_dim);
      {
      hier::IntVector smallest_box_to_refine(d_dim);
      // "false" argument: for_building_finer level = false
      getGriddingParameters(smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         *d_mb_hierarchy,
         level_number,
         false);
      }

      /*
       * If there is no level 0 in the patch hierarchy, then check
       * constraints on domain description.
       */

      if (!level_zero_exists) {

         int i = 0;
         for (hier::BoxList::Iterator itr(domain_boxes); itr; itr++, ++i) {
            hier::Box& test_box = *itr;
            for (int dir = 0; dir < d_dim.getValue(); dir++) {
               if (test_box.numberCells(dir) < smallest_patch(dir)) {
                  int error_coarsen_ratio =
                     d_tag_init_strategy->getErrorCoarsenRatio();
                  if (error_coarsen_ratio > 1) {
                     TBOX_ERROR(
                        d_object_name << ": "
                        << "\nBox " << i << ", "
                        << test_box
                        <<
                        ", from the input file violates"
                        <<
                        "the minimum patch size constraints."
                        <<
                        "\nVerify that boxes are larger than"
                        <<
                        "the maximum ghost width and/or"
                        <<
                        "\nthe specified minimum patch size."
                        <<
                        "\nNOTE: to assure the constraints are"
                        <<
                        "properly enforced during coarsening for"
                        <<
                        "\nerror computation, the minimum patch"
                        <<
                        "size is the smallest patch size multiplied"
                        <<
                        "\nby the error coarsen ratio, which is "
                        << error_coarsen_ratio
                        << " in this case."
                        << std::endl);
                  } else {
                     TBOX_ERROR(
                        d_object_name << ": "
                        << "\nBox " << i << ", "
                        << test_box
                        <<
                        ", from the input file violates"
                        <<
                        "the minimum patch size constraints."
                        <<
                        "\nVerify that boxes are larger than"
                        <<
                        "the maximum ghost width and/or"
                        <<
                        "\nthe specified minimum patch size."
                        << std::endl);
                  }
               }
            }
         }

         if (domain_boxes.boxesIntersect()) {
            TBOX_ERROR(
               d_object_name << ":  "
               << "Boxes specified for coarsest level "
               <<
               "contain intersections with each other!");
         }

         if ((d_mb_hierarchy->getMaxNumberOfLevels() > 1)
             && (!d_tag_init_strategy->coarsestLevelBoxesOK(domain_boxes))) {
            TBOX_ERROR(d_object_name << ":  "
                                     << "level gridding strategy encountered"
                                     << " a problem with the domain boxes!");
         }

      }

      /*
       * Apply load balancing algorithm to boxes describing domain
       * to determine configuration of patches on level 0.
       */

      hier::MappedBoxLevel domain_mapped_box_level(d_dim);
      {
         /* Initialize domain_mapped_box_level. */
         hier::MappedBoxSet domain_nodes;
         d_mb_hierarchy->getGridGeometry()->computePhysicalDomain(
            domain_nodes, one_vec, block_id);
         domain_mapped_box_level.swapInitialize(
            domain_nodes,
            hier::IntVector::getOne(d_dim),
            d_mb_hierarchy->getGridGeometry(),
            d_mb_hierarchy->getMPI(),
            hier::MappedBoxLevel::GLOBALIZED);
      }

      /*
       * domain_to_domain is a temporary peer connector
       * for domain used for:
       * - the required bias_attractor used by loadBalanceMappedBoxLevel.
       * This attractor is a kludge providing no useful information.
       * The load balancer interface should be modified not to
       * require it.
       * - a center for bridging the balanced mapped_box_level to itself.
       * Since domain is small, owned by just proc 0
       * and already globalized, it is fast and does not require
       * communication to findOverlaps for domain_to_domain.
       */

      domain_to_domain[nb].initialize(
         domain_mapped_box_level,
         domain_mapped_box_level,
         hier::IntVector::max(
            hier::IntVector::getOne(d_dim),
            d_mb_hierarchy->getRequiredConnectorWidth(0,0)));

      hier::OverlapConnectorAlgorithm oca;
      oca.findOverlaps(domain_to_domain[nb]);

      t_load_balance->barrierAndStart();

      t_load_balance_setup->start();

      const hier::IntVector patch_cut_factor(d_dim, 1);

      /*
       * FIXME: The code for generating the coarsest level's boxes is not
       * scalable because we involve the domain mapped_box_level in bridges and
       * Connector modifications, forcing proc 0 (which owns all domain
       * nodes) to do all of the searches.  This problem will show up in
       * the performance data if we load balance the coarsest level
       * regularly.
       *
       * A proper fix for this problem may be to use the existing coarse
       * level in Connector operations, since it is more evenly distributed.
       */

//      hier::MappedBoxLevel new_mapped_box_level(d_dim);
      new_mapped_box_level[nb] = domain_mapped_box_level;
      domain_to_new[nb].initialize(
         domain_mapped_box_level,
         new_mapped_box_level[nb],
         domain_to_domain[nb].getConnectorWidth(),
         domain_to_domain[nb].getNeighborhoodSets(),
         hier::MappedBoxLevel::DISTRIBUTED);
      domain_to_new[nb].setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      new_to_domain[nb].initialize(
         new_mapped_box_level[nb],
         domain_mapped_box_level,
         domain_to_domain[nb].getConnectorWidth(),
         domain_to_domain[nb].getNeighborhoodSets(),
         hier::MappedBoxLevel::DISTRIBUTED);
      new_to_domain[nb].setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      t_load_balance_setup->stop();
      d_load_balancer0->loadBalanceMappedBoxLevel(
         new_mapped_box_level[nb],
         new_to_domain[nb],
         domain_to_new[nb],
         d_mb_hierarchy,
         level_number,
         hier::Connector(),
         hier::Connector(),
         smallest_patch,
         largest_patch,
         domain_mapped_box_level,
         extend_ghosts,
         patch_cut_factor);
      t_load_balance->stop();

      if (d_sequentialize_patch_indices) {
         /*
          * Map to globalized sequential indices to interface with
          * current parts of the library that take patch indices.
          * This is a temporary measure needed to interface with
          * some parts of SAMRAI that requires knowing the patch
          * indices.  This step causes the patch indices to match
          * the corresponding local mapped_box indices.
          */
         sortNodes(new_mapped_box_level[nb],
            domain_to_new[nb],
            new_to_domain[nb],
            false,
            true);
      }

      if (s_check_connectors == 'y') {
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_domain[nb]) == 0);
         TBOX_ASSERT(oca.checkOverlapCorrectness(domain_to_new[nb]) == 0);
      }

      if (domain_mapped_box_level.getLocalNumberOfBoxes(0) ==
          (size_t)domain_mapped_box_level.getGlobalNumberOfBoxes()) {
         /*
          * If proc 0 owns all new boxes, it is faster find new<==>new by
          * globalizing the new boxes.
          *
          * The standard approach of bridging basically do the same,
          * but forces proc 0 to compute all the overlaps and send that
          * info to each processor one at a time.
          */
         t_find_new_to_new->barrierAndStart();
         new_to_new[nb].initialize(new_mapped_box_level[nb],
            new_mapped_box_level[nb],
            d_mb_hierarchy->getRequiredConnectorWidth(0,0));
         oca.findOverlaps(new_to_new[nb]);
         t_find_new_to_new->stop();
      } else {
         t_bridge_new_to_new->barrierAndStart();
         oca.bridge(new_to_new[nb],
            new_to_domain[nb],
            domain_to_new[nb],
            new_to_domain[nb],
            domain_to_new[nb]);
         t_bridge_new_to_new->stop();
         TBOX_ASSERT(new_to_new[nb].getConnectorWidth() ==
            d_mb_hierarchy->getRequiredConnectorWidth(0,0));
         TBOX_ASSERT(&new_to_new[nb].getBase() == &new_mapped_box_level[nb]);
         TBOX_ASSERT(&new_to_new[nb].getHead() == &new_mapped_box_level[nb]);
      }

      if (d_check_overlapping_patches != 'i') {
         checkOverlappingPatches(new_to_new[nb]);
      }
   }

   hier::MappedBoxSet full_new_mapped_boxes;
   for (unsigned int nb = 0; nb < nblocks; nb++) {

      const hier::MappedBoxSet& block_mapped_boxes =
         new_mapped_box_level[nb].getMappedBoxes();

      hier::MappedBoxSet::const_iterator ci;
      for (ci = block_mapped_boxes.begin();
           ci != block_mapped_boxes.end(); ++ci) {
         full_new_mapped_boxes.insert(*ci);
      }
   }

   hier::MappedBoxLevel full_new_mapped_box_level(
      full_new_mapped_boxes,
      hier::IntVector::getOne(d_dim),
      d_mb_hierarchy->getGridGeometry(),
      d_mb_hierarchy->getMPI());

   if (d_sequentialize_patch_indices) {
      sortNodes(full_new_mapped_box_level, true);
   }

   t_make_new->start();
   if (!level_zero_exists) {
       
      d_mb_hierarchy->makeNewPatchLevel(level_number,
         full_new_mapped_box_level);

      /*
       * Add computed Connectors to new MappedBoxLevel's cache.
       */
      tbox::Pointer<hier::PatchLevel> new_level(
         d_mb_hierarchy->getPatchLevel(level_number));
      new_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
      createConnector( 
         *new_level->getMappedBoxLevel(),
         d_mb_hierarchy->getRequiredConnectorWidth(0,0));
      if (s_print_mapped_box_level_hierarchy == 'y') {
         tbox::plog
         << "MultiblockGriddingAlgorithm::makeCoarsestLevel produced:\n";
         d_mb_hierarchy->recursivePrint(tbox::plog, "", 4);
      }

   } else {

      /*
       * Save old data before they are overwritten by the new mapped_box_level.
       */
      hier::MappedBoxLevel old_mapped_box_level = *d_mb_hierarchy->getMappedBoxLevel(0);

      old_level =
         d_mb_hierarchy->getPatchLevel(level_number);

      d_mb_hierarchy->removePatchLevel(level_number);

      d_mb_hierarchy->makeNewPatchLevel(level_number,
         full_new_mapped_box_level);

      if (s_print_mapped_box_level_hierarchy == 'y') {
         tbox::plog
         << "MultiblockGriddingAlgorithm::makeCoarsestLevel produced:\n";
         d_mb_hierarchy->recursivePrint(tbox::plog, "", 4);
      }

      tbox::Pointer<hier::PatchLevel> new_level(
         d_mb_hierarchy->getPatchLevel(level_number));
      for (unsigned int nb = 0; nb < nblocks; nb++) {

         /*
          * Compute old<==>new.  Doing it this way is not scalable,
          * but we only do this for the coarsest level.  The old approach
          * of bridging across the domain mapped_box_level is probably not scalable
          * anyway, because the domain is usually owned by just one
          * processor.
          */
         old_mapped_box_level.getPersistentOverlapConnectors().createConnector(
            *new_level->getMappedBoxLevel(),
            d_mb_hierarchy->getRequiredConnectorWidth(0, 0));
         new_level->getMappedBoxLevel()->getPersistentOverlapConnectors().createConnector(
            old_mapped_box_level,
            d_mb_hierarchy->getRequiredConnectorWidth(0, 0));


      }
   }
   t_make_new->stop();

   bool initial_time;
   tbox::Pointer<hier::PatchLevel> old_mb_level;
   if (!level_zero_exists) {
      initial_time = true;
      old_mb_level = (tbox::Pointer<hier::PatchLevel>)NULL;
   } else {
      initial_time = false;
      old_mb_level = old_level;
   }

   if (initial_time) {
      d_tag_init_strategy->initializeLevelData(d_mb_hierarchy, level_number,
         level_time,
         d_mb_hierarchy->levelCanBeRefined(level_number),
         initial_time,
         (tbox::Pointer<hier::PatchLevel>)NULL);
   } else {
      d_tag_init_strategy->initializeLevelData(d_mb_hierarchy, level_number,
         level_time,
         d_mb_hierarchy->levelCanBeRefined(level_number),
         false,
         old_mb_level);
   }

   if (level_zero_exists) {
      old_mb_level.setNull();
   }

   t_reset_hier->barrierAndStart();
   d_tag_init_strategy->resetHierarchyConfiguration(d_mb_hierarchy,
      level_number,
      level_number);
   t_reset_hier->barrierAndStop();

#ifdef MGA_RECORD_STATS
   recordStatistics( level_time );
#endif

   /*
    * Clear out Connectors to avoid dangling pointers when their
    * head/base MappedBoxLevels go out of scope first.
    */
   for (unsigned int nb = 0; nb < nblocks; nb++) {
      domain_to_domain[nb].clear();
      domain_to_new[nb].clear();
      new_to_domain[nb].clear();
   }
   delete[] domain_to_domain;
   delete[] domain_to_new;
   delete[] new_to_domain;
   delete[] new_to_new;

   d_mb_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
      *d_mb_hierarchy->getPatchLevel(0));

   t_make_coarsest->stop();
}

/*
 *************************************************************************
 *
 * Perform error estimation process on the finest hierarchy level to
 * determine if and where a new finest level is needed.  If it is
 * determined  that a new finest level should exist, it is created and
 * its patch data is allocated and initialized.  The algorithm is
 * summarized as follows:
 *
 * (1) Compute boxes whose union constitutes the region within the level
 *  in which the next finer level may reside (i.e., proper nesting).
 *
 * (2) Set tags on the level to indicate which cells should be refined.
 *
 * (3) Buffer the tags.  This prevents disturbances from moving off
 *  refined regions before the next remesh occurs.
 *
 * (4) Determine boxes for new patches that will cover the tagged cells.
 *  This step includes load balancing of the patches.
 *
 * (5) If there exist some regions to refine, construct the next finer
 *  level and insert it in the hierarchy.  Then, initialize the data
 *  on that level.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::makeFinerLevel(
   const double level_time,
   const bool initial_time,
   const int tag_buffer,
   const double regrid_start_time)
{
   TBOX_ASSERT(!(d_mb_hierarchy->getPatchLevel(
                    d_mb_hierarchy->getFinestLevelNumber()).isNull()));
   TBOX_ASSERT(tag_buffer >= 0);

   if (s_print_steps == 'y') {
      tbox::plog
      <<
      "MultiblockGriddingAlgorithm::makeFinerLevel entered with finest ln = "
      << d_mb_hierarchy->getFinestLevelNumber() << "\n";
   }

   t_make_finer->barrierAndStart();

   const unsigned int nblocks = d_mb_hierarchy->getGridGeometry()->getNumberBlocks();

   if (d_proper_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_proper_nesting_complement.size());
      d_proper_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_proper_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_proper_nesting_complement[nb].resize(
               d_mb_hierarchy->getMaxNumberOfLevels(), hier::MappedBoxLevel(d_dim));
         }
      }
   }
   if (d_to_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_to_nesting_complement.size());
      d_to_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_to_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_to_nesting_complement[nb].resize(d_mb_hierarchy->getMaxNumberOfLevels());
         }
      }
   }
   if (d_from_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_from_nesting_complement.size());
      d_from_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_from_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_from_nesting_complement[nb].resize(d_mb_hierarchy->getMaxNumberOfLevels());
         }
      }
   }

   if (static_cast<unsigned int>(d_domain_complement_tree.size()) < nblocks) {
      d_domain_complement_tree.resizeArray(nblocks, hier::MappedBoxTree(d_dim));
   }

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   const hier::MappedBoxLevelConnectorUtils dlbg_edge_utils;

   const int tag_ln = d_mb_hierarchy->getFinestLevelNumber();
   const int new_ln = tag_ln + 1;

   if (d_mb_hierarchy->levelCanBeRefined(tag_ln)) {
      t_make_finer_setup->start();

      /*
       * d_base_ln is used by private methods during regrid.
       */
      d_base_ln = tag_ln;

      if (d_mb_hierarchy->getNumberOfLevels() > d_base_ln) {

         for (unsigned int nb = 0; nb < nblocks; nb++) {
            /*
             * Compute nesting data at d_base_ln for use in constructing
             * level d_base_ln+1;
             */
            computeNestingData(d_base_ln, hier::BlockId(nb));
         }
      }

      const tbox::Pointer<hier::PatchLevel>
      tag_level = d_mb_hierarchy->getPatchLevel(tag_ln);

      tbox::Array<hier::BoxList> fine_boxes(nblocks, hier::BoxList(d_dim));

      tbox::Array<hier::MappedBoxLevel> new_mapped_box_level(
         nblocks,
         hier::MappedBoxLevel(d_dim));
      std::vector<hier::Connector> mb_new_to_new(nblocks);
      std::vector<hier::Connector> mb_new_to_tag(nblocks);
      std::vector<hier::Connector> mb_tag_to_new(nblocks);

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *   1) only user supplied refine boxes are used
       *   2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       */
      bool do_tagging = true;
      if (d_tag_init_strategy->refineUserBoxInputOnly()) do_tagging = false;
      t_make_finer_setup->stop();

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         t_make_finer_tagging->start();

         /*
          * Create communication schedule for buffer tags on this level.
          */
         t_bdry_fill_tags_create->start();
         d_bdry_sched_tags[tag_ln] =
            d_mblk_bdry_fill_tags->createSchedule(tag_level,
               d_mb_tagger_strategy);
         t_bdry_fill_tags_create->stop();

         /*
          * Initialize integer tag arrays on level to false.
          */

         tag_level->allocatePatchData(d_tag_indx);
         /*
          * FIXME: fillTagsFromMappedBoxLevel is more complex than needed here.
          * because we want to fill up the whole level, not fill selectively.
          */

         fillTagsFromMappedBoxLevel(
            d_false_tag,
            tag_level,
            d_tag_indx,
            d_mb_hierarchy->getConnector(tag_ln, tag_ln),
            true,
            zero_vector);

         /*
          * Perform pre-processing of error estimation data, if appropriate.
          */
         if (errorEstimationUsesTimeIntegration()) {
            d_tag_init_strategy->
            preprocessErrorEstimation(d_mb_hierarchy,
               tag_ln,
               level_time,
               regrid_start_time,
               initial_time);
         }

         /*
          * Determine cells needing refinement on level and set tags to true.
          * Because we are constructing a new level, not regridding the level,
          * the coarsest_sync_level argument is always false.
          */
         bool coarsest_sync_level = false;
         d_tag_init_strategy->
         tagCellsForRefinement(d_mb_hierarchy,
            tag_ln,
            level_time,
            d_tag_indx,
            initial_time,
            coarsest_sync_level,
            d_mb_hierarchy->levelCanBeRefined(tag_ln),
            regrid_start_time);

         /*
          * Check for user-tagged cells that violate proper nesting.
          * except if user specified that the violating tags be ignored.
          */
         if (d_check_nonrefined_tags != 'i') {
            for (unsigned int nb = 0; nb < nblocks; nb++) {
               checkNonrefinedTags(*tag_level, tag_ln, hier::BlockId(nb));
            }
         }

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next regrid
          * of the level occurs.
          */
         tag_level->allocatePatchData(d_buf_tag_indx);
         bufferTagsOnLevel(d_true_tag, tag_level, tag_buffer);
         tag_level->deallocatePatchData(d_buf_tag_indx);
         t_make_finer_tagging->stop();

         tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::commNull);
         tbox::Array<int> work_on_block(nblocks, 0);
         for (unsigned int nb = 0; nb < nblocks; nb++) {

            if (d_mb_hierarchy->getNumberOfLevels() > d_base_ln) {

               hier::Connector& tag_to_new = mb_tag_to_new[nb];
               hier::Connector& new_to_tag = mb_new_to_tag[nb];

               /*
                * Determine box array and processor mapping for new fine level.
                */
               findRefinementBoxes(new_mapped_box_level[nb],
                  tag_to_new,
                  new_to_tag,
                  tag_ln,
                  hier::BlockId(nb));

               work_on_block[nb] =
                  static_cast<int>(
                     new_mapped_box_level[nb].getLocalNumberOfCells());

               if (mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
                  mpi = d_mb_hierarchy->getMPI();
               }
            }
         }

         if (mpi.usingMPI()) {
            TBOX_ASSERT(mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull);
            mpi.AllReduce(work_on_block.getPointer(), nblocks, MPI_SUM);
         }

         tbox::Array<tbox::RankGroup> rank_group(nblocks,
                                                 tbox::RankGroup(mpi));

         if (nblocks > 1) {
            setRankGroups(rank_group, work_on_block, mpi);
         }

         for (unsigned int nb = 0; nb < nblocks; nb++) {

            hier::BlockId block_id(nb);
            if (d_mb_hierarchy->getNumberOfLevels() > d_base_ln) {

               hier::Connector& tag_to_new = mb_tag_to_new[nb];
               hier::Connector& new_to_tag = mb_new_to_tag[nb];
               hier::Connector& new_to_new = mb_new_to_new[nb];

               loadBalanceAndRefineBoxes(new_mapped_box_level[nb],
                  tag_to_new,
                  new_to_tag,
                  tag_ln,
                  rank_group[nb],
                  block_id);

               if (new_mapped_box_level[nb].isInitialized()) {

                  if (s_check_proper_nesting == 'y') {
                     /*
                      * Check that the new mapped_box_level nests inside the tag level.
                      *
                      * SAMRAI convention (my best understanding of it) supported
                      * (or should have been supported) by grid generation:
                      * - L0 must be equivalent to the domain.
                      * - L1 must nest in L0 by the max ghost width in L1 index space.
                      * - L(n) must nest in L(n-1) by getProperNestingBuffer(n-1)
                      */
                     hier::IntVector required_nesting(d_dim);
                     if (tag_ln > 0) {
                        required_nesting = hier::IntVector(d_dim, getProperNestingBuffer(tag_ln));
                        required_nesting *= getRatioToCoarserLevel(new_ln);
                     } else {
                        required_nesting =
                           d_mb_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_dim);
                     }
                     bool locally_nests = false;
                     const bool new_nests_in_tag =
                        dlbg_edge_utils.baseNestsInHead(
                           &locally_nests,
                           new_mapped_box_level[nb],
                           *d_mb_hierarchy->getMappedBoxLevel(tag_ln),
                           required_nesting,
                           hier::IntVector::getZero(d_dim),
                           hier::IntVector::getZero(d_dim),
                           &d_mb_hierarchy->getPeriodicDomainSearchTree(block_id));
                     if (!new_nests_in_tag) {
                        tbox::perr << "MultiblockGriddingAlgorithm::makeFinerLevel:\n"
                                   << "tag_ln=" << tag_ln << ":\n"
                                   << "new mapped_box_level does not properly nest\n"
                                   << "in tag mapped_box_level by the required nesting buffer of "
                                   << getProperNestingBuffer(tag_ln)
                                   << ".\nLocal nestingness: " << locally_nests
                                   << ".\nWriting MappedBoxLevels out to log file.\n"
                                   << "new_mapped_box_level of:\n" << new_mapped_box_level[nb].format("N->", 2)
                                   << "tag mapped_box_level of:\n" << d_mb_hierarchy->getMappedBoxLevel(tag_ln)->format("F->", 2)
                                   << "tag_to_new:\n" << tag_to_new.format("N->", 2)
                                   << "new_to_tag:\n" << new_to_tag.format("N->", 2);
                        hier::MappedBoxLevel violator(d_dim);
                        hier::Connector new_to_violator;
                        hier::MappedBoxLevelConnectorUtils edge_utils;
                        t_compute_external_parts->start();
                        edge_utils.computeExternalParts(
                           violator,
                           new_to_violator,
                           new_to_tag,
                           hier::IntVector(d_dim, -getProperNestingBuffer(tag_ln)),
                           d_mb_hierarchy->getDomainSearchTree(block_id));
                        t_compute_external_parts->stop();
                        TBOX_ERROR(
                           "Internal library error: Failed to produce proper nesting."
                           << "violator:\n" << violator.format("N->", 2)
                           << "new_to_violator:\n" << new_to_violator.format("N->", 2));
                     }
                  }

                  t_bridge_links->start();
                  t_bridge_new_to_new->start();
                  const hier::OverlapConnectorAlgorithm oca;
                  oca.bridge(new_to_new,
                     new_to_tag,
                     tag_to_new,
                     new_to_tag,
                     tag_to_new);
                  t_bridge_new_to_new->stop();
                  t_bridge_links->stop();
               }
            }
         }

         /*
          * Deallocate tag arrays and schedule -- no longer needed.
          */
         tag_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[tag_ln].setNull();

      } else {

         for (unsigned int nb = 0; nb < nblocks; nb++) {
            hier::Connector& tag_to_new = mb_tag_to_new[nb];
            hier::Connector& new_to_tag = mb_new_to_tag[nb];
            hier::Connector& new_to_new = mb_new_to_new[nb];

            /*
             * If tagging is not necessary (do_tagging = false) we simply
             * need to access the level boxes, either from a dumped file or
             * from user-supplied refine boxes, and load balance them before
             * constructing the finer level.
             */
            bool remove_old_fine_level = false;
            readLevelBoxes(fine_boxes[nb],
               new_mapped_box_level[nb],
               tag_to_new,
               new_to_tag,
               hier::BlockId(nb),
               tag_ln,
               level_time,
               remove_old_fine_level);
            if (new_mapped_box_level[nb].isInitialized()) {
               t_bridge_new_to_new->start();
               const hier::OverlapConnectorAlgorithm oca;
               oca.bridge(new_to_new,
                  new_to_tag,
                  tag_to_new,
                  new_to_tag,
                  tag_to_new);
               t_bridge_new_to_new->stop();
            }
         }
      }

      hier::MappedBoxSet full_new_mapped_boxes;
      for (unsigned int nb = 0; nb < nblocks; nb++) {

         const hier::MappedBoxSet& block_mapped_boxes =
            new_mapped_box_level[nb].getMappedBoxes();

         hier::MappedBoxSet::const_iterator ci;
         for (ci = block_mapped_boxes.begin();
              ci != block_mapped_boxes.end(); ++ci) {
            full_new_mapped_boxes.insert(*ci);
         }
      }

      hier::MappedBoxLevel full_new_mapped_box_level(
         full_new_mapped_boxes,
         getRatioToCoarserLevel(new_ln) * tag_level->getRatioToLevelZero(),
         d_mb_hierarchy->getGridGeometry(),
         d_mb_hierarchy->getMPI());

      if (d_sequentialize_patch_indices) {
         sortNodes(full_new_mapped_box_level, true);
      }

      int num_global_boxes = static_cast<int>(full_new_mapped_boxes.size());
      if (d_mb_hierarchy->getMPI().getSize() > 1) {
         d_mb_hierarchy->getMPI().AllReduce(&num_global_boxes, 1, MPI_SUM);
      }

      t_make_finer_create->start();
      /*
       * Make new finer level (new_ln == tag_ln+1),
       * if appropriate.
       */
      bool new_level_created = false;
      if (num_global_boxes) {
         new_level_created = true;
         d_mb_hierarchy->makeNewPatchLevel(new_ln,
            full_new_mapped_box_level);

         tbox::Pointer<hier::PatchLevel> new_level =
            d_mb_hierarchy->getPatchLevel(new_ln);

         const hier::Connector& new_to_new =
         new_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
            createConnector(
               *new_level->getMappedBoxLevel(),
               d_mb_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));

         if (d_check_overlapping_patches != 'i') {
            checkOverlappingPatches(new_to_new);
         }

         new_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
            createConnector(
               *tag_level->getMappedBoxLevel(),
               d_mb_hierarchy->getRequiredConnectorWidth(new_ln, tag_ln));
         tag_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
            createConnector(
               *new_level->getMappedBoxLevel(),
               d_mb_hierarchy->getRequiredConnectorWidth(tag_ln, new_ln));

         if (s_print_mapped_box_level_hierarchy == 'y') {
            tbox::plog
            << "MultiblockGriddingAlgorithm::makeFinerLevel produced:\n";
            d_mb_hierarchy->recursivePrint( tbox::plog, "", 4);
         }
      }

      t_make_finer_create->stop();

      if (new_level_created) {
         d_tag_init_strategy->initializeLevelData(d_mb_hierarchy,
            new_ln,
            level_time,
            d_mb_hierarchy->levelCanBeRefined(
               new_ln),
            initial_time,
            (tbox::Pointer<hier::PatchLevel>)NULL);

         t_reset_hier->barrierAndStart();
         d_tag_init_strategy->resetHierarchyConfiguration(d_mb_hierarchy,
            new_ln,
            new_ln);
         t_reset_hier->barrierAndStop();
      }

      d_base_ln = -1;

      if (d_mb_hierarchy->getNumberOfLevels() > new_ln) { 
         d_mb_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
            *d_mb_hierarchy->getPatchLevel(new_ln));
      }
   }  // if level cannot be refined, the routine drops through...

   t_make_finer->barrierAndStop();

}

/*
 *************************************************************************
 *
 * Regrid each level in the hierarchy which is finer than the specified
 * level.  First, we recursively compute proper nesting boxes for each
 * level that will be subject to regridding. If the regridding procedure
 * employs time integration, we perform any pre-processing necessary
 * to regrid the levels.  Then, each level finer than the specified
 * level is regridded from fine to coarse.  The recursive regridding
 * procedure is performed by the function regridFinerLevel().  Finally,
 * after the new hierarchy configuration is set, the application-
 * specific operations for resetting hierarchy-dependent infomation is
 * called.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::regridAllFinerLevels(
   const int level_number,
   const double regrid_time,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double> regrid_start_time,
   const bool level_is_coarsest_sync_level)
{
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_mb_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(!(d_mb_hierarchy->getPatchLevel(level_number).isNull()));
   TBOX_ASSERT(tag_buffer.getSize() >= level_number + 1);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif
   t_regrid_all_finer->barrierAndStart();

   if (d_mb_hierarchy->levelCanBeRefined(level_number)) {

      if (s_print_steps == 'y') {
         tbox::plog
         <<
         "MultiblockGriddingAlgorithm::regridAllFinerLevels: regridding finer than "
         << level_number << "\n";
      }

      /*
       * d_base_ln is used by private methods during regrid.
       */
      d_base_ln = level_number;

      t_process_error->start();
      /*
       * Perform pre-processing of error estimation data, if
       * appropriate.
       */
      if (errorEstimationUsesTimeIntegration()) {
         for (int ln = level_number;
              ln <= d_mb_hierarchy->getFinestLevelNumber(); ln++) {
            if (d_mb_hierarchy->levelCanBeRefined(ln)) {
               bool initial_time = false;
               double level_regrid_start_time = 0.;
               if (regrid_start_time.getSize() < ln + 1) {
                  tbox::IEEE::setNaN(level_regrid_start_time);
               } else {
                  level_regrid_start_time = regrid_start_time[ln];
               }

               d_tag_init_strategy->
               preprocessErrorEstimation(d_mb_hierarchy,
                  ln,
                  regrid_time,
                  level_regrid_start_time,
                  initial_time);
            }
         }
      }
      t_process_error->stop();

      /*
       * Recursively regrid each finer level.
       */
      const int finest_level_not_regridded = level_number;
      regridFinerLevel(
         level_number,
         regrid_time,
         finest_level_not_regridded,
         level_is_coarsest_sync_level,
         tag_buffer,
         regrid_start_time);

      /*
       * Invoke application-specific routines to reset information for those
       * levels which have been modified.
       */

      if (d_mb_hierarchy->getFinestLevelNumber() >= (level_number + 1)) {
         t_reset_hier->barrierAndStart();
         d_tag_init_strategy->
         resetHierarchyConfiguration(d_mb_hierarchy,
            level_number + 1,
            d_mb_hierarchy->
            getFinestLevelNumber());
         t_reset_hier->barrierAndStop();
      }

      d_base_ln = -1;

      if (s_print_steps == 'y') {
         tbox::plog
         <<
         "MultiblockGriddingAlgorithm::regridAllFinerLevels: regridded finer than "
         << level_number << "\n";
      }

   } //  if level cannot be refined, the routine drops through...

#ifdef MGA_RECORD_STATS
   // Verified that this does not use much time.
   recordStatistics( regrid_time );
#endif

   t_regrid_all_finer->barrierAndStop();

}

/*
 *************************************************************************
 *
 * Recursively, regrid each AMR hierarchy level finer than the specified
 * level (indexed by level_number).  The process is as follows:
 *
 * (1) Initialize tags to false on the level.
 *
 * (2) If a finer level exists, set tag to true on level for each cell
 *  that is refined.
 *
 * (3) Tag cells for refinement on level by applying application-
 *  specific error estimation routines.
 *
 * (4) If a finer level exists, invoke process recursively (i.e.,
 *  invoke step 1 on next finer level).
 *
 * (5) (Note we have popped out of recursion at this point).  Buffer
 *  true tags on current level to keep disturbances on fine grids
 *  until regridding occurs next.
 *
 * (6) Determine box configuration for new finer level, by calling
 *  findRefinementBoxes() function.
 *
 * (7) If a finer level should exist in the hierarchy, create its
 *  patches from the box description and initialize its data.  If
 *  necessary, discard old level.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::regridFinerLevel(
   const int tag_ln,
   const double regrid_time,
   const int finest_level_not_regridded,
   const bool level_is_coarsest_sync_level,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double>& regrid_start_time)
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_mb_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(!(d_mb_hierarchy->getPatchLevel(tag_ln).isNull()));
   TBOX_ASSERT(finest_level_not_regridded >= 0
      && finest_level_not_regridded <= tag_ln);
   TBOX_ASSERT(tag_buffer.getSize() >= tag_ln + 1);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *d_mb_hierarchy);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));
   const hier::IntVector neg1_vector(d_dim, -1);

   if (s_print_steps == 'y') {
      tbox::plog
      <<
      "MultiblockGriddingAlgorithm::regridFinerLevel: entered with tag_ln = "
      << tag_ln << "\n";
   }

   hier::MappedBoxLevelConnectorUtils dlbg_edge_utils;

   const unsigned int nblocks = d_mb_hierarchy->getGridGeometry()->getNumberBlocks();

   if (d_proper_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_proper_nesting_complement.size());
      d_proper_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_proper_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_proper_nesting_complement[nb].resize(
               d_mb_hierarchy->getMaxNumberOfLevels(), hier::MappedBoxLevel(d_dim));
         }
      }
   }
   if (d_to_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_to_nesting_complement.size());
      d_to_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_to_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_to_nesting_complement[nb].resize(d_mb_hierarchy->getMaxNumberOfLevels());
         }
      }
   }
   if (d_from_nesting_complement.size() < nblocks) {
      const int old_size = static_cast<int>(d_from_nesting_complement.size());
      d_from_nesting_complement.resize(nblocks);
      for (unsigned int nb = old_size; nb < nblocks; nb++) {
         if (d_from_nesting_complement[nb].size() !=
             (unsigned int)d_mb_hierarchy->getMaxNumberOfLevels()) {
            d_from_nesting_complement[nb].resize(d_mb_hierarchy->getMaxNumberOfLevels());
         }
      }
   }

   if (static_cast<unsigned int>(d_domain_complement_tree.size()) < nblocks) {
      d_domain_complement_tree.resizeArray(nblocks, hier::MappedBoxTree(d_dim));
   }

   if (d_mb_hierarchy->levelCanBeRefined(tag_ln)) {

      t_misc4->start();

      int new_ln = tag_ln + 1;

      tbox::Pointer<hier::PatchLevel>
      tag_level = d_mb_hierarchy->getPatchLevel(tag_ln);

      for (unsigned int nb = 0; nb < nblocks; nb++) {

         /*
          * Compute nesting data at tag_ln for use in constructing
          * level tag_ln+1;
          */
         if (d_mb_hierarchy->getNumberOfLevels() > tag_ln) {
            computeNestingData(tag_ln, hier::BlockId(nb));
         }
      }

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *   1) only user supplied refine boxes are used
       *   2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       *
       * The old level is generally removed when regridding, but
       * some circumstances may warrant keeping the old level.  For
       * example, if the refine region has not changed, there is no
       * need to regenerate the finer level.  The boolean
       * "remove_old_fine_level" specifies if the old level should
       * be removed.
       */
      bool do_tagging = true;
      if (d_tag_init_strategy->refineUserBoxInputOnly()) do_tagging = false;
      bool remove_old_fine_level = true;

      tbox::Array<hier::MappedBoxLevel> new_mapped_box_level(
         nblocks,
         hier::MappedBoxLevel(d_dim));
      std::vector<hier::Connector> mb_new_to_new(nblocks);
      std::vector<hier::Connector> mb_new_to_tag(nblocks);
      std::vector<hier::Connector> mb_tag_to_new(nblocks);

      t_misc4->stop();

      /*
       * tag_to_finer is [tag_ln]->[tag_ln+2].
       * finer_to_tag is [tag_ln]->[tag_ln+2].
       *
       * These are declared in this scope, computed and cached if
       * [tag_ln+2] exists.  They are used to tag the footprint of
       * [tag_ln+2] on level tag_ln.  Later on, they are used to
       * bridge for [new_ln] <-> [tag_ln+2].
       */
      hier::Connector tag_to_finer;
      hier::Connector finer_to_tag;

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         /*
          * Create communication schedule for buffer tags and set tags to
          * false.
          */

         t_misc1->start();

         t_bdry_fill_tags_create->start();
         d_bdry_sched_tags[tag_ln] =
            d_mblk_bdry_fill_tags->createSchedule(tag_level,
               d_mb_tagger_strategy);
         t_bdry_fill_tags_create->stop();

         tag_level->allocatePatchData(d_tag_indx);
         /*
          * FIXME: fillTagsFromMappedBoxLevel is more complex than needed here.
          * because we want to fill up the whole level, not fill selectively.
          */

         fillTagsFromMappedBoxLevel(
            d_false_tag,
            tag_level,
            d_tag_indx,
            d_mb_hierarchy->getConnector(tag_ln, tag_ln),
            true,
            zero_vector);

         /*
          * Set tags to true for cells that currently cover next finer level.
          * Note that this is not needed for all regridding strategies.  But
          * knowledge of currently refined cells is generally useful to avoid
          * repeated refining and coarsening of cells near boundaries of
          * refined regions where error estimates may hover around the error
          * tolerance. For example, the regridding scheme may require that
          * the error in a cell that is currently refined fall below a
          * certain tolerance  (generally different than the tolerance to
          * refine the cell in the first place) before the cell will be
          * de-refined.
          */

         if (d_mb_hierarchy->finerLevelExists(tag_ln)) {
            fillTagsFromMappedBoxLevel(
               d_true_tag,
               tag_level,
               d_tag_indx,
               d_mb_hierarchy->getConnector(tag_ln, tag_ln + 1),
               true,
               zero_vector);
         }

         /*
          * Determine cells needing refinement according to a specific
          * error estimation procedure and set to true.
          *
          * The "level_is_coarsest_sync_level" is provided as an argument
          * to this method.  Provide the additional check of whether the
          * level is not the coarsest level and that it is not a new level
          * in the d_mb_hierarchy.  If all three conditions are true, the
          * "coarsest_sync_level" argument passed into the tagCells method
          * will be true.  Otherwise, it will be false.
          */

         bool coarsest_sync_level =
            level_is_coarsest_sync_level &&
            tag_ln > 0 &&
            tag_ln <= finest_level_not_regridded;

         bool initial_time = false;
         double level_regrid_start_time = 0.;
         if (regrid_start_time.getSize() < tag_ln + 1) {
            tbox::IEEE::setNaN(level_regrid_start_time);
         } else {
            level_regrid_start_time = regrid_start_time[tag_ln];
         }

         t_misc1->barrierAndStop();

         t_tag_cells_for_refinement->barrierAndStart();
         d_tag_init_strategy->
         tagCellsForRefinement(d_mb_hierarchy,
            tag_ln,
            regrid_time,
            d_tag_indx,
            initial_time,
            coarsest_sync_level,
            d_mb_hierarchy->levelCanBeRefined(tag_ln),
            level_regrid_start_time);
         t_tag_cells_for_refinement->barrierAndStop();

         /*
          * Check for user-tagged cells that violate proper nesting.
          * except if user specified that the violating tags be ignored.
          */
         if (d_check_nonrefined_tags != 'i') {
            for (unsigned int nb = 0; nb < nblocks; nb++) {
               checkNonrefinedTags(*tag_level, tag_ln, hier::BlockId(nb));
            }
         }

         /*
          * Perform regridding recursively on finer levels, if appropriate.
          */
         if (d_mb_hierarchy->finerLevelExists(tag_ln)
             && d_mb_hierarchy->levelCanBeRefined(new_ln)) {
            regridFinerLevel(
               new_ln,
               regrid_time,
               finest_level_not_regridded,
               false,
               tag_buffer,
               regrid_start_time);
         }

         t_misc2->barrierAndStart();

         if (d_mb_hierarchy->finerLevelExists(new_ln)) {

            t_second_finer_tagging->start();
            const hier::IntVector fill_box_growth =
               getRatioToCoarserLevel(tag_ln + 2) * getProperNestingBuffer(
                  tag_ln + 1);

            tbox::Pointer<hier::PatchLevel> second_fine_level =
               d_mb_hierarchy->getPatchLevel(new_ln + 1);
            /*
             * We require [tag_ln+1]+tag_buffer to nest
             * [tag_ln+2]+overflow.
             * If the fine Connector width >= tag_buffer+overflow,
             * the [tag_ln] -> [tag_ln+2] should see all the edges
             * it needs to see for second_fine nesting.  (Note that
             * we currently allow zero overflow.)
             *
             * FIXME: potential problem.  tag_to_finer and finer_to_tag
             * may miss some edges to periodic images (I forgot
             * why--need to verify this statement).  This is not a
             * problem for tagging, but these Connectors are also used
             * to bridge for new_to_finer and finer_to_new later on.
             */

            if (d_mb_hierarchy->finerLevelExists(new_ln)) {

               const hier::Connector& tag_to_old = d_mb_hierarchy->getConnector(tag_ln,
                     new_ln);
               const hier::Connector& old_to_finer = d_mb_hierarchy->getConnector(new_ln,
                     new_ln + 1);
               const hier::Connector& finer_to_old = d_mb_hierarchy->getConnector(new_ln + 1,
                     new_ln);
               const hier::Connector& old_to_tag = d_mb_hierarchy->getConnector(new_ln,
                     tag_ln);
               const hier::OverlapConnectorAlgorithm oca;
               oca.bridge(tag_to_finer,
                  finer_to_tag,
                  tag_to_old,
                  old_to_finer,
                  finer_to_old,
                  old_to_tag);

#ifdef DEBUG_CHECK_ASSERTIONS
               oca.assertOverlapCorrectness(tag_to_finer, false, true, true);
               oca.assertOverlapCorrectness(finer_to_tag, false, true, true);
               TBOX_ASSERT(
                  tag_to_finer.getConnectorWidth()
                  * getRatioToCoarserLevel(
                     tag_ln + 1) * getRatioToCoarserLevel(
                     tag_ln + 2) >= fill_box_growth);
#endif

               fillTagsFromMappedBoxLevel(
                  d_true_tag,
                  tag_level,
                  d_tag_indx,
                  tag_to_finer,
                  true,
                  fill_box_growth);
            }
            t_second_finer_tagging->barrierAndStop();
         }

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next
          * regrid of the level occurs.
          */
         tag_level->allocatePatchData(d_buf_tag_indx);
         bufferTagsOnLevel(d_true_tag, tag_level, tag_buffer[tag_ln]);
         tag_level->deallocatePatchData(d_buf_tag_indx);

         t_misc2->stop();

         /*
          * Determine box array containing cells on level with a true tag
          * value.  The box array must be contained in array of proper
          * nesting boxes.
          */

         tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::commNull);
         tbox::Array<int> work_on_block(nblocks, 0);
         for (unsigned int nb = 0; nb < nblocks; nb++) {

            hier::Connector& tag_to_new = mb_tag_to_new[nb];
            hier::Connector& new_to_tag = mb_new_to_tag[nb];

            if (d_mb_hierarchy->getFinestLevelNumber() >= tag_ln) {

               findRefinementBoxes(new_mapped_box_level[nb],
                  tag_to_new,
                  new_to_tag,
                  tag_ln,
                  hier::BlockId(nb));

               if (s_print_steps == 'y') {
                  if (new_mapped_box_level[nb].isInitialized()) {
                     tbox::plog
                     <<
                     "MultiblockGriddingAlgorithm::regridFinerLevel got inititalized new_mapped_box_level\n";
                  } else {
                     tbox::plog
                     <<
                     "MultiblockGriddingAlgorithm::regridFinerLevel got un-inititalized new_mapped_box_level\n";
                  }
               }

               if (new_mapped_box_level[nb].isInitialized()) {
                  work_on_block[nb] =
                     static_cast<int>(
                        new_mapped_box_level[nb].getLocalNumberOfCells());

                  mpi = new_mapped_box_level[nb].getMPI();
               }
            }
         }

         if (mpi.usingMPI()) {
            TBOX_ASSERT(mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull);
            mpi.AllReduce(work_on_block.getPointer(), nblocks, MPI_SUM);
         }

         tbox::Array<tbox::RankGroup> rank_group(nblocks,
                                                 tbox::RankGroup(mpi));

         if (nblocks > 1) {
            setRankGroups(rank_group, work_on_block, mpi);
         }

         for (unsigned int nb = 0; nb < nblocks; nb++) {

            hier::Connector& tag_to_new = mb_tag_to_new[nb];
            hier::Connector& new_to_tag = mb_new_to_tag[nb];
            hier::Connector& new_to_new = mb_new_to_new[nb];

            if (d_mb_hierarchy->getFinestLevelNumber() >= tag_ln) {

               loadBalanceAndRefineBoxes(new_mapped_box_level[nb],
                  tag_to_new,
                  new_to_tag,
                  tag_ln,
                  rank_group[nb],
                  hier::BlockId(nb));

               if (new_mapped_box_level[nb].isInitialized()) {
                  if (s_print_steps == 'y') {
                     tbox::plog
                     <<
                     "MultiblockGriddingAlgorithm::findRefinementBoxes: bridge links\n";
                  }
                  t_bridge_links->start();
                  if (tag_ln == 0) {
                     t_bridge_new_to_new->start();
                     const hier::OverlapConnectorAlgorithm oca;
                     oca.bridge(new_to_new,
                        new_to_tag,
                        tag_to_new,
                        new_to_tag,
                        tag_to_new);
                     const hier::IntVector& gcw_new_to_new =
                        d_mb_hierarchy->getRequiredConnectorWidth(new_ln, new_ln);
                     oca.shrinkConnectorWidth(new_to_new, gcw_new_to_new);
                     t_bridge_new_to_new->stop();
                  } else {
                     /*
                      * I think we can bypass bridging to the coarser level
                      * first because the new mapped_box_level (after overflow restriction)
                      * and growth by the required gcw(new->new) still nests
                      * in tag^gcw(tag->new).  This is guaranteed if
                      * gcw(overflow) + gcw(new->new) <= gcw(tag->new).
                      */
                     const hier::IntVector& gcw_new_to_new =
                        d_mb_hierarchy->getRequiredConnectorWidth(new_ln, new_ln);
                     t_bridge_new_to_new->start();
                     const hier::OverlapConnectorAlgorithm oca;
                     oca.bridgeWithNesting(
                        new_to_new,
                        new_to_new,
                        new_to_tag,
                        tag_to_new,
                        new_to_tag,
                        tag_to_new,
                        zero_vector,
                        zero_vector,
                        gcw_new_to_new);
                     t_bridge_new_to_new->stop();
                  }
                  t_bridge_links->stop();
               }
            }
         }

         /*
          * Deallocate tag arrays and schedule; no longer needed on current
          * level.
          */

         tag_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[tag_ln].setNull();

      } else {
         TBOX_ERROR("Disabled because we don't support reading for new boxes.");
      }

      /*
       * Make new finer level (new_ln) if necessary, or remove
       * next finer level if it is no longer needed.
       */
      for (unsigned int nb = 0; nb < nblocks; nb++) {
         if (new_mapped_box_level[nb].isInitialized()) {
            hier::Connector& new_to_new = mb_new_to_new[nb];

            if (d_check_overlapping_patches != 'i') {
               checkOverlappingPatches(new_to_new);
            }
         }
      }

      hier::MappedBoxSet full_new_mapped_boxes;
      for (unsigned int nb = 0; nb < nblocks; nb++) {

         const hier::MappedBoxSet& block_mapped_boxes =
            new_mapped_box_level[nb].getMappedBoxes();

         hier::MappedBoxSet::const_iterator ci;
         for (ci = block_mapped_boxes.begin();
              ci != block_mapped_boxes.end(); ++ci) {
            full_new_mapped_boxes.insert(*ci);
         }
      }

      hier::MappedBoxLevel full_new_mapped_box_level(
         full_new_mapped_boxes,
         getRatioToCoarserLevel(new_ln) * tag_level->getRatioToLevelZero(),
         d_mb_hierarchy->getGridGeometry(),
         d_mb_hierarchy->getMPI());

      if (d_sequentialize_patch_indices) {
         sortNodes(full_new_mapped_box_level, true);
      }

      int num_global_boxes = static_cast<int>(full_new_mapped_boxes.size());
      if (d_mb_hierarchy->getMPI().getSize() > 1) {
         d_mb_hierarchy->getMPI().AllReduce(&num_global_boxes, 1, MPI_SUM);
      }

      bool initialize_level = false;
      tbox::Pointer<hier::PatchLevel> old_fine_mblk_level(NULL);
      if (d_mb_hierarchy->getNumberOfLevels() > new_ln) {
         old_fine_mblk_level = d_mb_hierarchy->getPatchLevel(new_ln);
      }

      tbox::Pointer<hier::Connector> old_to_new;
      tbox::Pointer<hier::Connector> new_to_old;

      if (num_global_boxes) {
         initialize_level = true;

         t_regrid_finer_create->start();
         /*
          * Either remove pre-existing fine level from hierarchy and make
          * a new level, or just make a new fine level for hierarchy.
          */

         tbox::Pointer<hier::PatchLevel> old_fine_level;

         tbox::ConstPointer<hier::MappedBoxLevel> old_mapped_box_level;
         const hier::Connector* old_to_tag = NULL;
         const hier::Connector* tag_to_old = NULL;

         if (d_mb_hierarchy->finerLevelExists(tag_ln)) {

            // Save a reference to old mapped_box_level before d_mb_hierarchy dumps it.
            old_mapped_box_level = d_mb_hierarchy->getMappedBoxLevel(new_ln);
            old_to_tag = &d_mb_hierarchy->getConnector(new_ln, tag_ln);
            tag_to_old = &d_mb_hierarchy->getConnector(tag_ln, new_ln);

            old_fine_level = old_fine_mblk_level;
            d_mb_hierarchy->removePatchLevel(new_ln);
         }

         if (d_mb_hierarchy->levelExists(new_ln + 1)) {
            if (s_check_proper_nesting == 'y') {
               for (unsigned int nb = 0; nb < nblocks; nb++) {
                  hier::BlockId block_id(nb);
                  /*
                   * Check that the new mapped_box_level nests the finer
                   * level (new_ln+1).  The mechanism for causing this
                   * nest to happen is that the footprint of
                   * mapped_box_level (new_ln+1) was used to tag for the
                   * new mapped_box_level.
                   */
                  const hier::IntVector required_nesting =
                     getRatioToCoarserLevel(new_ln
                        + 1) * getProperNestingBuffer(new_ln);
                  bool locally_nests = false;
                  const bool finer_nests_in_new =
                     dlbg_edge_utils.baseNestsInHead(
                        &locally_nests,
                        *d_mb_hierarchy->getMappedBoxLevel(new_ln + 1),
                        new_mapped_box_level[nb],
                        required_nesting,
                        hier::IntVector::getZero(d_dim),
                        hier::IntVector::getZero(d_dim),
                        &d_mb_hierarchy->getPeriodicDomainSearchTree(block_id));
                  if (!finer_nests_in_new) {
                     tbox::perr
                     << "MultiblockGriddingAlgorithm: new mapped_box_level\n"
                     << "at ln=" << new_ln
                     << " does not properly nest\n"
                     << "existing finer mapped_box_level at ln="
                     << new_ln + 1
                     << " by the required nesting buffer of "
                     << required_nesting
                     << ".\nLocal nestingness: " << locally_nests
                     << ".\nWriting MappedBoxLevels out to log file.\n"
                     << "new_mapped_box_level of\n" << new_mapped_box_level[nb].format("N->", 2)
                     << "finer mapped_box_level of\n" << d_mb_hierarchy->getMappedBoxLevel(new_ln + 1)->format("F->", 2);
                     hier::MappedBoxLevel external(d_dim);
                     hier::Connector finer_to_new, finer_to_external;
                     finer_to_new.initialize(
                        *d_mb_hierarchy->getMappedBoxLevel(new_ln + 1),
                        new_mapped_box_level[nb],
                        required_nesting);
                     const hier::OverlapConnectorAlgorithm oca;
                     oca.findOverlaps(finer_to_new);
                     tbox::plog << "Finer to new:\n" << finer_to_new.format("FN->", 3);
                     dlbg_edge_utils.computeExternalParts(
                        external,
                        finer_to_external,
                        finer_to_new,
                        -required_nesting,
                        d_mb_hierarchy->getDomainSearchTree(block_id));
                     tbox::plog << "External parts:\n" << finer_to_external.format("FE->", 3);
                     TBOX_ERROR(
                        "Internal library error: Failed to produce proper nesting.");
                  }
               }
            }
         }

         d_mb_hierarchy->makeNewPatchLevel(new_ln,
            full_new_mapped_box_level);

         /*
          * Cache Connectors for new level.
          */
         const tbox::Pointer<hier::PatchLevel>
         tag_level = d_mb_hierarchy->getPatchLevel(tag_ln);
         hier::Connector tag_to_new;
         hier::Connector new_to_tag;
         hier::Connector new_to_new;

         tbox::Pointer<hier::PatchLevel> new_level =
            d_mb_hierarchy->getPatchLevel(new_ln);

         new_to_new =
         new_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
         createConnector(
            *new_level->getMappedBoxLevel(),
            d_mb_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));
      
         new_to_tag =
         new_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
         createConnector(
            *tag_level->getMappedBoxLevel(),
            d_mb_hierarchy->getRequiredConnectorWidth(new_ln, tag_ln));

         tag_to_new =
         tag_level->getMappedBoxLevel()->getPersistentOverlapConnectors().
         createConnector(
            *new_level->getMappedBoxLevel(),
            d_mb_hierarchy->getRequiredConnectorWidth(tag_ln, new_ln));

         if (d_mb_hierarchy->levelExists(new_ln + 1)) {

            /*
             * Connect the new level to the finer level.
             */
            hier::Connector new_to_finer;
            hier::Connector finer_to_new;
            t_bridge_new_to_finer->start();

            /*
             * FIXME: We should use gcw_limit in this bridge to make sure
             * new_to_finer and finer_to_new do not have more edges than
             * needed (and resulting in degraded performance.
             */
            const hier::OverlapConnectorAlgorithm oca;
            oca.bridgeWithNesting(
               new_to_finer,
               finer_to_new,
               new_to_tag,
               tag_to_finer,
               finer_to_tag,
               tag_to_new,
               zero_vector,
               neg1_vector,
               d_mb_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1));
               t_bridge_new_to_finer->stop();

#ifdef DEBUG_CHECK_ASSERTIONS
            oca.assertOverlapCorrectness(new_to_finer, false, true, true);
            oca.assertOverlapCorrectness(finer_to_new, false, true, true);
#endif
            tbox::Pointer<hier::PatchLevel> finer_level =
               d_mb_hierarchy->getPatchLevel(new_ln + 1);
            new_level->getMappedBoxLevel()->getPersistentOverlapConnectors()
            .createConnector(
               *finer_level->getMappedBoxLevel(),
               new_to_finer.getConnectorWidth(),
               new_to_finer.getNeighborhoodSets());
            finer_level->getMappedBoxLevel()->getPersistentOverlapConnectors()
            .createConnector(
               *new_level->getMappedBoxLevel(),
               finer_to_new.getConnectorWidth(),
               finer_to_new.getNeighborhoodSets());
         }

         old_to_new = new hier::Connector;
         new_to_old = new hier::Connector;
         if (!old_mapped_box_level.isNull()) {
            t_bridge_new_to_old->start();
            const hier::OverlapConnectorAlgorithm oca;
            oca.bridgeWithNesting(
               *(old_to_new),
               *(new_to_old),
               *old_to_tag,
               d_mb_hierarchy->getConnector(tag_ln, tag_ln + 1),
               d_mb_hierarchy->getConnector(tag_ln + 1, tag_ln),
               *tag_to_old,
               zero_vector,
               zero_vector,
               d_mb_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));
            t_bridge_new_to_old->stop();

            new_level->getMappedBoxLevel()->getPersistentOverlapConnectors()
            .createConnector(
               *old_fine_level->getMappedBoxLevel(),
               new_to_old->getConnectorWidth(),
               new_to_old->getNeighborhoodSets());
            old_fine_level->getMappedBoxLevel()->
            getPersistentOverlapConnectors().createConnector(
               *new_level->getMappedBoxLevel(),
               old_to_new->getConnectorWidth(),
               old_to_new->getNeighborhoodSets());
         }

         t_regrid_finer_create->barrierAndStop();

      } else {
         /*
          * Here if a block had patches at level new_ln before this
          * regrid but will not after this regrid, we remove the old
          * level from the PatchHierarchy for this block.
          */
         if (d_mb_hierarchy->finerLevelExists(tag_ln)) {
            d_mb_hierarchy->removePatchLevel(new_ln);
         }
      }

      if (d_mb_hierarchy->getNumberOfLevels() > new_ln) {
         d_mb_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
            *d_mb_hierarchy->getPatchLevel(new_ln));
      }

      if (initialize_level) {

         // "false" argument": const bool initial_time = false;
         d_tag_init_strategy->initializeLevelData(d_mb_hierarchy,
            new_ln,
            regrid_time,
            d_mb_hierarchy->levelCanBeRefined(
               new_ln),
            false,
            old_fine_mblk_level);

         /*
          * Destroy old patch level, if such a level existed prior to regrid.
          */
         if (!old_fine_mblk_level.isNull()) {
            tbox::Pointer<hier::PatchLevel> old_fine_level =
               old_fine_mblk_level;
            old_fine_level.setNull();

            old_fine_mblk_level.setNull();
         }
      }

      if (!initialize_level) {

         /*
          * If there are no boxes for the new fine level, remove the
          * pre-existing fine level if it existed.
          */

         if (d_mb_hierarchy->finerLevelExists(tag_ln)
             && remove_old_fine_level) {
            d_mb_hierarchy->removePatchLevel(new_ln);
         }

      } // if we are not re-gridding the level


   } //  if level cannot be refined, the routine drops through...

}

/*
 *************************************************************************
 *                                                                       *
 *                                                                       *
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::recordStatistics(
   double current_time)
{
#ifdef MGA_RECORD_STATS
// MGA_RECORD_STATS is defined in MultiblockGriddingAlgorithm.h
/*
 * For statistics, record number of cells and patches on new level.
 */
   for (int ln = 0; ln < d_mb_hierarchy->getNumberOfLevels(); ++ln) {

      int level_gridcells = 0;
      int level_local_patches = 0;

      tbox::Pointer<hier::PatchLevel> mb_patch_level =
         d_mb_hierarchy->getPatchLevel(ln);

      level_gridcells += mb_patch_level->getLocalNumberOfCells();
      level_local_patches += mb_patch_level->getLocalNumberOfPatches();

      d_boxes_stat[ln]->recordProcStat(double(level_local_patches));
      d_cells_stat[ln]->recordProcStat(double(level_gridcells));
      d_timestamp_stat[ln]->recordProcStat(double(current_time));
   }
#endif
}

/*
 *************************************************************************
 *                                                                       *
 *                                                                       *
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::printStatistics(
   std::ostream& s) const
{
#ifdef MGA_RECORD_STATS
   tbox::Pointer<hier::PatchHierarchy> hierarchy0(
      d_mb_hierarchy);
   tbox::Pointer<hier::PatchLevel> level0(
      hierarchy0->getPatchLevel(0));

   const tbox::SAMRAI_MPI& mpi(level0->getMappedBoxLevel()->getMPI());
   /*
    * Output statistics.
    */
   // Collect statistic on mesh size.
   tbox::Statistician* statn = tbox::Statistician::getStatistician();

   statn->finalize();
   // statn->printLocalStatData(s);
   if (mpi.getRank() == 0) {
      // statn->printAllGlobalStatData(s);
      for (int ln = 0; ln < d_mb_hierarchy->getMaxNumberOfLevels(); ++ln) {
         tbox::Statistic& cstat = *d_cells_stat[ln];
         tbox::Statistic& bstat = *d_boxes_stat[ln];
         tbox::Statistic& tstat = *d_timestamp_stat[ln];
         s << "statistic " << cstat.getName() << ":" << std::endl;
         if (0) {
            s << "Global: \n";
            statn->printGlobalProcStatDataFormatted(cstat.getInstanceId(), s);
         }
         s
         <<
         " Seq#  SimTime       C-Sum      C-Avg      C-Min   ->    C-Max  C-NormDiff  B-Sum B-Avg B-Min -> B-Max B-NormDiff C/B-Avg\n";
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
         for (int sn = 0; sn < cstat.getStatSequenceLength(); ++sn) {
            double csum = statn->getGlobalProcStatSum(cstat.getInstanceId(), sn);
            double cmax = statn->getGlobalProcStatMax(cstat.getInstanceId(), sn);
            double cmin = statn->getGlobalProcStatMin(cstat.getInstanceId(), sn);
            double cdiffnorm = cmax != 0 ? 1.0 - cmin / cmax : 0;
            double bsum = statn->getGlobalProcStatSum(bstat.getInstanceId(), sn);
            double bmax = statn->getGlobalProcStatMax(bstat.getInstanceId(), sn);
            double bmin = statn->getGlobalProcStatMin(bstat.getInstanceId(), sn);
            double bdiffnorm = bmax != 0 ? 1.0 - bmin / bmax : 0;
            double stime = statn->getGlobalProcStatMin(
                  tstat.getInstanceId(), sn);
            s << std::setw(3) << sn << "  "
              << std::scientific << std::setprecision(6) << std::setw(12)
              << stime
              << " "
              << std::fixed << std::setprecision(0)
              << std::setw(10) << csum << " "
              << std::setw(10) << csum / mpi.getSize() << " "
              << std::setw(10) << cmin << " -> "
              << std::setw(10) << cmax
              << "  " << std::setw(4) << std::setprecision(4) << cdiffnorm
              << "  "
              << std::fixed << std::setprecision(0)
              << std::setw(6) << bsum << " "
              << std::setw(5) << bsum / mpi.getSize() << " "
              << std::setw(5) << bmin << "  ->"
              << std::setw(5) << bmax
              << "   " << std::setw(4) << std::setprecision(4) << bdiffnorm
              << std::setw(10) << std::setprecision(0)
              << (bsum != 0 ? csum / bsum : 0)
              << std::endl;
         }
      }
   }
#endif
}

/*
 *************************************************************************
 * All tags reside in the given level.  But due to nesting restrictions,
 * not all cells in the level are allowed to be refined.  We look for the
 * portions of the level that would violate nesting if refined.  Any tags
 * there are nonnesting tags.
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::checkNonrefinedTags(
   const hier::PatchLevel& level,
   int tag_ln,
   const hier::BlockId& block_id) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, level, *d_mb_hierarchy);

   const tbox::Pointer<hier::PatchHierarchy> patch_hierarchy(d_mb_hierarchy);
   const tbox::Pointer<hier::PatchLevel> tag_level = patch_hierarchy->getPatchLevel(tag_ln);
   const hier::MappedBoxLevel& tag_mapped_box_level = *tag_level->getMappedBoxLevel();
   hier::MappedBoxLevel violator(d_dim);
   hier::Connector tag_to_violator;
   const hier::Connector& tag_mapped_box_level_to_self = patch_hierarchy->getConnector(tag_ln,
         tag_ln);
   computeNestingViolator(tag_mapped_box_level,
      violator,
      tag_to_violator,
      tag_mapped_box_level_to_self,
      tag_mapped_box_level_to_self,
      tag_ln + 1,
      block_id);

   /*
    * Check for user-tagged cells in the violating parts of the tag level.
    */
   math::PatchCellDataBasicOps<int> dataop;
   math::PatchCellDataOpsInteger dataopi;
   int maxval = 0;
   const hier::NeighborhoodSet& tag_eto_violator = tag_to_violator.getNeighborhoodSets();
   for (hier::NeighborhoodSet::const_iterator ei = tag_eto_violator.begin();
        ei != tag_eto_violator.end(); ++ei) {
      const hier::MappedBoxId &mapped_box_id = ei->first;
      tbox::Pointer<hier::Patch> patch = level.getPatch(mapped_box_id);
      tbox::Pointer<pdat::CellData<int> > tag_data =
         patch->getPatchData(d_tag_indx);
      const NeighborSet& nabrs = (*ei).second;
      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end(); ++na) {
         const hier::MappedBox& vio_mapped_box = *na;
         maxval = dataop.max(tag_data, vio_mapped_box.getBox());
         if (maxval > 0) {
            // std::cout << "violator found in vio_mapped_box " << vio_mapped_box << std::endl;
            // dataopi.printData( tag_data, tag_data->getBox(), std::cout );
            break;
         }
      }
      if (maxval > 0) {
         // std::cout << "Violator found" << std::endl;
         break;
      }
   }
   const tbox::SAMRAI_MPI mpi(violator.getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&maxval, 1, MPI_MAX);
   }

   if (maxval > 0) {
      if (d_check_nonrefined_tags == 'w') {
         TBOX_WARNING("User code has tagged cells in\n"
            << "violation of nesting requirements.\n"
            << "Violating tags will be discarded.\n"
            << "See MultiblockGriddingAlgorithm::checkNonrefinedTags()\n");
      } else if (d_check_nonrefined_tags == 'e') {
         TBOX_ERROR("User code has tagged cells in\n"
            << "violation of nesting requirements.\n"
            << "See MultiblockGriddingAlgorithm::checkNonrefinedTags()\n");
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::checkOverlappingPatches(
   const hier::Connector& mapped_box_level_to_self) const
{
   int has_overlap = 0;
   const hier::MappedBoxLevel& mapped_box_level = mapped_box_level_to_self.getBase();
   const hier::NeighborhoodSet& edges = mapped_box_level_to_self.getNeighborhoodSets();
   for (hier::NeighborhoodSet::const_iterator ei = edges.begin();
        ei != edges.end(); ++ei) {
      const hier::MappedBox& mapped_box = *mapped_box_level.getMappedBoxStrict(ei->first);
      const NeighborSet& nabrs = ei->second;
      for (NeighborSet::const_iterator na = nabrs.begin();
           na != nabrs.end();
           ++na) {
         const hier::MappedBox& nabr = *na;
         if (nabr != mapped_box &&
             nabr.getBlockId() == mapped_box.getBlockId()) {
            if (nabr.getBox().intersects(mapped_box.getBox())) {
               has_overlap = 1;
               break;
            }
         }
      }
      if (has_overlap) {
         break;
      }
   }
   const tbox::SAMRAI_MPI mpi(mapped_box_level.getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&has_overlap, 1, MPI_MAX);
   }

   if (has_overlap > 0) {
      if (d_check_overlapping_patches == 'w') {
         TBOX_WARNING("PatchLevel has patches which overlap in index space\n"
            << "See MultiblockGriddingAlgorithm::checkOverlappingPatches()\n");
      } else if (d_check_overlapping_patches == 'e') {
         TBOX_ERROR("PatchLevel has patches which overlap in index space\n"
            << "See MultiblockGriddingAlgorithm::checkOverlappingPatches()\n");
      }
   }
}

/*
 *************************************************************************
 *
 * For cases where tagging is not performed read the new level boxes
 * either from user input or from stored level boxes.
 *
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::readLevelBoxes(
   hier::BoxList& new_level_boxes,
   hier::MappedBoxLevel& new_mapped_box_level,
   hier::Connector& coarser_to_new,
   hier::Connector& new_to_coarser,
   const hier::BlockId& block_id,
   const int level_number,
   const double regrid_time,
   bool& remove_old_fine_level)
{
   (void)new_level_boxes;
   const tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);

   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, *hierarchy, new_mapped_box_level);

   const hier::MappedBoxLevel& coarser_mapped_box_level = *hierarchy->getMappedBoxLevel(level_number);

   int fine_level_number = level_number + 1;
   hier::BoxList boxes_to_refine(d_dim);

   /*
    * Access the user supplied refine boxes.  The
    * "new_level_has_new_boxes" boolean specifies whether the
    * level boxes have changed from the last time
    * getUserSuppliedRefineBoxes() was called.  If they have changed,
    * it returns true.  If they are unchanged, it returns false.
    */
   bool new_level_has_new_boxes = true;
   if (d_tag_init_strategy->refineUserBoxInputOnly()) {

      new_level_has_new_boxes = d_tag_init_strategy->
         getUserSuppliedRefineBoxes(boxes_to_refine,
            level_number,
            regrid_time);
      if (boxes_to_refine.size() && block_id.getBlockValue() > 0) {
         new_level_has_new_boxes = true;
      }
   }

   /*
    * If "new_level_has_new_boxes" is false we wish to keep the
    * existing fine level intact.  Avoid further work by setting
    * the parameter "compute_load_balanced_level_boxes" to false
    * and indicate that we want to avoid removing the old fine level
    * by setting "remove_old_fine_level" to false.
    */
   bool compute_load_balanced_level_boxes = true;
   if (!new_level_has_new_boxes) {
      compute_load_balanced_level_boxes = false;
      remove_old_fine_level = false;
   }

   /*
    * If we are using the nonuniform load balance option, we
    * still need to redo the load balance and construct a new level,
    * even if the level boxes have not changed.
    */

   if (d_load_balancer->getLoadBalanceDependsOnPatchData(fine_level_number)
       && boxes_to_refine.getNumberOfBoxes() > 0) {
      compute_load_balanced_level_boxes = true;
      remove_old_fine_level = true;
   }

   /*
    * If the boxes_to_refine are empty, this implies that no
    * refinement is desired so a new finer level will NOT be
    * constructed.  In this case, avoid load balance steps and
    * specify that we want to remove the old fine level.
    */
   if (boxes_to_refine.getNumberOfBoxes() == 0) {
      compute_load_balanced_level_boxes = false;
      remove_old_fine_level = true;
   }

   if (compute_load_balanced_level_boxes) {

      // hier::IntVector patch_cut_factor(getRatioToCoarserLevel(fine_level_number));
      hier::IntVector patch_cut_factor(d_dim, 1);

      hier::IntVector smallest_patch(d_dim);
      hier::IntVector largest_patch(d_dim);
      hier::IntVector extend_ghosts(d_dim);
      {
      hier::IntVector smallest_box_to_refine(d_dim);
      // "false" argument: for_building_finer level = false
      getGriddingParameters(smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         *hierarchy,
         fine_level_number,
         false);
      }

      hier::MappedBoxLevel unbalanced_mapped_box_level(d_dim);
      unbalanced_mapped_box_level.initialize(
         coarser_mapped_box_level.getRefinementRatio(),
         coarser_mapped_box_level.getGridGeometry(),
         coarser_mapped_box_level.getMPI(),
         hier::MappedBoxLevel::GLOBALIZED);
      hier::LocalId i(0);
      for (hier::BoxList::Iterator itr(boxes_to_refine); itr; itr++, ++i) {
         hier::MappedBox unbalanced_mapped_box(*itr, i, 0, block_id);
         unbalanced_mapped_box_level.addMappedBox(unbalanced_mapped_box);
      }

      // Get the max gcw required for peer and cross Connectors.
      hier::IntVector gcw =
         hierarchy->getRequiredConnectorWidth(level_number + 1, level_number + 1);
      gcw =
         hier::IntVector::ceiling(gcw, getRatioToCoarserLevel(fine_level_number));
      gcw.max(hierarchy->getRequiredConnectorWidth(level_number, level_number + 1));
      gcw.max(hierarchy->getRequiredConnectorWidth(level_number, level_number));

      new_mapped_box_level = unbalanced_mapped_box_level;
      coarser_to_new.initialize(
         coarser_mapped_box_level,
         new_mapped_box_level,
         gcw);
      new_to_coarser.initialize(
         new_mapped_box_level,
         coarser_mapped_box_level,
         gcw);
      const hier::OverlapConnectorAlgorithm oca;
      oca.findOverlaps(coarser_to_new);
      oca.findOverlaps(new_to_coarser);

      t_load_balance->start();
      d_load_balancer0->loadBalanceMappedBoxLevel(
         new_mapped_box_level,
         new_to_coarser,
         coarser_to_new,
         hierarchy,
         level_number,
         hier::Connector(),
         hier::Connector(),
         smallest_patch,
         largest_patch,
         d_singleblock_domain_mapped_box_level[block_id.getBlockValue()],
         extend_ghosts,
         patch_cut_factor);
      t_load_balance->stop();

      const hier::IntVector& ratio =
         getRatioToCoarserLevel(fine_level_number);
      refineNewMappedBoxLevel(new_mapped_box_level,
         coarser_to_new,
         new_to_coarser,
         ratio);
      if (d_sequentialize_patch_indices) {
         sortNodes(new_mapped_box_level,
            coarser_to_new,
            new_to_coarser,
            false,
            true);
      }
      const hier::Connector& coarser_to_coarser =
         hierarchy->getConnector(level_number,
            level_number);
      const hier::MappedBoxLevelConnectorUtils dlbg_edge_utils;
      dlbg_edge_utils.addPeriodicImagesAndRelationships(
         new_mapped_box_level,
         new_to_coarser,
         coarser_to_new,
         hierarchy->getDomainSearchTree(block_id),
         coarser_to_coarser);
   }
}

/*
 *************************************************************************
 *
 * Set each integer value in specified tag array to tag_value where
 * patch level intersects given box array.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::fillTagsFromMappedBoxLevel(
   const int tag_value,
   const tbox::Pointer<hier::PatchLevel> level,
   const int tag_index,
   const hier::Connector& level_to_fill_mapped_box_level,
   const bool interior_only,
   const hier::IntVector& fill_box_growth) const
{
   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(!(level.isNull()));
   TBOX_ASSERT(tag_index == d_tag_indx || tag_index == d_buf_tag_indx);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, *level, fill_box_growth);

   t_fill_tags_from_mapped_box_level->start();

   const hier::OverlapConnectorAlgorithm oca;

   hier::IntVector ratio =
      level_to_fill_mapped_box_level.getHead().getRefinementRatio()
      / level_to_fill_mapped_box_level.getBase().getRefinementRatio();
   TBOX_ASSERT(ratio == level_to_fill_mapped_box_level.getRatio());

   const hier::IntVector gcw =
      hier::IntVector::ceiling(fill_box_growth,
         level_to_fill_mapped_box_level.getRatio());

   for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch> patch = *ip;

      tbox::Pointer<pdat::CellData<int> >
      tag_data = patch->getPatchData(tag_index);

      TBOX_ASSERT(!(tag_data.isNull()));

      const hier::MappedBoxId& mapped_box_id(patch->getMappedBox().getId());
      const hier::BlockId& block_id(mapped_box_id.getBlockId());

      hier::Connector::NeighborSet neighbors;

      oca.extractNeighbors(
         neighbors,
         level_to_fill_mapped_box_level,
         mapped_box_id,
         gcw);

      for (hier::Connector::NeighborSet::const_iterator
           ni = neighbors.begin(); ni != neighbors.end(); ++ni) {
         if (ni->getBlockId() == block_id) {
            hier::Box box = (*ni).getBox();
            box.grow(fill_box_growth);
            box.coarsen(ratio);
            if (interior_only) {
               box = box * tag_data->getBox();
            }
            tag_data->fill(tag_value, box);
         }
      }

   }
   t_fill_tags_from_mapped_box_level->stop();
}

/*
 *************************************************************************
 *
 * Buffer each integer tag with given value on the patch level by the
 * specified buffer size.  Note that the patch data indexed by
 * d_buf_tag_indx is used temporarily to buffer the tag data. The
 * communication of ghost cell (i.e., buffer) information forces all
 * tags on all patch interiors to represent a consistent buffering of
 * the original configuration of tagged cells.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::bufferTagsOnLevel(
   const int tag_value,
   const tbox::Pointer<hier::PatchLevel> level,
   const int buffer_size) const
{
   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(!(level.isNull()));
   TBOX_ASSERT(buffer_size >= 0);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *level);
   /*
    * Start timer for this method.
    */
   t_buffer_tags->start();

   /*
    * Set temporary buffered tags based on buffer width and
    * distance from actual tags.
    */
   const int not_tag = ((tag_value == d_true_tag) ? d_false_tag : d_true_tag);
   for (hier::PatchLevel::Iterator ip1(level); ip1; ip1++) {
      tbox::Pointer<hier::Patch> patch = *ip1;

      tbox::Pointer<pdat::CellData<int> >
      buf_tag_data = patch->getPatchData(d_buf_tag_indx);
      tbox::Pointer<pdat::CellData<int> >
      tag_data = patch->getPatchData(d_tag_indx);

      buf_tag_data->fillAll(not_tag);

      hier::Box interior = patch->getBox();

      for (int bc = buffer_size; bc >= 0; bc--) {

         int fill_val = buffer_size - bc + 1;

         for (pdat::CellIterator ic(interior); ic; ic++) {
            if ((*tag_data)(ic()) == tag_value) {
               hier::Box buf_box(ic() - bc, ic() + bc);
               buf_tag_data->fill(fill_val, buf_box * interior);
            }
         }
      }
   }

   /*
    * Communicate boundary data for buffered tag array so that tags
    * near patch boundaries will become buffered properly.
    */
   const double dummy_time = 0.0;

   t_bdry_fill_tags_comm->start();
   d_bdry_sched_tags[level->getLevelNumber()]->fillData(dummy_time, false);
   t_bdry_fill_tags_comm->stop();

   /*
    * Buffer tags on patch interior according to buffered tag data.
    */
   for (hier::PatchLevel::Iterator ip2(level); ip2; ip2++) {
      tbox::Pointer<hier::Patch> patch = *ip2;

      tbox::Pointer<pdat::CellData<int> > buf_tag_data =
         patch->getPatchData(d_buf_tag_indx);
      tbox::Pointer<pdat::CellData<int> > tag_data =
         patch->getPatchData(d_tag_indx);

      hier::Box buf_tag_box = buf_tag_data->getGhostBox();
      hier::Box tag_box = tag_data->getBox();

      /*
       * Set all interior tags to tag value where buffer tags non-zero.
       */
      for (pdat::CellIterator ic(tag_box); ic; ic++) {
         (*tag_data)(ic()) = ((*buf_tag_data)(ic()) ? tag_value : not_tag);
      }

      /*
       * Set all interior tags in buffers around tags in ghosts.
       */
      for (pdat::CellIterator ic2(buf_tag_box); ic2; ic2++) {
         int tval = (*buf_tag_data)(ic2());
         if (tval > 1) {
            int buf_size = tval - 1;
            hier::Box buf_box(ic2() - buf_size, ic2() + buf_size);
            tag_data->fill(tag_value, buf_box);
         }
      }
   }

   t_buffer_tags->stop();
}

/*
 *************************************************************************
 *
 * Given a patch level, determine an appropriate array of boxes from
 * which a new finer level may be constructed.  That is, find an array
 * of boxes that covers all tags having the specified tag value.  Note
 * that it is assumed that the integer tag arrays have been set
 * properly; i.e., cells have been tagged through error estimation and
 * the tags have been buffered to ensure disturbances remain on fine
 * level until next regrid occurs.  Note that load balancing is
 * performed once an appropriate list of boxes containing the tags is
 * found.  This procedure massages the list of boxes further and then
 * assigns each to a single processor (i.e., the mapping).
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::findRefinementBoxes(
   hier::MappedBoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const int tag_ln,
   const hier::BlockId& block_id) const
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_mb_hierarchy->getFinestLevelNumber()));
   const tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT(!(hierarchy->getPatchLevel(tag_ln).isNull()));
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, new_mapped_box_level, *hierarchy);

   if (s_print_steps == 'y') {
      tbox::plog
      <<
      "MultiblockGriddingAlgorithm::findRefinementBoxes entered with tag_ln = "
      << tag_ln << "\n";
   }

   tbox::SAMRAI_MPI mpi(hierarchy->getMPI());

   const hier::OverlapConnectorAlgorithm oca;
   const hier::MappedBoxLevelConnectorUtils dlbg_edge_utils;

   /*
    * Start timer for this method.
    */
   t_find_refinement->barrierAndStart();

   const hier::MappedBoxLevel& tag_mapped_box_level = *hierarchy->getMappedBoxLevel(tag_ln);

   const int new_ln = tag_ln + 1;


   /*
    * Construct list of boxes covering the true tags on the level.
    * Note that box list will be contained in the bounding box
    * but will not be contained in the list of proper nesting boxes,
    * in general.  So we intersect the box list against the list of
    * nesting boxes.  Note that this may produce boxes which are too
    * small.  Thus, boxes are regrown later.
    */

   hier::IntVector smallest_box_to_refine(d_dim);
   hier::IntVector extend_ghosts(d_dim);
   hier::IntVector largest_patch(d_dim);
   hier::IntVector smallest_patch(d_dim);
   // "true" argument: for_building_finer level = true
   getGriddingParameters(smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts,
      *hierarchy,
      new_ln,
      true);

   const hier::IntVector extend_ghosts_in_tag_space =
      hier::IntVector::ceiling(extend_ghosts,
                               hierarchy->getRatioToCoarserLevel(new_ln));

   tbox::Pointer<hier::PatchLevel> level =
      hierarchy->getPatchLevel(tag_ln);

   const hier::MappedBoxSet& global_mapped_boxes =
      level->getMappedBoxLevel()->getGlobalizedVersion().getGlobalMappedBoxes();

   /*
    * Determine single smallest bounding box for all nesting boxes.
    */
   hier::Box bounding_box(d_dim);
   for (hier::MappedBoxSetSingleBlockIterator mbi(global_mapped_boxes, block_id); mbi.isValid(); ++mbi) {
      bounding_box += mbi->getBox();
   } 

   if (s_print_steps == 'y') {
      tbox::plog
      << "MultiblockGriddingAlgorithm::findRefinementBoxes: clustering\n";
   }

   if (!bounding_box.empty()) {
      t_find_boxes_containing_tags->barrierAndStart();
      hier::IntVector ratio = getRatioToCoarserLevel(new_ln);
      /*
       * Get the max gcw required for cross Connectors from tag to new level.
       * Because new <-> new has a smaller gcw than this ammount, we have
       * enough gcw to compute new <-> new later.
       */
      hier::IntVector gcw =
         hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln + 1);
      d_box_generator->findBoxesContainingTags(
         new_mapped_box_level,
         tag_to_new,
         new_to_tag,
         level, d_tag_indx, d_true_tag, bounding_box,
         smallest_box_to_refine,
         getEfficiencyTolerance(tag_ln),
         getCombineEfficiency(tag_ln),
         gcw,
         block_id);
      t_find_boxes_containing_tags->stop();
   } else {
      /*
       * Initialize empty new_mapped_box_level.
       */
      new_mapped_box_level.initialize(
         level->getMappedBoxLevel()->getRefinementRatio(),
         level->getGridGeometry(),
         level->getMappedBoxLevel()->getMPI());
   }

   if (new_mapped_box_level.getGlobalNumberOfBoxes() > 0) {

      if (s_check_connectors == 'y') {
         /*
          * At this stage, there are no edges to periodic images yet, so
          * don't check for them.
          */
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag,
               false,
               true,
               true) == 0);
         TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new,
               false,
               true,
               true) == 0);
      }

      t_box_massage->start();

      {
         if (s_print_steps == 'y') {
            tbox::plog
            <<
            "MultiblockGriddingAlgorithm::findRefinementBoxes: enforcing overflow nesting\n";
         }

         t_limit_overflow->start();
         /*
          * Do not allow the new mapped_box_level to overflow the tag mapped_box_level.
          * If we want to allow the overflow, we have to add the
          * overflow ammount to gcw(tag->new).  Such additions
          * may make the ABR algorithm slower, because more
          * non-contributing processors would have to be included
          * in the contributing group (unless ABR keep track of
          * the non-contributing processors and don't seek tag
          * histogram from them).
          */
         const hier::IntVector allowed_overflow(d_dim, 0);
         hier::MappedBoxLevel nested_mapped_box_level(d_dim);
         hier::Connector unnested_to_nested;
         makeOverflowNestingMap(new_mapped_box_level,
            new_to_tag,
            nested_mapped_box_level,
            unnested_to_nested,
            block_id,
            allowed_overflow);
         t_use_overflow_map->barrierAndStart();
         t_modify_connector->start();
         const hier::MappingConnectorAlgorithm mca;
         mca.modify(tag_to_new,
            new_to_tag,
            unnested_to_nested,
            &new_mapped_box_level);
         t_modify_connector->stop();
         t_use_overflow_map->stop();
         t_limit_overflow->stop();

         if (s_check_overflow_nesting == 'y') {
            bool locally_nested = false;
            bool nested = dlbg_edge_utils.baseNestsInHead(
                  &locally_nested,
                  new_mapped_box_level,
                  tag_mapped_box_level,
                  allowed_overflow,
                  hier::IntVector::getZero(d_dim),
                  hier::IntVector::getZero(d_dim),
                  &hierarchy->getDomainSearchTree(block_id));
            if (!nested) {
               TBOX_ERROR(
                  "Error found."
                  << "Failed overflow nesting: new mapped_box_level does not nest in tagged mapped_box_level.\n"
                  << "allowed_overflow = " << allowed_overflow << std::endl
                  << "Local nestedness = " << locally_nested << std::endl
                  << "tag_mapped_box_level:\n" << tag_mapped_box_level.format("", 2)
                  << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
                  << "tag_to_new:\n" << tag_to_new.format("", 2)
                  << "new_to_tag:\n" << new_to_tag.format("", 2));
            }
         }
      }

      if (d_enforce_proper_nesting) {

         if (s_print_steps == 'y') {
            tbox::plog
            <<
            "MultiblockGriddingAlgorithm::findRefinementBoxes: enforcing proper nesting\n";
         }
         t_enforce_nesting->barrierAndStart();

         hier::Connector unnested_to_nested;
         /*
          * Copy the input new_mapped_box_level so we can later connect
          * new_mapped_box_level to its properly nested version.
          */
         hier::MappedBoxLevel nested_mapped_box_level(d_dim);
         makeProperNestingMap(new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            new_ln,
            block_id,
            nested_mapped_box_level,
            unnested_to_nested);
         t_use_nesting_map->start();
         t_modify_connector->start();
         const hier::MappingConnectorAlgorithm mca;
         mca.modify(tag_to_new,
            new_to_tag,
            unnested_to_nested,
            &new_mapped_box_level);
         t_modify_connector->stop();
         t_use_nesting_map->stop();

         t_enforce_nesting->stop();

         if (tag_ln == d_base_ln && s_check_proper_nesting == 'y') {
            /*
             * Tag level will be regridded when we exit the current
             * recursion if tag_ln is not d_base_ln, so do not check
             * proper nesting in that case.
             *
             * Check that the new mapped_box_level nest in the tag
             * level (tag_ln).
             */
            hier::IntVector required_nesting(d_dim);
            if (tag_ln > 0) {
               required_nesting =
                  hier::IntVector(d_dim, getProperNestingBuffer(tag_ln));
               required_nesting *= getRatioToCoarserLevel(new_ln);
            } else {
               required_nesting =
                  hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_dim);
            }
            bool locally_nests = false;
            const bool new_nests_in_tag =
               dlbg_edge_utils.baseNestsInHead(
                  &locally_nests,
                  new_mapped_box_level,
                  tag_mapped_box_level,
                  required_nesting,
                  hier::IntVector::getZero(d_dim),
                  hier::IntVector::getZero(d_dim),
                  &hierarchy->getPeriodicDomainSearchTree(block_id));
            if (!new_nests_in_tag) {
               tbox::perr
               << "MultiblockGriddingAlgorithm: new mapped_box_level\n"
               << "at ln=" << new_ln
               << " does not properly nest in\n"
               << "tag level at tag_ln=" << tag_ln
               << " by the required nesting buffer of "
               << required_nesting
               << ".\nLocal nestingness: " << locally_nests
               << ".\nWriting MappedBoxLevels out to log file."
               << std::endl;
               tbox::plog << "Proper nesting violation with new_mapped_box_level of\n" << new_mapped_box_level.format("N->", 2);
               tbox::plog << "Proper nesting violation with tag mapped_box_level of\n" << tag_mapped_box_level.format("T->", 2);
               hier::MappedBoxLevel external(d_dim);
               hier::Connector tmp_new_to_tag(
                  new_mapped_box_level,
                  tag_mapped_box_level,
                  required_nesting);
               oca.findOverlaps(tmp_new_to_tag);
               tbox::plog << "tmp_new_to_tag:\n" << tmp_new_to_tag.format("NT->", 3);
               hier::Connector new_to_external;
               dlbg_edge_utils.computeExternalParts(
                  external,
                  new_to_external,
                  tmp_new_to_tag,
                  -required_nesting,
                  hierarchy->getDomainSearchTree(block_id));
               tbox::plog << "External parts:\n" << new_to_external.format("NE->", 3);
               TBOX_ERROR(
                  "Internal library error: Failed to produce proper nesting.");
            }
         }
      }

      if (d_extend_to_domain_boundary) {

         if (s_print_steps == 'y') {
            tbox::plog
            <<
            "MultiblockGriddingAlgorithm::findRefinementBoxes: extending nodes\n";
         }

         t_extend_to_domain_boundary->barrierAndStart();
         extendMappedBoxesToDomainBoundary(
            new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            hierarchy->getConnector(tag_ln, tag_ln),
            level->getPhysicalDomain(block_id),
            extend_ghosts_in_tag_space,
            *hierarchy,
            tag_ln);
         t_extend_to_domain_boundary->stop();
      }

      bool allow_patches_smaller_than_minimum_size_to_prevent_overlaps
         = hierarchy->allowPatchesSmallerThanMinimumSize();

      // BTNG: these if-else blocks can be significantly simplified by factoring.
      if (!allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
         if (s_print_steps == 'y') {
            tbox::plog
            <<
            "MultiblockGriddingAlgorithm::findRefinementBoxes: growing boxes\n";
         }
         t_extend_within_domain->start();
         hier::MappedBoxLevel grown_mapped_box_level(d_dim);
         hier::Connector new_to_grown;
         makeBoxGrowingMap(
            new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            tag_ln,
            block_id,
            grown_mapped_box_level,
            new_to_grown,
            smallest_box_to_refine);
         t_modify_connector->start();
         const hier::MappingConnectorAlgorithm mca;
         mca.modify(tag_to_new,
            new_to_tag,
            new_to_grown,
            &new_mapped_box_level);
         t_modify_connector->stop();
         t_extend_within_domain->stop();
      } else {
         const hier::IntVector& periodic_dirs(
            hierarchy->getGridGeometry()->getPeriodicShift(
               hier::IntVector::getOne(d_dim)));

         bool need_to_grow = false;
         hier::IntVector min_size(hier::IntVector::getOne(d_dim));
         for (int i = 0; i < d_dim.getValue(); i++) {
            if (periodic_dirs(i)) {
               need_to_grow = true;
               min_size(i) = smallest_box_to_refine(i);
            }
         }

         if (need_to_grow) {
            t_extend_within_domain->start();
            hier::MappedBoxLevel grown_mapped_box_level(d_dim);
            hier::Connector new_to_grown;
            makeBoxGrowingMap(
               new_mapped_box_level,
               tag_to_new,
               new_to_tag,
               tag_ln,
               block_id,
               grown_mapped_box_level,
               new_to_grown,
               min_size);
            t_modify_connector->start();
            const hier::MappingConnectorAlgorithm mca;
            mca.modify(tag_to_new,
               new_to_tag,
               new_to_grown,
               &new_mapped_box_level);
            t_modify_connector->stop();
            t_extend_within_domain->stop();
         }
      }

      t_box_massage->stop();
   }

   mpi.Barrier();

   t_find_refinement->stop();
}

/*
 *******************************************************************
 * Load balance a level and refine it to the new index space.
 *******************************************************************
 */

void MultiblockGriddingAlgorithm::loadBalanceAndRefineBoxes(
   hier::MappedBoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const int tag_ln,
   const tbox::RankGroup& rank_group,
   const hier::BlockId& block_id) const
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_mb_hierarchy->getFinestLevelNumber()));
   const tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);
   TBOX_ASSERT(!(hierarchy.isNull()));
   TBOX_ASSERT(!(hierarchy->getPatchLevel(tag_ln).isNull()));
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, new_mapped_box_level, *hierarchy);

   const hier::OverlapConnectorAlgorithm oca;

   tbox::SAMRAI_MPI mpi(hierarchy->getMPI());

   const int new_ln = tag_ln + 1;

   /*
    * Start timer for this method.
    */
   t_find_refinement->barrierAndStart();

   const hier::MappedBoxLevel& tag_mapped_box_level =
      *hierarchy->getMappedBoxLevel(tag_ln);

   hier::IntVector ratio(getRatioToCoarserLevel(new_ln));

   if (new_mapped_box_level.getGlobalNumberOfBoxes() > 0) {

      if (d_load_balance) {
         if (s_print_steps == 'y') {
            tbox::plog
            <<
            "MultiblockGriddingAlgorithm::findRefinementBoxes: load balancing\n";
         }

         hier::IntVector smallest_box_to_refine(d_dim);
         hier::IntVector extend_ghosts(d_dim);
         hier::IntVector largest_patch(d_dim);
         hier::IntVector smallest_patch(d_dim);
         // "true" argument: for_building_finer level = true
         getGriddingParameters(smallest_patch,
            smallest_box_to_refine,
            largest_patch,
            extend_ghosts,
            *hierarchy,
            new_ln,
            true);

         const hier::IntVector smallest_patch_in_tag_space =
            hier::IntVector::ceiling(smallest_patch,
                                     hierarchy->getRatioToCoarserLevel(new_ln));
         const hier::IntVector largest_patch_in_tag_space =
            largest_patch/hierarchy->getRatioToCoarserLevel(new_ln);

         const hier::IntVector extend_ghosts_in_tag_space =
            hier::IntVector::ceiling(extend_ghosts,
                                     hierarchy->getRatioToCoarserLevel(new_ln));

         t_load_balance->barrierAndStart();
         t_load_balance_setup->start();

         const hier::IntVector& patch_cut_factor =
            hier::IntVector::getOne(d_dim);

         t_load_balance_setup->stop();

         d_load_balancer->loadBalanceMappedBoxLevel(
            new_mapped_box_level,
            new_to_tag,
            tag_to_new,
            hierarchy,
            new_ln,
            new_to_tag, // FIXME: try using finer as the attractor.
            tag_to_new, // FIXME: try using finer as the attractor.
            smallest_patch_in_tag_space,
            largest_patch_in_tag_space,
            d_singleblock_domain_mapped_box_level[block_id.getBlockValue()],
            extend_ghosts_in_tag_space,
            patch_cut_factor,
            rank_group);

         t_load_balance->stop();

         if (s_check_connectors == 'y') {
            tbox::plog << "MultiblockGriddingAlgorithm checking new-tag"
                       << std::endl;
            int errs = 0;
            if (oca.checkOverlapCorrectness(new_to_tag, false, true, true)) {
               ++errs;
               tbox::perr << "Error found in new_to_tag!\n";
            }
            if (oca.checkOverlapCorrectness(tag_to_new, false, true, true)) {
               ++errs;
               tbox::perr << "Error found in tag_to_new!\n";
            }
            if (new_to_tag.checkTransposeCorrectness(tag_to_new)) {
               ++errs;
               tbox::perr << "Error found in new-tag transpose!\n";
            }
            if (errs != 0) {
               TBOX_ERROR(
                  "Errors found after using load balance map."
                  << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
                  << "tag_mapped_box_level:\n" << tag_mapped_box_level.format("", 2)
                  << "new_to_tag:\n" << new_to_tag.format("", 2)
                  << "tag_to_new:\n" << tag_to_new.format("", 2));
            }
         }
      }

      if (d_sequentialize_patch_indices) {
         if (s_check_connectors == 'y') {
            tbox::plog << "MultiblockGriddingAlgorithm begin sorting nodes."
                       << std::endl;
         }
         /*
          * Sequentializing the indices is not required.
          * But while transitioning between the old global-boxes
          * approach and the new DLBG approaches, this makes the
          * indices the same as the patch indices, which
          * helps in code development and debugging.
          */
         sortNodes(new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            false,
            true);
         if (s_check_connectors == 'y') {
            tbox::plog << "MultiblockGriddingAlgorithm end sorting nodes."
                       << std::endl;
         }
      }

      /*
       * Add periodic image MappedBoxes to new_mapped_box_level and add edges
       * incident on those nodes.
       */
      if (s_check_connectors == 'y') {
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag) == 0);
         TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new) == 0);
      }

      /*
       * We have been working with new_mapped_box_level in the tag_mapped_box_level's index space.
       * Now, refine it so we can build the new level.
       */
      refineNewMappedBoxLevel(new_mapped_box_level,
         tag_to_new,
         new_to_tag,
         ratio);

      if (s_check_connectors == 'y') {
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag) == 0);
         TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new) == 0);
      }

   } else {

      /*
       * On return, new_mapped_box_level should be initialized if we generated boxes,
       * deallocated if we didnt.
       */
      new_mapped_box_level.clear();

   }

   mpi.Barrier();

   t_find_refinement->stop();
}


/*
 ***********************************************************************
 * Set RankGroups, allocating ranks to blocks based on a fraction of work
 ***********************************************************************
 */
void MultiblockGriddingAlgorithm::setRankGroups(
   tbox::Array<tbox::RankGroup>& rank_groups,
   const tbox::Array<int>& work_on_block,
   const tbox::SAMRAI_MPI& mpi) const
{
   const unsigned int nblocks = rank_groups.size();
   const int num_ranks = mpi.getSize();

   int active_blocks = 0;
   int total_work = 0;
   unsigned int max_active_block = tbox::MathUtilities<int>::getMax();
   for (unsigned int nb = 0; nb < nblocks; nb++) {
      if (work_on_block[nb]) {
         active_blocks++;
         total_work += work_on_block[nb];
         max_active_block = nb;
      }
   }

   if (active_blocks > 1 && num_ranks > active_blocks) {

      int min_rank = 0;
      int max_rank;
      for (unsigned int nb = 0; nb < nblocks; nb++) {

         if (work_on_block[nb]) {
            if (nb != max_active_block) {
               double work_fraction = (double)work_on_block[nb]/
                                      (double)total_work;

               double est_ranks = (work_fraction * num_ranks);
               int procs_for_block = (int)(est_ranks + 0.5) ;
               if (procs_for_block == 0) {
                  procs_for_block = 1;
               }

               max_rank = min_rank + procs_for_block - 1;
               if (max_rank >= num_ranks) {
                  max_rank = num_ranks - 1;
               }
            } else {
               max_rank = num_ranks - 1;
            }

            if (min_rank > max_rank) {
               min_rank = max_rank;
            }

            rank_groups[nb].setMinMax(min_rank, max_rank);

            min_rank = max_rank + 1;
         }
      }
   }
}


/*
 ***********************************************************************
 ***********************************************************************
 */
void MultiblockGriddingAlgorithm::sortNodes(
   hier::MappedBoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   bool sort_by_corners,
   bool sequentialize_global_indices) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, new_mapped_box_level);

   t_sort_nodes->start();

   const hier::OverlapConnectorAlgorithm oca;
   const hier::MappingConnectorAlgorithm mca;

   hier::Connector sorting_map;
   hier::MappedBoxLevel seq_mapped_box_level(d_dim);
   hier::MappedBoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      seq_mapped_box_level,
      sorting_map,
      new_mapped_box_level,
      sort_by_corners /* sort nodes by corners */,
      sequentialize_global_indices /* do not sequentialize indices globally */);

   if (0) {
      // Check sorting_map before using it.
      int errs = 0;
      if (mca.findMappingErrors(sorting_map, 'n')) {
         ++errs;
         tbox::perr << "sorting_map is not a valid mapping!\n";
      }
      if (oca.checkOverlapCorrectness(sorting_map, false, true, true)) {
         ++errs;
         tbox::perr << "Error found in sorting_map!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "unsequentialized mapped_box_level:\n" << new_mapped_box_level.format("", 2)
            << "sequentialized mapped_box_level:\n" << seq_mapped_box_level.format("", 2)
            << "sorting_map:\n" << sorting_map.format("", 2));
      }
   }
   t_modify_connector->start();
   mca.modify(tag_to_new,
      new_to_tag,
      sorting_map,
      &new_mapped_box_level);
   t_modify_connector->stop();
   if (0) {
      // Check result of mapping.
      int errs = 0;
      if (oca.checkOverlapCorrectness(tag_to_new, true, false)) {
         ++errs;
         tbox::perr << "Error found in tag_to_new!\n";
      }
      if (oca.checkOverlapCorrectness(new_to_tag, true, false)) {
         ++errs;
         tbox::perr << "Error found in new_to_tag!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors found after sequentializing nodes."
            << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
            << "tag mapped_box_level:\n" << tag_to_new.getBase().format("", 2)
            << "tag_to_new:\n" << tag_to_new.format("", 2)
            << "new_to_tag:\n" << new_to_tag.format("", 2));
      }
   }

#if 0
   TBOX_ASSERT(new_to_tag.checkOverlapCorrectness(false, true) == 0);
   TBOX_ASSERT(tag_to_new.checkOverlapCorrectness(false, true) == 0);
#endif

   t_sort_nodes->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MultiblockGriddingAlgorithm::sortNodes(
   hier::MappedBoxLevel& new_mapped_box_level,
   bool sequentialize_global_indices) const
{
   (void)sequentialize_global_indices;

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, new_mapped_box_level);
   t_sort_nodes->start();

   const hier::MappedBoxSet& global_mapped_boxes =
      new_mapped_box_level.getGlobalizedVersion().getGlobalMappedBoxes();

   const hier::LocalId& initial_sequential_index = hier::LocalId::getZero();
   hier::LocalId last_index = initial_sequential_index - 1;

   const int rank = new_mapped_box_level.getRank();

   hier::MappedBoxSet new_mapped_boxes;
   for (hier::MappedBoxSet::const_iterator bi =
        global_mapped_boxes.begin(); bi != global_mapped_boxes.end(); ++bi)
   {
      const hier::MappedBox& cur_mapped_box = *bi;
      if (cur_mapped_box.getOwnerRank() == rank) {
         const hier::MappedBox new_mapped_box(cur_mapped_box.getBox(),
                                        ++last_index,
                                        cur_mapped_box.getOwnerRank(),
                                        cur_mapped_box.getBlockId(),
                                        cur_mapped_box.getPeriodicId());
         new_mapped_boxes.insert(new_mapped_boxes.end(), new_mapped_box);
      } else {
         ++last_index;
      }
   }

   new_mapped_box_level.swapInitialize(
      new_mapped_boxes,
      new_mapped_box_level.getRefinementRatio(),
      new_mapped_box_level.getGridGeometry(),
      new_mapped_box_level.getMPI(),
      hier::MappedBoxLevel::DISTRIBUTED);

   t_sort_nodes->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MultiblockGriddingAlgorithm::refineNewMappedBoxLevel(
   hier::MappedBoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const hier::IntVector& ratio) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, new_mapped_box_level, ratio);

   hier::MappedBoxSet refined_nodes;
   new_mapped_box_level.getMappedBoxes().refine(refined_nodes, ratio);
   new_mapped_box_level.swapInitialize(
      refined_nodes,
      ratio * new_mapped_box_level.getRefinementRatio(),
      new_mapped_box_level.getGridGeometry(),
      new_mapped_box_level.getMPI());

   new_to_tag.initialize(
      new_mapped_box_level,
      new_to_tag.getHead(),
      ratio * new_to_tag.getConnectorWidth(),
      new_to_tag.getNeighborhoodSets(),
      hier::MappedBoxLevel::DISTRIBUTED);

   hier::NeighborhoodSet tag_eto_new;
   tag_to_new.swapInitialize(
      tag_to_new.getBase(),
      new_mapped_box_level,
      tag_to_new.getConnectorWidth(),
      tag_eto_new,
      hier::MappedBoxLevel::DISTRIBUTED);
   for (hier::NeighborhoodSet::const_iterator ei = tag_eto_new.begin();
        ei != tag_eto_new.end(); ++ei) {
      ei->second.refine(tag_eto_new[ei->first], ratio);
   }
   tag_to_new.swapInitialize(
      tag_to_new.getBase(),
      new_mapped_box_level,
      tag_to_new.getConnectorWidth(),
      tag_eto_new,
      hier::MappedBoxLevel::DISTRIBUTED);
#if 0
   const hier::OverlapConnectorAlgorithm oca;
   TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new) == 0);
   TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag) == 0);
#endif
}

/*
 *************************************************************************
 *
 * Extend nodes to domain boundary if they are too close.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::extendMappedBoxesToDomainBoundary(
   hier::MappedBoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const hier::Connector& tag_to_tag,
   const hier::BoxList& physical_domain_list,
   const hier::IntVector& extend_ghosts,
   const hier::PatchHierarchy& hierarchy,
   const int new_ln) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, new_mapped_box_level);

   (void)hierarchy;
   (void)new_ln;

   tbox::SAMRAI_MPI mpi(new_mapped_box_level.getMPI());

   /*
    * Extend boxes to domain boundary if they are too close to one.  I
    * think there is no need to modify connectivities when a box is
    * extended to the domain boundary.  There is no increased overlap
    * to a finer level, because the finer level already nests in the
    * nodes being extended.  There is increased overlap with the
    * coarser level, but it should not result in additional edges,
    * because there would not be enough room for an unseen coarse MappedBox
    * to live in the small gap across which the MappedBox is being extended.
    */
   const hier::MappedBoxSet& before_nodes =
      new_mapped_box_level.getMappedBoxes();

   hier::MappedBoxLevel after_mapped_box_level(d_dim);
   after_mapped_box_level.initialize(
      new_mapped_box_level.getRefinementRatio(),
      new_mapped_box_level.getGridGeometry(),
      new_mapped_box_level.getMPI());
   hier::NeighborhoodSet before_eto_after, after_eto_before;
   int n_local_nodes_extended = 0;
   for (hier::MappedBoxSet::const_iterator
        nn = before_nodes.begin(); nn != before_nodes.end(); ++nn) {
      const hier::MappedBox& before_mapped_box = *nn;
      hier::MappedBox after_mapped_box = before_mapped_box;
      bool extended =
         hier::BoxUtilities::extendBoxToDomainBoundary(
            after_mapped_box.getBox(),
            physical_domain_list,
            extend_ghosts);
      if (extended) {
         ++n_local_nodes_extended;
      }
      after_mapped_box_level.addMappedBox(after_mapped_box);
      before_eto_after[
         before_mapped_box.getId()].insert(after_mapped_box);
      after_eto_before[
         after_mapped_box.getId()].insert(before_mapped_box);
   }
   hier::Connector before_to_after, after_to_before;
   before_to_after.initialize(
      new_mapped_box_level,
      after_mapped_box_level,
      hier::IntVector::getZero(d_dim),
      before_eto_after);
   before_to_after.setConnectorType(hier::Connector::BASE_GENERATED);
   after_to_before.initialize(
      after_mapped_box_level,
      new_mapped_box_level,
      hier::IntVector::getZero(d_dim),
      after_eto_before);
   after_to_before.setConnectorType(hier::Connector::BASE_GENERATED);

   int n_global_nodes_extended;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(
         &n_local_nodes_extended,
         &n_global_nodes_extended,
         1,
         MPI_INT,
         MPI_SUM);
   } else {
      n_global_nodes_extended = n_local_nodes_extended;
   }
   if (n_global_nodes_extended > 0) {
      /*
       * If there are changed nodes, modify the tag<->new Connectors.
       * Because the modified nodes were extended and the new nodes
       * may not sufficiently nest in the tag mapped_box_level, tag<->new may be
       * missing some edges.  Fix this using a peer bridge.
       */
      const hier::MappingConnectorAlgorithm mca;
      mca.modify(tag_to_new,
         new_to_tag,
         before_to_after,
         &new_mapped_box_level);
      hier::Connector tmp_new_to_tag = new_to_tag;
      hier::Connector tmp_tag_to_new = tag_to_new;
      const hier::OverlapConnectorAlgorithm oca;
      oca.bridge(new_to_tag,
         tag_to_new,
         tmp_new_to_tag,
         tag_to_tag,
         tag_to_tag,
         tmp_tag_to_new);
      oca.shrinkConnectorWidth(new_to_tag, tmp_new_to_tag.getConnectorWidth());
      oca.shrinkConnectorWidth(tag_to_new, tmp_tag_to_new.getConnectorWidth());
#if 0
      TBOX_WARNING("Performing extensive error checking due to using new code!");
      TBOX_ASSERT(new_to_tag.checkOverlapCorrectness() == 0);
      TBOX_ASSERT(tag_to_new.checkOverlapCorrectness() == 0);
#endif
   }
}

/*
 *************************************************************************
 * Make a map that can be used to enforce overflow nesting.
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::makeOverflowNestingMap(
   const hier::MappedBoxLevel& unnested_mapped_box_level,
   const hier::Connector& unnested_to_nominal,
   hier::MappedBoxLevel& nested_mapped_box_level,
   hier::Connector& unnested_to_nested,
   const hier::BlockId& block_id,
   const hier::IntVector& allowed_overflow) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim,
      unnested_mapped_box_level,
      nested_mapped_box_level,
      allowed_overflow);

   (void)unnested_mapped_box_level;

   tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);

   t_make_overflow_map->start();

   hier::MappedBoxLevel violator_mapped_box_level(d_dim);
   hier::Connector unnested_to_violator;
   t_make_overflow_map_compute->start();
   hier::MappedBoxLevelConnectorUtils edge_utils;
   t_compute_external_parts->start();
   edge_utils.computeExternalParts(
      violator_mapped_box_level,
      unnested_to_violator,
      unnested_to_nominal,
      allowed_overflow,
      hierarchy->getDomainSearchTree(block_id));
   t_compute_external_parts->stop();
   t_make_overflow_map_compute->stop();

   TBOX_ASSERT(unnested_to_violator.isLocal());

   t_make_overflow_map_convert->start();
   edge_utils.makeRemainderMap(
      nested_mapped_box_level,
      unnested_to_nested,
      unnested_to_violator);
   t_make_overflow_map_convert->stop();
   t_make_overflow_map->stop();
}

/*
 *************************************************************************
 * Make a map that can be used to enforce proper nesting.
 * @param unnested_ln Level number for refinenement ratio of unnested_mapped_box_level.
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::makeProperNestingMap(
   const hier::MappedBoxLevel& unnested_mapped_box_level,
   const hier::Connector& hierarchy_to_unnested,
   const hier::Connector& unnested_to_hierarchy,
   const int unnested_ln,
   const hier::BlockId& block_id,
   hier::MappedBoxLevel& nested_mapped_box_level,
   hier::Connector& unnested_to_nested) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim,
      unnested_mapped_box_level,
      nested_mapped_box_level,
      nested_mapped_box_level);

   t_make_nesting_map->start();

   hier::Connector unnested_to_violator;
   hier::MappedBoxLevel violator(d_dim);
   t_make_nesting_map_compute->start();
   computeNestingViolator(unnested_mapped_box_level,
      violator,
      unnested_to_violator,
      unnested_to_hierarchy,
      hierarchy_to_unnested,
      unnested_ln,
      block_id);
   t_make_nesting_map_compute->stop();


   /*
    * unnested_to_violator is the Connector from the nodes
    * that violate nesting to their violating parts.
    * Convert it to the mapping from unnested to nested.
    */
   const hier::MappedBoxLevelConnectorUtils edge_utils;
   t_make_nesting_map_convert->start();
   edge_utils.makeRemainderMap(
      nested_mapped_box_level,
      unnested_to_nested,
      unnested_to_violator);
   t_make_nesting_map_convert->stop();

   t_make_nesting_map->stop();
}

/*
 *************************************************************************
 * Make a map from a mapped_box_level to parts of that mapped_box_level
 * that violate proper nesting.
 *
 * The violating mapped_boxes are found by comparing candidate
 * mapped_boxes to d_to_nesting_complement's head MappedBoxLevel.  Boxes
 * inside the nesting complement violate nesting.
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::computeNestingViolator(
   const hier::MappedBoxLevel& candidate,
   hier::MappedBoxLevel& violator,
   hier::Connector& candidate_to_violator,
   const hier::Connector& candidate_to_hierarchy,
   const hier::Connector& hierarchy_to_candidate,
   const int candidate_ln,
   const hier::BlockId& block_id) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, candidate, violator);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   // Check requirements on arguments.
   TBOX_ASSERT(candidate_to_hierarchy.getRatio() ==
      hier::IntVector::getOne(d_dim));
   TBOX_ASSERT(hierarchy_to_candidate.getRatio() ==
      hier::IntVector::getOne(d_dim));

   t_compute_nesting_violator->start();

   const int block_number = block_id.getBlockValue();
   const tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);

   const hier::MappedBoxLevelConnectorUtils edge_utils;
   const hier::OverlapConnectorAlgorithm oca;

   const hier::MappedBoxSet& candidate_mapped_boxes = candidate.getMappedBoxes();

   /*
    * Bridge candidate to d_nesting_complement.
    */
   hier::Connector candidate_to_complement, complement_to_candidate;
   oca.bridge(candidate_to_complement,
      complement_to_candidate,
      candidate_to_hierarchy,
      d_to_nesting_complement[block_number][candidate_ln - 1],
      d_from_nesting_complement[block_number][candidate_ln - 1],
      hierarchy_to_candidate);

   edge_utils.computeInternalParts(
      violator,
      candidate_to_violator,
      candidate_to_complement,
      zero_vector,
      hierarchy->getDomainSearchTree(block_id));
   hier::NeighborhoodSet candidate_eto_violator;
   candidate_to_violator.swapInitialize(
      candidate_to_violator.getBase(),
      candidate_to_violator.getHead(),
      candidate_to_violator.getConnectorWidth(),
      candidate_eto_violator);
   /*
    * Above step ignored the domain complement components of nesting
    * definition (by necessity).  Where the candidate overlaps
    * with the domain complement, it violates nesting.
    */
   tbox::Pointer<hier::MappedBoxTree> refined_domain_complement_tree =
      d_domain_complement_tree[block_number].createRefinedTree(candidate.getRefinementRatio());
   for (hier::MappedBoxSet::const_iterator ni = candidate_mapped_boxes.begin();
        ni != candidate_mapped_boxes.end(); ++ni) {
      const hier::MappedBox& cmb = *ni;
      if (cmb.getBlockId() != block_id) {
         continue;
      }
      hier::BoxList addl_violators(cmb.getBox());
      addl_violators.intersectBoxes(*refined_domain_complement_tree);
      if (!addl_violators.isEmpty()) {
         if (candidate_eto_violator.find(cmb.getId()) !=
             candidate_eto_violator.end()) {
            /*
             * Remove parts that we already know are non-nesting.
             * Leftovers are non-nesting parts not found using
             * candidate_to_complement.
             */
            hier::Connector::NeighborSet& current_violators =
               candidate_eto_violator[cmb.getId()];
            for (hier::Connector::NeighborSet::const_iterator na =
                    current_violators.begin();
                 na != current_violators.end() && !addl_violators.isEmpty();
                 ++na) {
               addl_violators.removeIntersections(na->getBox());
            }
            if (!addl_violators.isEmpty()) {
               for (hier::BoxList::Iterator bi(addl_violators); bi; bi++) {
                  hier::MappedBoxSet::iterator new_violator = violator.addBox(
                        *bi, hier::BlockId(block_number));
                  current_violators.insert(*new_violator);
               }
            }
         }
      }
   }

   candidate_to_violator.swapInitialize(
      candidate_to_violator.getBase(),
      candidate_to_violator.getHead(),
      candidate_to_violator.getConnectorWidth(),
      candidate_eto_violator);

   t_compute_nesting_violator->stop();
}

/*
 *************************************************************************
 * Precompute data used to define proper nesting.  Data is associated
 * with level number ln, to be used for constructing level number ln+1.
 *
 * If ln > d_base_ln, assume data at ln-1 is already set.
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::computeNestingData(
   const int ln,
   const hier::BlockId& block_id)
{
   TBOX_ASSERT(d_base_ln >= 0 && ln >= d_base_ln);

   const int block_number = block_id.getBlockValue();
   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   const tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);

   const hier::MappedBoxTree& domain_search_tree = hierarchy->getDomainSearchTree(block_id);

   const hier::MappedBoxLevelConnectorUtils edge_utils;
   const hier::OverlapConnectorAlgorithm oca;

   hier::IntVector tmp_gcw(d_dim);

   if (ln == d_base_ln) {
      /*
       * At the base level, nesting domain is domain of d_base_ln,
       * shrunken by d_proper_nesting_buffer[d_base_ln].
       */
      hier::MappedBoxLevel& proper_nesting_complement =
         d_proper_nesting_complement[block_number][ln];
      const hier::Connector& peer_connector =
         hierarchy->getConnector(ln, ln);

      edge_utils.computeExternalParts(
         proper_nesting_complement,
         d_to_nesting_complement[block_number][ln],
         peer_connector,
         hier::IntVector(d_dim, -getProperNestingBuffer(ln)),
         hierarchy->getDomainSearchTree(block_id));

      d_from_nesting_complement[block_number][ln].initializeToLocalTranspose(
         d_to_nesting_complement[block_number][ln]);

      /*
       * computeExternalParts neglected parts of the proper nesting
       * complement that are outside the domain (by necessity).  We
       * must also build a search tree of the domain's complement which
       * will be used to detect non-nesting cells lying outside the
       * domain.
       */
      hier::Box universe(d_singleblock_domain_mapped_box_level[block_id.getBlockValue()].getGlobalBoundingBox(block_id.getBlockValue()));
      universe.grow(hier::IntVector::getOne(d_dim));
      hier::BoxList domain_complement_list(universe);
      domain_complement_list.removeIntersections(domain_search_tree);
      std::vector<hier::MappedBox> domain_complement_vector;
      hier::MappedBoxContainerUtils::convertBoxListToMappedBoxVector(
         domain_complement_list,
         domain_complement_vector,
         block_id);
      d_domain_complement_tree[block_number].generateTree(domain_complement_vector);

   } else {

      TBOX_ASSERT(d_to_nesting_complement[block_number][ln - 1].isInitialized());

      /*
       * How to build d_proper_nesting_complement[ln] and connect it to level ln:
       *
       *                  (new
       *               Connectors)
       |
       *            Mapped  |   Proper
       *              box   |   nesting
       *            levels  | complements
       |    |     |
       *               v    v     v
       *              ln <-----> ln
       *               ^        ^
       |       /
       |      /
       |     /
       * (existing     |    / <--(temporary
       * Connectors)--> |   /     Connectors)
       |  /
       | /
       |/
       *               v
       *            ln-1 <-----> ln-1
       *                    ^
       |
       *                (existing
       *               Connectors)
       *
       * We have existing Connectors between ln and ln-1 and also between
       * level ln-1 and the nesting complement at ln-1.
       *
       * 1. Build the complement at ln from the complement at ln-1.
       *
       * 2. Build the temporary Connector from level ln-1 to
       *   complements at ln by using the fact that the complement at ln
       *   is similar to the one at ln-1.
       *
       * 3. Bridge for the new Connector from level ln to the
       *   complement at ln, using the temporary Connector.
       */

      /*
       * 1. Build d_proper_nesting_complement[ln] from d_proper_nesting_complement[ln-1].
       */
      d_proper_nesting_complement[block_number][ln].initialize(
         hierarchy->getMappedBoxLevel(ln)->getRefinementRatio(),
         hierarchy->getMappedBoxLevel(ln)->getGridGeometry(),
         d_to_nesting_complement[block_number][ln - 1].getMPI());
      const hier::MappedBoxSet& lnm1_complement_mapped_boxes =
         d_proper_nesting_complement[block_number][ln - 1].getMappedBoxes();
      for (hier::MappedBoxSet::const_iterator ni =
              lnm1_complement_mapped_boxes.begin();
           ni != lnm1_complement_mapped_boxes.end(); ++ni) {
         hier::MappedBox tmp_mapped_box = *ni;
         TBOX_ASSERT(!tmp_mapped_box.isPeriodicImage());
         tmp_mapped_box.getBox().refine(getRatioToCoarserLevel(ln));
         tmp_mapped_box.getBox().grow(
            hier::IntVector(d_dim, getProperNestingBuffer(ln)));
         d_proper_nesting_complement[block_number][ln].addMappedBox(
            tmp_mapped_box);
      }

      /*
       * 2. Temporarily connect level ln-1 and d_proper_nesting_complement[ln].
       */
      tmp_gcw = d_to_nesting_complement[block_number][ln - 1].getConnectorWidth();
      tmp_gcw -= hier::IntVector::ceiling(
            hier::IntVector(d_dim, getProperNestingBuffer(ln - 1)),
            getRatioToCoarserLevel(ln));
      hier::NeighborhoodSet lnm1_eto_ln_complement;
      const hier::NeighborhoodSet& lnm1_eto_lnm1_complement =
         d_to_nesting_complement[block_number][ln - 1].getNeighborhoodSets();
      for (hier::NeighborhoodSet::const_iterator ei = lnm1_eto_lnm1_complement.begin();
           ei != lnm1_eto_lnm1_complement.end(); ++ei) {
         const hier::Connector::NeighborSet& lnm1_nabrs = ei->second;
         hier::Connector::NeighborSet& ln_nabrs =
            lnm1_eto_ln_complement[ei->first];
         for (hier::Connector::NeighborSet::const_iterator na =
                 lnm1_nabrs.begin();
              na != lnm1_nabrs.end(); ++na) {
            hier::MappedBox tmp_mapped_box = *na;
            tmp_mapped_box.getBox().refine(getRatioToCoarserLevel(ln));
            tmp_mapped_box.getBox().grow(
               hier::IntVector(d_dim, getProperNestingBuffer(ln)));
            ln_nabrs.insert(ln_nabrs.end(), tmp_mapped_box);
         }
      }
      tmp_gcw *= getRatioToCoarserLevel(ln);
      hier::Connector lnm1_to_ln_complement;
      lnm1_to_ln_complement.swapInitialize(
         *hierarchy->getMappedBoxLevel(ln - 1),
         d_proper_nesting_complement[block_number][ln],
         zero_vector,
         lnm1_eto_ln_complement);
      lnm1_to_ln_complement.setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      hier::Connector ln_complement_to_lnm1(
         d_proper_nesting_complement[block_number][ln],
         *hierarchy->getMappedBoxLevel(ln - 1),
         zero_vector,
         d_from_nesting_complement[block_number][ln - 1].getNeighborhoodSets(),
         hier::MappedBoxLevel::DISTRIBUTED);
      ln_complement_to_lnm1.setConnectorType(hier::Connector::COMPLETE_OVERLAP);

      /*
       * 3. Bridge for Connector between level ln and d_proper_nesting_complement[ln].
       */
      oca.bridge(d_to_nesting_complement[block_number][ln],
         d_from_nesting_complement[block_number][ln],
         hierarchy->getConnector(ln, ln - 1),
         lnm1_to_ln_complement,
         ln_complement_to_lnm1,
         hierarchy->getConnector(ln - 1, ln),
         hierarchy->getRequiredConnectorWidth(ln - 1, ln));
   }
}

/*
 *************************************************************************
 * Make a map that can be used to grow a box within the domain where
 * necessary.
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::makeBoxGrowingMap(
   const hier::MappedBoxLevel& orig_mapped_box_level,
   const hier::Connector& hierarchy_to_orig,
   const hier::Connector& orig_to_hierarchy,
   const int ln,
   const hier::BlockId& block_id,
   hier::MappedBoxLevel& grown_mapped_box_level,
   hier::Connector& orig_to_grown,
   const hier::IntVector& min_size) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim,
      orig_mapped_box_level,
      grown_mapped_box_level,
      min_size);

   const hier::OverlapConnectorAlgorithm oca;

   const hier::IntVector& zero_vec = hier::IntVector::getZero(d_dim);
   const hier::IntVector& one_vec = hier::IntVector::getOne(d_dim);
   const int block_number = block_id.getBlockValue();

   const hier::MappedBoxSet& orig_mapped_boxes =
      orig_mapped_box_level.getMappedBoxes();
   const hier::Connector& hierarchy_to_nesting_complement =
      d_to_nesting_complement[block_number][ln];
   const hier::Connector& nesting_complement_to_hierarchy =
      d_from_nesting_complement[block_number][ln];

   hier::Connector orig_to_nesting_complement;
   oca.bridge(
      orig_to_nesting_complement,
      orig_to_hierarchy,
      hierarchy_to_nesting_complement,
      nesting_complement_to_hierarchy,
      hierarchy_to_orig);

   grown_mapped_box_level.initialize(
      orig_mapped_box_level.getRefinementRatio(),
      orig_mapped_box_level.getGridGeometry(),
      orig_mapped_box_level.getMPI());

   hier::NeighborhoodSet orig_eto_grown;

   /*
    * Above step ignored the domain complement components of nesting
    * definition (by necessity).  We must add the domain's complement
    * to neighbors in candidate_to_complement.
    */
   tbox::Pointer<hier::MappedBoxTree> refined_domain_complement_tree =
      d_domain_complement_tree[block_number].createRefinedTree(
         orig_mapped_box_level.getRefinementRatio());

   std::vector<hier::MappedBox> tmp_mapped_box_vector;
   tmp_mapped_box_vector.reserve(10);

   for (hier::MappedBoxSet::const_iterator ni = orig_mapped_boxes.begin();
        ni != orig_mapped_boxes.end(); ++ni) {
      const hier::MappedBox& omb = *ni;
      TBOX_ASSERT(!omb.isPeriodicImage());

      /*
       * Build the domain complement as seen in omb's locality,
       * consisting of omb's neighbors in orig_to_nesting_complement
       * and the domain's complement.
       */

      hier::BoxList local_nesting_complement;
      if (orig_to_nesting_complement.hasNeighborSet(omb.getId())) {
         orig_to_nesting_complement.getNeighborSet(omb.getId()).convertToBoxList(
            local_nesting_complement);
      }

      hier::Box box = omb.getBox();
      box.grow(one_vec);
      refined_domain_complement_tree->findOverlapMappedBoxes(
         tmp_mapped_box_vector,
         box);
      if (!tmp_mapped_box_vector.empty()) {
         hier::MappedBoxContainerUtils::convertMappedBoxVectorToBoxList(
               tmp_mapped_box_vector,
               local_nesting_complement);
         tmp_mapped_box_vector.clear();
      }

      hier::MappedBox grown_mapped_box = omb;
      hier::BoxUtilities::growBoxWithinDomain(
         grown_mapped_box.getBox(),
         local_nesting_complement,
         min_size);

      if (omb.getBox() != grown_mapped_box.getBox()) {
         grown_mapped_box_level.addMappedBox(grown_mapped_box);
         orig_eto_grown[omb.getId()].insert(grown_mapped_box);
      } else {
         grown_mapped_box_level.addMappedBox(omb);
      }

   }

   orig_to_grown.swapInitialize(
      orig_mapped_box_level,
      grown_mapped_box_level,
      zero_vec,
      orig_eto_grown);
   orig_to_grown.setConnectorType(hier::Connector::MAPPING);
}

/*
 *************************************************************************
 *
 * Set patch size and ghost cell information needed to create new
 * mesh levels.  The maximum number of ghost cells over all variables
 * is used to compute the smallest patch size allowed and the extent to
 * which patches may be extended to touch the physical boundary.  This
 * avoids problems in setting ghost cell data that may occur when ghost
 * cell regions intersect the physical boundary in strange ways.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::getGriddingParameters(
   hier::IntVector& smallest_patch,
   hier::IntVector& smallest_box_to_refine,
   hier::IntVector& largest_patch,
   hier::IntVector& extend_ghosts,
   const hier::PatchHierarchy& hierarchy,
   const int level_number,
   const bool for_building_finer) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(d_dim,
      smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts);

   TBOX_ASSERT((level_number >= 0) && (level_number < d_mb_hierarchy->getMaxNumberOfLevels()));

   /*
    * Determine maximum ghost cell width needed over all variables
    * currently known to the patch descriptor, and set the smallest
    * patch size.  The maximum number of ghosts is multiplied by the
    * error coarsen ratio (which should always be 1 unless regridding
    * uses error estimation).  This assures that when levels are
    * coarsened during error estimation, the coarser level patches
    * will meet the ghost cell constraint.
    */
   bool allow_patches_smaller_than_ghost_width =
      hierarchy.allowPatchesSmallerThanGhostWidth();
   smallest_patch = getSmallestPatchSize(level_number);
   hier::IntVector max_ghosts(
      hierarchy.getPatchDescriptor()->getMaxGhostWidth(d_dim));
   max_ghosts = max_ghosts * d_tag_init_strategy->getErrorCoarsenRatio();
   smallest_patch = getSmallestPatchSize(level_number);
   if (!allow_patches_smaller_than_ghost_width) {
      smallest_patch.max(max_ghosts);
   } else {
      const hier::IntVector& periodic_dirs(
         hierarchy.getGridGeometry()->getPeriodicShift(hier::IntVector::getOne(
               d_dim)));

      for (int i = 0; i < d_dim.getValue(); i++) {
         if (periodic_dirs(i)) {
            smallest_patch(i) =
               tbox::MathUtilities<int>::Max(smallest_patch(i), max_ghosts(i));
         }
      }
   }

   /*
    * Set largest patch size.
    */
   largest_patch = getLargestPatchSize(level_number);

   /*
    * Following if-check prevents changing a negative largest_patch
    * bacause TreeLoadBalancer interprets the non-negative value as
    * dissabling the upper limit on patch size.
    */
   if (largest_patch > hier::IntVector::getZero(d_dim)) {
      largest_patch.max(smallest_patch);
   }

   /*
    * Set the smallest box to refine based on the number of cells that
    * coarsened patches must accomodate to meet ghost cell needs of variables.
    * On the finest level, the smallest box to refine is the smallest patch.
    * On coarser levels, it is a function of the error coarsen ratio and
    * the ratio to the next finer level.
    *
    * If we are accessing gridding parameters for a level that is being
    * reconstructed, the smallest box to refine is not applicable so we
    * set it to -1 to indicate an invalid entry in case it is used.
    */
   if (for_building_finer) {

      smallest_box_to_refine = smallest_patch;

      smallest_box_to_refine /= getRatioToCoarserLevel(level_number);
      const hier::IntVector den(
         getRatioToCoarserLevel(level_number)
         / d_tag_init_strategy->getErrorCoarsenRatio());
      const hier::IntVector sz(hier::IntVector::ceiling(max_ghosts, den));
      smallest_box_to_refine.max(sz);

   } else {

      smallest_box_to_refine = hier::IntVector(d_dim, -1);

   }

   /*
    * Determine number of cells box may be extended to physical
    * domain boundary to accomodate ghost cells.
    */
   extend_ghosts = max_ghosts;

}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::warnIfDomainTooSmallInPeriodicDir(
   const hier::IntVector& smallest_patch_size,
   const hier::IntVector& domain_bounding_box_size) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim,
      smallest_patch_size,
      domain_bounding_box_size);

   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(d_dim);
   if (shift_catalog->isPeriodic()) {
      hier::IntVector min_nonzero_shift(d_dim, tbox::MathUtilities<int>::getMax());
      for (int n = 1; n < shift_catalog->getNumberOfShifts(); ++n) {
         hier::PeriodicId tmp(n);
         const hier::IntVector& shift(shift_catalog->shiftNumberToShiftDistance(tmp));
         for (int d = 0; d < d_dim.getValue(); ++d) {
            const int shift_mag = tbox::MathUtilities<int>::Abs(shift(d));
            if (min_nonzero_shift(d) < shift_mag) {
               min_nonzero_shift = hier::IntVector(d_dim, shift_mag);
            }
         }
      }
      TBOX_ASSERT(domain_bounding_box_size >= smallest_patch_size);
      if (domain_bounding_box_size == smallest_patch_size) {
         TBOX_ERROR("Domains should be bigger than smallest patch size in\n"
            << "any periodic index direction or bad things will happen\n"
            << "with Connector data and overlap computations.\n"
            << "Domain bounding box size = " << domain_bounding_box_size
            << "\nSmallest patch size = " << smallest_patch_size << "\n");
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::makeSingleBlockMappedBoxLevel(
   hier::MappedBoxLevel &singleblock_mapped_box_level,
   const hier::MappedBoxLevel &multiblock_mapped_box_level,
   const hier::BlockId &block_id) const
{
   const hier::MappedBoxSet &multiblock_mapped_boxes(multiblock_mapped_box_level.getGlobalMappedBoxes());
   hier::MappedBoxSet singleblock_mapped_boxes;

   for ( hier::MappedBoxSetSingleBlockIterator bi(multiblock_mapped_boxes,block_id);
         bi.isValid(); ++bi ) {
      singleblock_mapped_boxes.insert( singleblock_mapped_boxes.end(), *bi );
   }

   singleblock_mapped_box_level.initialize(
      singleblock_mapped_boxes,
      multiblock_mapped_box_level.getRefinementRatio(),
      multiblock_mapped_box_level.getGridGeometry(),
      multiblock_mapped_box_level.getMPI(),
      hier::MappedBoxLevel::GLOBALIZED );

   return;
}

/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::printClassData(
   std::ostream& os) const
{
   os << "\nMultiblockGriddingAlgorithm::printClassData..." << std::endl;
   os << "   static data members:" << std::endl;
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; d++) {
      os << "      (*s_tag_indx)[" << d << "] = "
      << (*s_tag_indx)[d] << std::endl;
      os << "      (*s_buf_tag_indx)[" << d << "] = "
      << (*s_buf_tag_indx)[d] << std::endl;
   }
   os << "MultiblockGriddingAlgorithm: this = "
      << (MultiblockGriddingAlgorithm *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_tag_init_strategy = "
      << (TagAndInitializeStrategy *)d_tag_init_strategy << std::endl;
   os << "d_box_generator = "
      << (BoxGeneratorStrategy *)d_box_generator << std::endl;
   os << "d_load_balancer = "
      << (LoadBalanceStrategy *)d_load_balancer << std::endl;
   os << "d_load_balancer0 = "
      << (LoadBalanceStrategy *)d_load_balancer0 << std::endl;
   os << "d_tag = " << d_tag.getPointer() << std::endl;
   os << "d_tag_indx = " << d_tag_indx << std::endl;
   os << "d_buf_tag_indx = " << d_buf_tag_indx << std::endl;
   os << "d_true_tag = " << d_true_tag << std::endl;
   os << "d_false_tag = " << d_false_tag << std::endl;
   os << "d_mb_hierarchy->getMaxNumberOfLevels() = " << d_mb_hierarchy->getMaxNumberOfLevels() << std::endl;

   int ln;

   os << "d_efficiency_tolerance..." << std::endl;
   for (ln = 0; ln < d_efficiency_tolerance.getSize(); ln++) {
      os << "    d_efficiency_tolerance[" << ln << "] = "
      << d_efficiency_tolerance[ln] << std::endl;
   }
   os << "d_combine_efficiency..." << std::endl;
   for (ln = 0; ln < d_combine_efficiency.getSize(); ln++) {
      os << "    d_combine_efficiency[" << ln << "] = "
      << d_combine_efficiency[ln] << std::endl;
   }

}

/*
 *************************************************************************
 *
 * Write out class version number and data members to database.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
   TBOX_ASSERT(!db.isNull());

   db->putInteger("ALGS_GRIDDING_ALGORITHM_VERSION",
      ALGS_GRIDDING_ALGORITHM_VERSION);

   db->putInteger("d_true_tag", d_true_tag);
   db->putInteger("d_false_tag", d_false_tag);
   db->putInteger("d_mb_hierarchy->getMaxNumberOfLevels()", d_mb_hierarchy->getMaxNumberOfLevels());

   db->putDoubleArray("d_efficiency_tolerance", d_efficiency_tolerance);
   db->putDoubleArray("d_combine_efficiency", d_combine_efficiency);

   db->putBool("d_sequentialize_patch_indices", d_sequentialize_patch_indices);
}

/*
 *************************************************************************
 *
 * If simulation is not from restart, read data from input database.
 * Otherwise, override data members initialized from restart with
 * values in the input database.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
   TBOX_ASSERT(!db.isNull());

   NULL_USE(is_from_restart);

   if (s_check_connectors == 0) {
      s_check_overflow_nesting = 'n';
      s_check_proper_nesting = 'n';
      s_check_connectors = 'n';
      s_print_mapped_box_level_hierarchy = 'n';
      s_print_steps = 'n';
   }

   s_check_overflow_nesting =
      db->getCharWithDefault("check_overflow_nesting", s_check_overflow_nesting);
   s_check_proper_nesting =
      db->getCharWithDefault("check_proper_nesting", s_check_proper_nesting);
   s_check_connectors =
      db->getCharWithDefault("check_connectors", s_check_connectors);
   s_print_mapped_box_level_hierarchy =
      db->getCharWithDefault("print_mapped_box_level_hierarchy",
         s_print_mapped_box_level_hierarchy);
   s_print_steps =
      db->getCharWithDefault("print_steps", s_print_steps);

   /*
    * Read input for efficiency tolerance.
    */

   if (db->keyExists("efficiency_tolerance")) {
      d_efficiency_tolerance = db->getDoubleArray("efficiency_tolerance");

      for (int ln = 0; ln < d_efficiency_tolerance.getSize(); ln++) {
         if ((d_efficiency_tolerance[ln] <= 0.0e0)
             || (d_efficiency_tolerance[ln] >= 1.0e0)) {
            TBOX_ERROR(
               d_object_name << ":  "
               <<
               "Key data `efficiency_tolerance' has values"
               << " out of range 0.0 < tol < 1.0.");

         }
      }

   }

   /*
    * Read input for combine efficiency.
    */

   if (db->keyExists("combine_efficiency")) {
      d_combine_efficiency = db->getDoubleArray("combine_efficiency");

      for (int ln = 0; ln < d_combine_efficiency.getSize(); ln++) {
         if ((d_combine_efficiency[ln] <= 0.0e0)
             || (d_combine_efficiency[ln] >= 1.0e0)) {
            TBOX_ERROR(
               d_object_name << ":  "
               <<
               "Key data `combine_efficiency' has values"
               << " out of range 0.0 < tol < 1.0.");

         }
      }

   }


   //d_proper_nesting_complement.resize(d_mb_hierarchy->getHierarchy(0)->getMaxNumberOfLevels(), hier::MappedBoxLevel(d_dim));
   d_proper_nesting_complement.resize(0);
   d_to_nesting_complement.resize(0);
   d_from_nesting_complement.resize(0);
   d_domain_complement_tree.resizeArray(0, hier::MappedBoxTree(d_dim));

   if (d_mb_hierarchy->getMaxNumberOfLevels() > 1) {
      tbox::Pointer<hier::PatchHierarchy> hierarchy(d_mb_hierarchy);
      tbox::Array<hier::IntVector> ratio_to_coarser(hierarchy->getMaxNumberOfLevels(),
                                                    hier::IntVector::getOne(d_dim));
      for ( int ln=0; ln<hierarchy->getMaxNumberOfLevels(); ++ln ) {
         ratio_to_coarser[ln] = hierarchy->getRatioToCoarserLevel(ln);
      }
      d_tag_init_strategy->checkCoarsenRatios(ratio_to_coarser);
   }

   std::string tmp_str;

   if ( db->isChar("check_nonrefined_tags") ) {
      // Temporary backward compatibility.
      tmp_str = db->getChar("check_nonrefined_tags");
   } else {
      tmp_str = db->getStringWithDefault("check_nonrefined_tags",
                                         std::string("WARN"));
   }
   d_check_nonrefined_tags = char(tolower(*tmp_str.c_str()));
   if (d_check_nonrefined_tags != 'i' &&
       d_check_nonrefined_tags != 'w' &&
       d_check_nonrefined_tags != 'e') {
      TBOX_ERROR("GriddingAlgorithm: input parameter check_nonrefined_tags\n"
         << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
   }

   if ( db->isChar("check_overlapping_patches") ) {
      // Temporary backward compatibility.
      tmp_str = db->getChar("check_overlapping_patches");
   } else {
      tmp_str = db->getStringWithDefault("check_overlapping_patches",
                                         std::string("IGNORE"));
   }
   d_check_overlapping_patches = char(tolower(*tmp_str.c_str()));
   if (d_check_overlapping_patches != 'i' &&
       d_check_overlapping_patches != 'w' &&
       d_check_overlapping_patches != 'e') {
      TBOX_ERROR(
         "GriddingAlgorithm: input parameter check_overlapping_patches\n"
         << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
   }

   d_sequentialize_patch_indices =
      db->getBoolWithDefault("sequentialize_patch_indices",
         d_sequentialize_patch_indices);

   d_enforce_proper_nesting =
      db->getBoolWithDefault("enforce_proper_nesting", d_enforce_proper_nesting);
   d_extend_to_domain_boundary =
      db->getBoolWithDefault("extend_to_domain_boundary",
         d_extend_to_domain_boundary);
   d_load_balance =
      db->getBoolWithDefault("load_balance", d_load_balance);

}

/*
 *************************************************************************
 *
 * Gets the database in the root database that corresponds to the object
 * name.  This method then checks to make sure that the version number
 * of the class is that same as the version number in the restart file.
 * If these values are equal, the data members are read in from the
 * restart database.
 *
 *************************************************************************
 */

void MultiblockGriddingAlgorithm::getFromRestart()
{
   tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if (root_db->isDatabase(d_object_name)) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file.");
   }

   int ver = db->getInteger("ALGS_GRIDDING_ALGORITHM_VERSION");
   if (ver != ALGS_GRIDDING_ALGORITHM_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
         <<
         "Restart file version different than class version.");
   }

   //d_proper_nesting_complement.resize(d_mb_hierarchy->getHierarchy(0)->getMaxNumberOfLevels(), hier::MappedBoxLevel(d_dim));
   d_proper_nesting_complement.resize(0);

   d_efficiency_tolerance = db->getDoubleArray("d_efficiency_tolerance");
   d_combine_efficiency = db->getDoubleArray("d_combine_efficiency");

   d_sequentialize_patch_indices = db->getBool("d_sequentialize_patch_indices");

}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::startupCallback()
{
   s_tag_indx = new tbox::Array<int>(tbox::Dimension::MAXIMUM_DIMENSION_VALUE, -1);
   s_buf_tag_indx = new tbox::Array<int>(tbox::Dimension::MAXIMUM_DIMENSION_VALUE, -1);
}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::shutdownCallback()
{
   delete s_tag_indx;
   delete s_buf_tag_indx;
}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::initializeCallback()
{
   /*
    * Timers:  for gathering performance information about box
    * calculus and other regridding operations.
    */
   t_load_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::load_balance");
   t_load_balance0 = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::load_balance0");
   t_load_balance_setup = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::load_balance_setup");
   t_bdry_fill_tags_create = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bdry_fill_tags_create");
   t_make_coarsest = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeCoarsestLevel()");
   t_make_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeFinerLevel()");
   t_make_finer_setup = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeFinerLevel()_setup");
   t_make_finer_tagging = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeFinerLevel()_tagging");
   t_make_finer_create = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeFinerLevel()_create");
   t_regrid_all_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::regridAllFinerLevels()");
   t_regrid_finer_create = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::regridFinerLevel()_create");
   t_fill_tags_from_mapped_box_level = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::fillTagsFromMappedBoxLevel()");
   t_tag_cells_for_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::tag_cells_for_refinement");
   t_buffer_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bufferTagsOnLevel()");
   t_second_finer_tagging = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::second_finer_tagging");
   t_bdry_fill_tags_comm = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bdry_fill_tags_comm");
   t_find_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::findRefinementBoxes()");
   t_find_boxes_containing_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::find_boxes_containing_tags");
   t_enforce_nesting = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::enforce_nesting");
   t_make_nesting_map = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeProperNestingMap()");
   t_make_nesting_map_compute = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeProperNestingMap()_compute");
   t_make_nesting_map_convert = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeProperNestingMap()_convert");
   t_use_nesting_map = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::use_nesting_map");
   t_make_overflow_map = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeOverflowNestingMap()");
   t_make_overflow_map_compute = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeOverflowNestingMap()_compute");
   t_make_overflow_map_convert = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeOverflowNestingMap()_convert");
   t_use_overflow_map = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::use_overflow_map");
   t_compute_external_parts = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::compute_external_parts");
   t_compute_nesting_violator = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::computeNestingViolator()");
   t_extend_to_domain_boundary = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::extend_to_domain_boundary");
   t_extend_within_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::extend_within_domain");
   t_grow_boxes_within_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::grow_boxes_within_domain");
   t_sort_nodes = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::sortNodes()");
   t_find_new_to_new = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::find_new_to_new");
   t_bridge_new_to_new = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bridge_new_to_new");
   t_bridge_new_to_coarser = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bridge_new_to_coarser");
   t_bridge_new_to_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bridge_new_to_finer");
   t_bridge_new_to_old = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bridge_new_to_old");
   t_bridge_links = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::bridge_links");
   t_modify_connector = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::modify_connector");
   t_make_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeCoarsestLevel()_make_domain");
   t_get_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::get_balance");
   t_use_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::use_balance");
   t_make_new = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::makeCoarsestLevel()_make_new");
   t_process_error = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::process_error");
   t_reset_hier = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::reset_hierarchy_config");
   t_misc1 = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::misc1");
   t_misc2 = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::misc2");
   t_misc3 = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::misc3");
   t_misc4 = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::misc4");
   t_misc5 = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::misc5");
   t_limit_overflow = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::limit_overflow");
   t_box_massage = tbox::TimerManager::getManager()->
      getTimer("mesh::MultiblockGriddingAlgorithm::box_massage");
}

/*
 *************************************************************************
 *************************************************************************
 */
void MultiblockGriddingAlgorithm::finalizeCallback()
{
   t_find_domain_complement.setNull();
   t_load_balance.setNull();
   t_load_balance0.setNull();
   t_load_balance_setup.setNull();
   t_bdry_fill_tags_create.setNull();
   t_make_coarsest.setNull();
   t_make_finer.setNull();
   t_make_finer_setup.setNull();
   t_make_finer_tagging.setNull();
   t_make_finer_create.setNull();
   t_regrid_all_finer.setNull();
   t_regrid_finer_create.setNull();
   t_bridge_links.setNull();
   t_fill_tags_from_mapped_box_level.setNull();
   t_tag_cells_for_refinement.setNull();
   t_buffer_tags.setNull();
   t_bdry_fill_tags_comm.setNull();
   t_second_finer_tagging.setNull();
   t_find_refinement.setNull();
   t_bridge_new_to_new.setNull();
   t_find_new_to_new.setNull();
   t_bridge_new_to_coarser.setNull();
   t_bridge_new_to_finer.setNull();
   t_bridge_new_to_old.setNull();
   t_find_boxes_containing_tags.setNull();
   t_enforce_nesting.setNull();
   t_make_nesting_map.setNull();
   t_make_nesting_map_compute.setNull();
   t_make_nesting_map_convert.setNull();
   t_use_nesting_map.setNull();
   t_make_overflow_map.setNull();
   t_make_overflow_map_compute.setNull();
   t_make_overflow_map_convert.setNull();
   t_use_overflow_map.setNull();
   t_compute_external_parts.setNull();
   t_compute_nesting_violator.setNull();
   t_extend_to_domain_boundary.setNull();
   t_extend_within_domain.setNull();
   t_grow_boxes_within_domain.setNull();
   t_sort_nodes.setNull();
   t_modify_connector.setNull();
   t_misc1.setNull();
   t_misc2.setNull();
   t_misc3.setNull();
   t_misc4.setNull();
   t_misc5.setNull();
   t_make_domain.setNull();
   t_get_balance.setNull();
   t_use_balance.setNull();
   t_make_new.setNull();
   t_process_error.setNull();
   t_limit_overflow.setNull();
   t_reset_hier.setNull();
   t_box_massage.setNull();
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
