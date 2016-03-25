//
// File:        MultiblockTester.C
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 1.6 $
// Modified:    $Date: 2004/02/11 23:46:08 $
// Description: Manager class for patch data communication tests.
//

#include "MultiblockTester.h"

#include "BergerRigoutsos.h"
#include "CoarsenOperator.h"
#include "PatchMultiblockTestStrategy.h"
#include "MultiblockPatchLevel.h"
#include "StandardTagAndInitialize.h"
#include "MultiblockGriddingAlgorithm.h"
#include "MultiblockGriddingTagger.h"
#include "RefineOperator.h"
#include "LoadBalancer.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

using namespace SAMRAI;

/*
*************************************************************************
*									*
* The constructor initializes object state.  The destructor is empty.   * 
*									*
*************************************************************************
*/

MultiblockTester::MultiblockTester(
   const string& object_name,
   tbox::Pointer<tbox::Database> main_input_db,
   PatchMultiblockTestStrategy* data_test,
   bool do_refine,
   bool do_coarsen,
   const string& refine_option)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!main_input_db.isNull());
   assert(data_test != (PatchMultiblockTestStrategy*)NULL);
#endif

   d_object_name = object_name;
   d_data_test_strategy = data_test;
   d_patch_hierarchy = NULL;

   d_fake_time = 0.0;

   d_is_reset = false;

   d_do_refine = do_refine;
   d_do_coarsen = false;
   if (!do_refine) { 
      d_do_coarsen = do_coarsen; 
   }

   d_refine_option = refine_option;
   if ( !( (d_refine_option == "INTERIOR_FROM_SAME_LEVEL")
          || (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL")) ) {
      TBOX_ERROR(object_name << " input error: illegal refine_option = "
                             << d_refine_option << endl);
   }

   d_refine_algorithm = new xfer::RefineAlgorithm<NDIM>();

   d_patch_data_components.clrAllFlags();
   d_refine_schedule.resizeArray(0);
   d_coarsen_schedule.resizeArray(0);

   d_source = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("SOURCE");
   d_destination = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("DESTINATION");
   d_refine_scratch = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("REFINE_SCRATCH");

   d_reset_source = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("SOURCE");
   d_reset_destination = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("DESTINATION");
   d_reset_refine_scratch = 
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("REFINE_SCRATCH");

   d_data_test_strategy->registerVariables(this);
}

MultiblockTester::~MultiblockTester()
{

}

/*
*************************************************************************
*                                                                       *
* Add variable with associated attributes to set of test variables.     *
*                                                                       *
*************************************************************************
*/

void MultiblockTester::registerVariable(
   const tbox::Pointer< hier::Variable<NDIM> > src_variable,
   const tbox::Pointer< hier::Variable<NDIM> > dst_variable,
   const hier::IntVector<NDIM>& src_ghosts,
   const hier::IntVector<NDIM>& dst_ghosts,
   const tbox::Pointer< xfer::Geometry<NDIM> > xfer_geom,
   const string& operator_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!src_variable.isNull());
   assert(!dst_variable.isNull());
   assert(!xfer_geom.isNull());
   assert(!operator_name.empty());
#endif

   hier::VariableDatabase<NDIM>* variable_db =
      hier::VariableDatabase<NDIM>::getDatabase();

   int src_id = variable_db->registerVariableAndContext(src_variable,
                                                        d_source,
			                                src_ghosts);

   int dst_id = variable_db->registerVariableAndContext(dst_variable,
		                                        d_destination,
					                dst_ghosts);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( src_id != -1 );
   assert( dst_id != -1 );
#endif

   d_patch_data_components.setFlag(src_id);
   d_patch_data_components.setFlag(dst_id);

   tbox::Pointer< xfer::RefineOperator<NDIM> > refine_operator = NULL;
   tbox::Pointer< xfer::CoarsenOperator<NDIM> > coarsen_operator = NULL;

   if (d_do_refine) {
      refine_operator = xfer_geom->lookupRefineOperator(src_variable,
                                                        operator_name);

      hier::IntVector<NDIM> scratch_ghosts =
         hier::IntVector<NDIM>::max(src_ghosts, dst_ghosts);
      scratch_ghosts.max(hier::IntVector<NDIM>(1));
      if (!refine_operator.isNull()) {
         scratch_ghosts.max(refine_operator->getStencilWidth());
      }
      int scratch_id = 
         variable_db->registerVariableAndContext(src_variable,
                                                 d_refine_scratch,
                                                 scratch_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( scratch_id != -1 );
#endif

      d_patch_data_components.setFlag(scratch_id);

      d_refine_algorithm->registerRefine(dst_id,
                                        src_id,
                                        scratch_id,
                                        refine_operator);

   } else if (d_do_coarsen) {
      coarsen_operator = xfer_geom->lookupCoarsenOperator(src_variable,
                                                          operator_name);
      d_coarsen_algorithm->registerCoarsen(dst_id,
                                          src_id,
                                          coarsen_operator);
   }

   registerVariableForReset(src_variable, dst_variable,
			    src_ghosts, dst_ghosts, xfer_geom,
                            operator_name);

   d_mblk_refine_alg =
      new mblk::MultiblockRefineAlgorithm<NDIM>(d_refine_algorithm,
                                                d_patch_hierarchy);
}

void MultiblockTester::registerVariableForReset(
   const tbox::Pointer< hier::Variable<NDIM> > src_variable,
   const tbox::Pointer< hier::Variable<NDIM> > dst_variable,
   const hier::IntVector<NDIM>& src_ghosts,
   const hier::IntVector<NDIM>& dst_ghosts,
   const tbox::Pointer< xfer::Geometry<NDIM> > xfer_geom,
   const string& operator_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!src_variable.isNull());
   assert(!dst_variable.isNull());
   assert(!xfer_geom.isNull());
   assert(!operator_name.empty());
#endif

   hier::VariableDatabase<NDIM>* variable_db =
      hier::VariableDatabase<NDIM>::getDatabase();

   int src_id = variable_db->registerVariableAndContext(src_variable,
                                                        d_reset_source,
			                                src_ghosts);

   int dst_id = variable_db->registerVariableAndContext(dst_variable,
		                                        d_reset_destination,
					                dst_ghosts);

   d_patch_data_components.setFlag(src_id);
   d_patch_data_components.setFlag(dst_id);

   tbox::Pointer< xfer::RefineOperator<NDIM> > refine_operator = NULL;
   tbox::Pointer< xfer::CoarsenOperator<NDIM> > coarsen_operator = NULL;

   if (d_do_refine) {
      refine_operator = xfer_geom->lookupRefineOperator(src_variable,
                                                        operator_name);

      hier::IntVector<NDIM> scratch_ghosts =
         hier::IntVector<NDIM>::max(src_ghosts, dst_ghosts);

      scratch_ghosts.max(hier::IntVector<NDIM>(1));
      if (!refine_operator.isNull()) {
         scratch_ghosts.max(refine_operator->getStencilWidth());
      }
      int scratch_id = 
         variable_db->registerVariableAndContext(src_variable,
                                                 d_reset_refine_scratch,
                                                 scratch_ghosts);

      d_patch_data_components.setFlag(scratch_id);

      d_reset_refine_algorithm.registerRefine(dst_id,
                                              src_id,
                                              scratch_id,
                                              refine_operator);

   } else if (d_do_coarsen) {
      coarsen_operator = xfer_geom->lookupCoarsenOperator(src_variable,
                                                          operator_name);
      d_reset_coarsen_algorithm.registerCoarsen(dst_id,
                                                src_id,
                                                coarsen_operator);
   }

}

/*
*************************************************************************
*                                                                       *
* Create refine and coarsen communication schedules for hierarchy.      * 
*                                                                       *
*************************************************************************
*/

void MultiblockTester::createRefineSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0)
           && (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif

   tbox::Pointer< mblk::MultiblockPatchLevel<NDIM> > level =
      d_patch_hierarchy->getPatchLevel(level_number);

   if (d_do_refine) {

      d_refine_schedule.resizeArray(
         d_patch_hierarchy->getFinestLevelNumber()+1);
      d_refine_schedule[level_number].setNull();

      if (   (level_number == 0)
          || (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") ) {
         d_refine_schedule[level_number] =
            d_mblk_refine_alg->createSchedule(level,
                                              level_number-1,
                                              d_patch_hierarchy,
                                              this);
      } else if (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") {
         d_refine_schedule[level_number] =
            d_mblk_refine_alg->createSchedule(level,
                                              NULL,
                                              level_number-1,
                                              d_patch_hierarchy,
                                              this);
      } 

   }

}

void MultiblockTester::resetRefineSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0)
           && (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif

   if (d_do_refine) {

      d_reset_refine_algorithm.resetSchedule(d_refine_schedule[level_number]);

   }

   d_is_reset = true;
}

void MultiblockTester::createCoarsenSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0)
           && (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif

   if (d_do_coarsen && (level_number > 0)) {

      d_coarsen_schedule.resizeArray(
         d_patch_hierarchy->getFinestLevelNumber()+1);
      d_coarsen_schedule[level_number].setNull();

      tbox::Pointer<hier::PatchLevel<NDIM> > level = 
         d_patch_hierarchy->getPatchLevel(level_number);
      tbox::Pointer<hier::PatchLevel<NDIM> > coarser_level = 
         d_patch_hierarchy->getPatchLevel(level_number-1);
      d_coarsen_schedule[level_number] =
         d_coarsen_algorithm->createSchedule(coarser_level, level, this);

   }

}

void MultiblockTester::resetCoarsenSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (level_number >= 0)
           && (level_number <= d_patch_hierarchy->getFinestLevelNumber()) );
#endif

   if (d_do_coarsen && (level_number > 0)) {

      d_reset_coarsen_algorithm.resetSchedule(
         d_coarsen_schedule[level_number]);

   }

   d_is_reset = true;
}


/*
*************************************************************************
*                                                                       *
* Perform data refine and coarsen operations.                           *
*                                                                       *
*************************************************************************
*/

void MultiblockTester::performRefineOperations(
   const int level_number)
{
   if (d_do_refine) {
      if (d_is_reset) {
         d_data_test_strategy->setDataContext(d_reset_refine_scratch);
      } else {
         d_data_test_strategy->setDataContext(d_destination);
      }
      if (!d_refine_schedule[level_number].isNull()) {
         d_refine_schedule[level_number]->fillData(d_fake_time);
      }
      d_data_test_strategy->clearDataContext();
   }
}

void MultiblockTester::performCoarsenOperations(
   const int level_number) 
{
   if (d_do_coarsen) {
      if (d_is_reset) {
         d_data_test_strategy->setDataContext(d_reset_source);
      } else {
         d_data_test_strategy->setDataContext(d_source);
      }
      if (!d_coarsen_schedule[level_number].isNull()) {
         d_coarsen_schedule[level_number]->coarsenData();
      }
      d_data_test_strategy->clearDataContext();
   }
}

/*
*************************************************************************
*                                                                       *
* Verify results of communication operations.                           *
*                                                                       *
*************************************************************************
*/

void MultiblockTester::verifyCommunicationResults() const
{
   if (d_is_reset) {
      d_data_test_strategy->setDataContext(d_reset_destination);
   } else {
      d_data_test_strategy->setDataContext(d_destination);
   }
   for (int ln = 0;
        ln <= d_patch_hierarchy->getFinestLevelNumber(); ln++) {
      tbox::Pointer< mblk::MultiblockPatchLevel<NDIM> > level =
         d_patch_hierarchy->getPatchLevel(ln);

      for (int nb = 0; nb < level->getNumberBlocks(); nb++) {
         tbox::Pointer<hier::PatchLevel<NDIM> > patch_level =
            level->getPatchLevelForBlock(nb);

         if (!patch_level.isNull()) {
            for (hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++) {
               tbox::Pointer<hier::Patch<NDIM> > patch =
                  patch_level->getPatch(p());

               d_data_test_strategy->verifyResults(*patch, d_patch_hierarchy,
                                                   ln, nb);
            }
         }
      }
   }
   d_data_test_strategy->clearDataContext();
}

/*
*************************************************************************
*                                                                       *
* Cell tagging and patch level data initialization routines declared    *
* in the GradientDetectorStrategy interface.  They are used to          *
* construct the hierarchy initially.                                    *
*                                                                       *
*************************************************************************
*/

void MultiblockTester::initializeLevelData(
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level,
   const bool allocate_data)
{
   (void) can_be_refined;
   (void) initial_time;
   (void) old_level;
   (void) allocate_data;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(hierarchy.isNull()));
   assert(!(hierarchy->getPatchLevel(level_number).isNull()));
   assert(level_number >= 0);
#endif

   tbox::Pointer< mblk::MultiblockPatchHierarchy<NDIM> > mblk_hierarchy =
      hierarchy;

   tbox::Pointer< mblk::MultiblockPatchLevel<NDIM> > level =
      hierarchy->getPatchLevel(level_number);

   level->allocatePatchData(d_patch_data_components, time);

   for (int nb = 0; nb < level->getNumberBlocks(); nb++) {
      tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
         mblk_hierarchy->getHierarchy(nb);
      tbox::Pointer<hier::PatchLevel<NDIM> > patch_level =
         level->getPatchLevelForBlock(nb);

      if (!patch_level.isNull()) {
         for (hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++) {
            tbox::Pointer<hier::Patch<NDIM> > patch =
               patch_level->getPatch(p());

            int level_num = patch_level->getLevelNumber();

            d_data_test_strategy->setDataContext(d_source);
            d_data_test_strategy->initializeDataOnPatch(*patch,
                                                        patch_hierarchy,
                                                        level_num, nb,
                                                        's');
            d_data_test_strategy->clearDataContext();

            d_data_test_strategy->setDataContext(d_reset_source);
            d_data_test_strategy->initializeDataOnPatch(*patch,
                                                        patch_hierarchy,
                                                        level_num, nb,
                                                        's');
            d_data_test_strategy->clearDataContext();

            if (d_do_coarsen) {

               d_data_test_strategy->setDataContext(d_destination);
               d_data_test_strategy->initializeDataOnPatch(*patch,
                                                           patch_hierarchy,
                                                           level_num, nb,
                                                           'd');
               d_data_test_strategy->clearDataContext(); 

               d_data_test_strategy->setDataContext(d_reset_destination);
               d_data_test_strategy->initializeDataOnPatch(*patch,
                                                           patch_hierarchy,
                                                           level_num, nb,
	         					   'd');
               d_data_test_strategy->clearDataContext(); 

            }

         }
      }
   }

}

void MultiblockTester::resetHierarchyConfiguration(
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy, 
   const int coarsest_level, 
   const int finest_level)
{
   (void) hierarchy; 
   (void) coarsest_level; 
   (void) finest_level;
}

void MultiblockTester::applyGradientDetector(
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double dt_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   (void) dt_time;
   (void) initial_time;
   (void) uses_richardson_extrapolation_too;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(hierarchy.isNull()));
   assert(!(hierarchy->getPatchLevel(level_number).isNull()));
#endif

   tbox::Pointer<hier::PatchLevel<NDIM> > level =
      hierarchy->getPatchLevel(level_number);

   d_data_test_strategy->setDataContext(d_source);

   for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
      tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

      d_data_test_strategy->tagCellsToRefine(*patch,
                                             hierarchy,
                                             level_number,
                                             tag_index);
   }

   d_data_test_strategy->clearDataContext();

}

/*
*************************************************************************
*                                                                       *
* Physical boundary condition and user-defined coarsen and refine       *
* operations declared in RefinePatchStrategy and CoarsenPatchStrategy.  *
* They are passed off to patch data test object.                        *
*                                                                       *
*************************************************************************
*/

void MultiblockTester::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch, 
   const double time,
   const hier::IntVector<NDIM>& gcw)
{
   (void) time;
   tbox::Pointer<hier::VariableContext> save_context =
      d_data_test_strategy->getDataContext();

   if (d_filling_coarse_scratch) {
      d_data_test_strategy->setDataContext(d_refine_scratch);
   } else {
      d_data_test_strategy->setDataContext(d_destination);
   }

   d_data_test_strategy->setPhysicalBoundaryConditions(patch,
                                                       d_fake_time,
                                                       gcw);

   d_data_test_strategy->setDataContext(save_context);

}

void MultiblockTester::fillSingularityBoundaryConditions(
   hier::Patch<NDIM>& patch,
   tbox::List<mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      singularity_patches,
   const double time,
   const hier::Box<NDIM>& fill_box,
   const hier::BoundaryBox<NDIM>& boundary_box)
{

   (void) time;
   tbox::Pointer<hier::VariableContext> save_context =
      d_data_test_strategy->getDataContext();

   if (d_filling_coarse_scratch) {
      d_data_test_strategy->setDataContext(d_refine_scratch);
   } else {
      d_data_test_strategy->setDataContext(d_destination);
   }

   d_data_test_strategy->fillSingularityBoundaryConditions(
      patch,
      singularity_patches,
      fill_box,
      boundary_box);

   d_data_test_strategy->setDataContext(save_context);
}

hier::IntVector<NDIM> MultiblockTester::getRefineOpStencilWidth() const
{
   return (hier::IntVector<NDIM>(1));
}

void MultiblockTester::preprocessRefine(
   hier::Patch<NDIM>& fine, 
   const hier::Patch<NDIM>& coarse, 
   const hier::Box<NDIM>& fine_box, 
   const hier::IntVector<NDIM>& ratio)
{
   d_data_test_strategy->preprocessRefine(fine, coarse, d_refine_scratch,
                                          fine_box, ratio);
}

void MultiblockTester::postprocessRefine(
   hier::Patch<NDIM>& fine, 
   const hier::Patch<NDIM>& coarse, 
   const hier::Box<NDIM>& fine_box, 
   const hier::IntVector<NDIM>& ratio)
{
   d_data_test_strategy->postprocessRefine(fine, coarse, d_refine_scratch,
                                           fine_box, ratio);
}

hier::IntVector<NDIM> MultiblockTester::getCoarsenOpStencilWidth() const
{
   return (hier::IntVector<NDIM>(0));
}

void MultiblockTester::preprocessCoarsen(
   hier::Patch<NDIM>& coarse, 
   const hier::Patch<NDIM>& fine, 
   const hier::Box<NDIM>& coarse_box, 
   const hier::IntVector<NDIM>& ratio)
{
   d_data_test_strategy->preprocessCoarsen(coarse, fine, NULL,
                                           coarse_box, ratio);
}

void MultiblockTester::postprocessCoarsen(
   hier::Patch<NDIM>& coarse, 
   const hier::Patch<NDIM>& fine, 
   const hier::Box<NDIM>& coarse_box, 
   const hier::IntVector<NDIM>& ratio)
{
   d_data_test_strategy->postprocessCoarsen(coarse, fine, NULL,
                                            coarse_box, ratio);
}

/*
*************************************************************************
*                                                                       *
* Create and configure gridding objects used to build the hierarchy.    *
* Then, create hierarchy and initialize data.  Note this routine        *
* must be called after variables are registered with this tester object.*
*                                                                       *
*************************************************************************
*/

void MultiblockTester::setupHierarchy(
   tbox::Pointer<tbox::Database> main_input_db,
   tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > cell_tagger)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!main_input_db.isNull());
#endif

   tbox::Pointer<tbox::Database> mult_db =
      main_input_db->getDatabase("Multiblock");

   int num_blocks = mult_db->getInteger("num_blocks");

   tbox::Array< tbox::Pointer<hier::PatchHierarchy<NDIM> > >
      hierarchies(num_blocks);

   char hier_name[32];
   for (int b = 0; b < num_blocks; b++) {
      sprintf(hier_name, "PatchHierarchy%d", b);

      hierarchies[b] =
         new hier::PatchHierarchy<NDIM>(hier_name,
                            d_data_test_strategy->getGridGeometry(b));
   }

   d_patch_hierarchy =
      new mblk::MultiblockPatchHierarchy<NDIM>("MultiblockPatchHierarchy",
                                               mult_db, hierarchies, true);

   tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator =
      new mesh::BergerRigoutsos<NDIM>();

   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
      new mesh::LoadBalancer<NDIM>(
         "LoadBalancer",
         main_input_db->getDatabase("LoadBalancer"));

   tbox::Pointer<mblk::MultiblockGriddingAlgorithm<NDIM> > gridding_alg =
      new mblk::MultiblockGriddingAlgorithm<NDIM>(
         "GriddingAlgorithm",
         main_input_db->getDatabase("GriddingAlgorithm"),
         d_patch_hierarchy,
         cell_tagger,
         box_generator,
         load_balancer,
         (mblk::MultiblockGriddingTagger<NDIM>*)NULL,
         true);

   int fake_tag_buffer = 0;

   gridding_alg->makeCoarsestLevel(d_patch_hierarchy, d_fake_time);

   bool initial_time = true;
   for (int ln = 0; gridding_alg->levelCanBeRefined(ln); ln++) {
      gridding_alg->makeFinerLevel(d_patch_hierarchy, d_fake_time,
                                   initial_time, fake_tag_buffer,
                                   d_fake_time);
   }

}


