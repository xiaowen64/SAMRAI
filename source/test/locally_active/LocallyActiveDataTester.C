//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/test/locally_active/LocallyActiveDataTester.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2195 $
// Modified:    $LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description: Class to test locally-active data communication
//

#include "LocallyActiveDataTester.h"

#include <stdlib.h>

#include "LocallyActiveVariableDatabase.h"
#include "LocallyActiveDataPatchLevelManager.h"

/*
 * Header files for SAMRAI library classes
 */
#include "BoundaryBox.h"
#include "BoxArray.h"
#include "BoxTree.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenOperator.h"
#include "GridGeometry.h"
#include "Index.h"
#include "LoadBalancer.h"
#include "ProcessorMapping.h"
#include "RefineOperator.h"
#include "CellDoubleConstantRefine.h"
#include "CartesianCellDoubleLinearRefine.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "tbox/TimerManager.h"
#include "tbox/InputDatabase.h"

/*
 * Declarations for fortran functions that operate on patch data.
 */
extern "C" {

void laplacian_(const int&, const int&,
                const int&, const int&,
#if (NDIM == 3)
                const int&, const int&,
#endif
                const int&, const int&,
#if (NDIM == 3)
                const int&, 
#endif
                const int&, const int&,
#if (NDIM == 3)
                const int&, 
#endif
                const double*,
                double*, double*);

void init_(const int&, const int&,
           const int&, const int&,
#if (NDIM == 3)
           const int&, const int&,
#endif
           const int&, const int&,
#if (NDIM == 3)
           const int&, 
#endif
           const double*, const double*,
           const double*, const double&,
           double*);

void bdry_(const int&, const int&,
           const int&, const int&,
#if (NDIM == 3)
           const int&, const int&,
#endif
           const int&, const int&,
#if (NDIM == 3)
           const int&, 
#endif
           const int&, const int&,
           const int&, const int&,
#if (NDIM == 3)
           const int&, const int&,
#endif
           const double*, const double*,
           const double*, const double&,
           double*);

}

#include "tbox/Array.C"
template class tbox::Array< LocallyActiveDataTester::Function >;

/*************************************************************************
 *
 * Constructor and Destructor for LocallyActiveDataTester class.
 *
 ************************************************************************/

LocallyActiveDataTester::LocallyActiveDataTester(const string& object_name,
                                 tbox::Pointer<tbox::Database> input_db,
                                 tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                 tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!hierarchy.isNull());
#endif

   d_object_name   = object_name;
   d_input_db = input_db;
   d_hierarchy  = hierarchy; 

   /*
    * Default values that may change via input file values.
    */
   d_num_test_functions = 1;
   d_laplacian_iterations = 1;

   d_nghosts_function = hier::IntVector<NDIM>(1);
   d_nghosts_solution = hier::IntVector<NDIM>(0);

   d_check_function_definitions = false;

   d_num_levels = 2;
   d_npatches_on_coarsest = hier::IntVector<NDIM>(1);

   d_check_results = true;
   d_function_distribution = UNDEFINED_FCN_DIST;
   d_function_specification = UNDEFINED_FCN_SPEC;

   getFromInput();

   tbox::TimerManager* tm = tbox::TimerManager::getManager();

   if (d_test_to_run == "LAPLACIAN_TEST") {

      d_refine_comm_setup_timers.resizeArray(d_num_levels); 
      d_coarsen_comm_setup_timers.resizeArray(d_num_levels); 
      d_refine_comm_execute_timers.resizeArray(d_num_levels); 
      d_coarsen_comm_execute_timers.resizeArray(d_num_levels);
      d_laplacian_execute_timers.resizeArray(d_num_levels);

      d_main_laplacian_timer = 
         tm->getTimer("apps::LocallyActiveDataTester::main_laplacian_timer");

      d_sum_refine_comm_setup_timer = 
         tm->getTimer("apps::LocallyActiveDataTester::sum_refine_comm_setup_timer");
      d_sum_refine_comm_execute_timer = 
         tm->getTimer("apps::LocallyActiveDataTester::sum_refine_comm_execute_timer");

   } else {

      if (d_test_to_run == "REFINE_TEST") {
         d_refine_comm_setup_timers.resizeArray(d_num_levels); 
         d_refine_comm_execute_timers.resizeArray(d_num_levels); 

         d_sum_refine_comm_setup_timer = 
            tm->getTimer("apps::LocallyActiveDataTester::sum_refine_comm_setup_timer");
         d_sum_refine_comm_execute_timer = 
            tm->getTimer("apps::LocallyActiveDataTester::sum_refine_comm_execute_timer");

      }

      if (d_test_to_run == "COARSEN_TEST") {
         d_coarsen_comm_setup_timers.resizeArray(d_num_levels);
         d_coarsen_comm_execute_timers.resizeArray(d_num_levels);
      }

   }

   for (int ln = 0; ln < d_num_levels; ln++) {
      string lev_string = tbox::Utilities::intToString(ln) + "\0";

      if (d_refine_comm_setup_timers.size() > ln) {
         d_refine_comm_setup_timers[ln] = 
            tm->getTimer(
               "apps::LocallyActiveDataTester::refine_comm_setup_" + lev_string);
      }
      if (d_refine_comm_execute_timers.size() > ln) {
         d_refine_comm_execute_timers[ln] = 
            tm->getTimer(
               "apps::LocallyActiveDataTester::refine_comm_execute_" + lev_string);
      }
      if (d_coarsen_comm_setup_timers.size() > ln) {
         d_coarsen_comm_setup_timers[ln] = 
            tm->getTimer(
               "apps::LocallyActiveDataTester::coarsen_comm_setup_" + lev_string);
      }
      if (d_coarsen_comm_execute_timers.size() > ln) {
         d_coarsen_comm_execute_timers[ln] = 
            tm->getTimer(
               "apps::LocallyActiveDataTester::coarsen_comm_execute_" + lev_string);
      }
      if (d_laplacian_execute_timers.size() > ln) {
         d_laplacian_execute_timers[ln] = 
            tm->getTimer(
               "apps::LocallyActiveDataTester::laplacian_execute_" + lev_string);
      }
   }

   registerVariables(visit_writer);

}

LocallyActiveDataTester::~LocallyActiveDataTester()
{
   d_bdry_fill_sched.resizeArray(0);
   d_coarsen_sched.resizeArray(0);
}

/*************************************************************************
 *
 * Create variables and register with variable database.
 *
 ************************************************************************/

void LocallyActiveDataTester::registerVariables(
   tbox::Pointer< appu::VisItDataWriter<NDIM> > visit_writer)
{
   hier::LocallyActiveVariableDatabase<NDIM>* var_db =
      hier::LocallyActiveVariableDatabase<NDIM>::getDatabase();

   for (int fn = 0; fn < d_functions.getSize(); fn++) {
      int depth = 1;

      string func_name = d_functions[fn].d_name;

      if (d_test_to_run != "COARSEN_TEST") {
         tbox::Pointer<hier::Variable<NDIM> > func = 
            new pdat::CellVariable<NDIM,double>(func_name, depth); 
         int func_id = var_db->registerVariable(func, d_nghosts_function);

         d_functions[fn].d_func = func;
         d_functions[fn].d_func_data_index = func_id;

         if (!visit_writer.isNull()) {
            visit_writer->registerPlotQuantity(func_name, "SCALAR", func_id);
         }
      }

      if (d_test_to_run != "REFINE_TEST") {
         string soln_name = func_name + ":soln";
         tbox::Pointer<hier::Variable<NDIM> > soln = 
            new pdat::CellVariable<NDIM,double>(soln_name, depth);
         int soln_id = var_db->registerVariable(soln, d_nghosts_solution);

         d_functions[fn].d_soln = soln;

         d_functions[fn].d_soln_data_index = soln_id;

         if (!visit_writer.isNull()) {
            visit_writer->registerPlotQuantity(soln_name, "SCALAR", soln_id);
         }
      }

   } // loop over functions
}

/*************************************************************************
 *
 * Establish patches on the hierarchy where the variables are 
 * defined, and allocate data for all mf variables.
 *
 ************************************************************************/

void LocallyActiveDataTester::setActivePatchesOnHierarchy()
{
   d_total_work_units.resizeArray(d_hierarchy->getNumberOfLevels());
   for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
      d_total_work_units[ln] = 0;
   }

   for (int fn = 0; fn < d_functions.getSize(); fn++) {
      for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
         tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

         tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
            hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
               getLocallyActiveDataPatchLevelManager(level);

         for (int ip = 0; ip < level->getNumberOfPatches(); ip++) {
            if (patchIsActiveForFunction(fn, level, ip)) {
               d_total_work_units[ln]++; 
               if (d_test_to_run != "COARSEN_TEST") {
                  level_mgr->setPatchDataActive(
                     hier::PatchDataId(d_functions[fn].d_func_data_index), 
                     hier::PatchNumber(ip) );
               }
               if (d_test_to_run != "REFINE_TEST") {
                  level_mgr->setPatchDataActive(
                     hier::PatchDataId(d_functions[fn].d_soln_data_index), 
                     hier::PatchNumber(ip) );
               }
            }

         } // loop over patches

      } // loop over levels

   } // loop over functions

}

bool LocallyActiveDataTester::patchIsActiveForFunction(
   int fn, tbox::Pointer<hier::PatchLevel<NDIM> > level, int ip)
{
   bool patch_is_active = false;

   if (level->getLevelNumber() == 0) {

      patch_is_active = true;

   } else {

      const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > ggeom =
         level->getGridGeometry();
      const double* lev0_xlo = ggeom->getXLower();
      const double* lev0_xup = ggeom->getXUpper();
      const double* lev0_dx  = ggeom->getDx();

      const hier::Box<NDIM>& pbox = level->getBoxForPatch(ip);
      const hier::IntVector<NDIM>& ratio_to_lev0 = level->getRatio();

      const hier::IntVector<NDIM> domain_per_shift = ggeom->getPeriodicShift(ratio_to_lev0);

      int idim;

      tbox::Array<double> centroid(NDIM);
      tbox::Array<double> domain_length(NDIM);
      for (idim = 0; idim < NDIM; idim++) {
         domain_length[idim] = lev0_xup[idim] - lev0_xlo[idim];
      }

      /*
       * Check possible periodic shifts of domain to determine active patches when
       * domain is periodic.
       */

      int shift_case = -1;
      bool do_shift = false;
      tbox::Array<bool> shift_in_dim(NDIM);

      while ( !patch_is_active && (shift_case <= NDIM*(NDIM-1)) ) {

         for (idim = 0; idim < NDIM; idim++) {
            shift_in_dim[idim] = false;
         }

         if (shift_case < 0) {

            do_shift = true;

         } else {

            if (shift_case < NDIM) {

               for (idim = 0; idim < NDIM; idim++) {
                  if (idim == shift_case) {
                     shift_in_dim[idim] = true;
                  } else {
                     shift_in_dim[idim] = false;
                  }
               }

            } else {

               switch (shift_case) {
                  case NDIM*(NDIM-1): {
                     for (idim = 0; idim < NDIM; idim++) {
                        shift_in_dim[idim] = true;
                     }
                     break;
                  }
#if (NDIM == 3)
                  case 3: {
                     shift_in_dim[0] = true;
                     shift_in_dim[1] = true;
                     shift_in_dim[2] = false;
                     break;
                  }
                  case 4: {
                     shift_in_dim[0] = true;
                     shift_in_dim[1] = false;
                     shift_in_dim[2] = true;
                     break;
                  }
                  case 5: {
                     shift_in_dim[0] = false;
                     shift_in_dim[1] = true;
                     shift_in_dim[2] = true;
                     break;
                  }
#endif
                  default: {
                  }
               } // switch (shift_case)

            } // else shift_case >= NDIM

            do_shift = true;
            for (idim = 0; idim < NDIM; idim++) {
               if (shift_in_dim[idim]) {
                  do_shift = do_shift && (domain_per_shift[idim] != 0);
               }
            }

         } // else shift_case >= 0

         if (do_shift) {
            bool test_patch_is_active = true; 
            for (idim = 0; idim < NDIM; idim++) {
               centroid[idim] = d_functions[fn].d_centroid[idim];
               double level_dx = lev0_dx[idim]/(double)ratio_to_lev0(idim);
               double p_xlo = lev0_xlo[idim] + level_dx * (double)(pbox.lower(idim));
               double p_xhi = lev0_xlo[idim] + level_dx * (double)(pbox.upper(idim)+1);
               if (shift_in_dim[idim]) {
                  if (centroid[idim] < p_xlo) {
                     centroid[idim] = centroid[idim] + domain_length[idim];
                  } else {
                     if (p_xhi < centroid[idim]) {
                        centroid[idim] = centroid[idim] - domain_length[idim];
                     }
                  }
               }
               double func_lo = centroid[idim] - d_functions[fn].d_radius;
               double func_hi = centroid[idim] + d_functions[fn].d_radius;
               test_patch_is_active = test_patch_is_active &&
                     ( tbox::MathUtilities<double>::Max(func_lo, p_xlo) <=
                       tbox::MathUtilities<double>::Min(func_hi, p_xhi) ); 
            }
            patch_is_active = test_patch_is_active; 
      
         } // if (do_shift)

         shift_case++;

      } // while !patch_is_active and possible shifts remain

   } // else level number > 0

   return (patch_is_active);
}

/*************************************************************************
 *
 * Initialize data for each multi-function variable.
 *
 ************************************************************************/
void LocallyActiveDataTester::initializeFunctionData()
{
   if (d_check_function_definitions) {
      initializeFunctionDataForVisitCheck(); 
   } else {
      for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
         tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

         tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
            hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
               getLocallyActiveDataPatchLevelManager(level);

         level_mgr->allocateAllPatchData();
      
         for (int fn = 0; fn < d_functions.getSize(); fn++) {
   
            int data_id = d_functions[fn].d_func_data_index; 
            if (d_test_to_run == "COARSEN_TEST") {
               data_id = d_functions[fn].d_soln_data_index;
            }
               
            hier::LocallyActiveDataPatchLevelManager<NDIM>::Iterator p =
               level_mgr->getIterator( hier::PatchDataId(data_id) );
            for ( ; p; p++) {

               tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

               if (d_test_to_run == "LAPLACIAN_TEST") {
                  setGaussianDataOnPatch(patch, fn);
               }

               if (d_test_to_run == "REFINE_TEST") {
                  setRefineTestDataOnPatch(patch, fn);
               }

               if (d_test_to_run == "COARSEN_TEST") {
                  setCoarsenTestDataOnPatch(patch, fn);
               }

            } // iterate over active patches

         } // iterate over functions

      } // iterate over levels

   }

}

void LocallyActiveDataTester::setGaussianDataOnPatch(
   tbox::Pointer<hier::Patch<NDIM> > patch, int fn)
{
   const hier::Box<NDIM>& pbox = patch->getBox();

   const hier::Index<NDIM>& ifirst = pbox.lower();
   const hier::Index<NDIM>& ilast  = pbox.upper();
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = 
      patch->getPatchGeometry();
   const double* dx = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();

   tbox::Pointer< pdat::CellData<NDIM,double> > func_data = 
      patch->getPatchData(d_functions[fn].d_func_data_index);

   func_data->fillAll(tbox::MathUtilities<double>::getSignalingNaN());

   init_(ifirst(0),ilast(0),
         ifirst(1),ilast(1),
#if (NDIM == 3)
         ifirst(2),ilast(2),      
#endif
         d_nghosts_function(0), d_nghosts_function(1),
#if (NDIM == 3)
         d_nghosts_function(2),      
#endif
         xlo,dx,
         d_functions[fn].d_centroid,
         d_functions[fn].d_alpha,
         func_data->getPointer());

   tbox::Pointer< pdat::CellData<NDIM,double> > soln_data = 
      patch->getPatchData(d_functions[fn].d_soln_data_index);
   soln_data->fillAll(tbox::MathUtilities<double>::getSignalingNaN());
}

void LocallyActiveDataTester::setRefineTestDataOnPatch(
   tbox::Pointer<hier::Patch<NDIM> > patch, int fn)
{
   int func_data_id = d_functions[fn].d_func_data_index;

   double func_init_val = getRefineTestValue(patch->getPatchLevelNumber(),
                                             patch->getPatchNumber(),
                                             func_data_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > func_data = 
      patch->getPatchData(func_data_id);
   func_data->fillAll(tbox::MathUtilities<double>::getSignalingNaN());
   func_data->fill(func_init_val, patch->getBox());
}

void LocallyActiveDataTester::setCoarsenTestDataOnPatch(
   tbox::Pointer<hier::Patch<NDIM> > patch, int fn)
{
   int soln_data_id = d_functions[fn].d_soln_data_index;

   double soln_init_val = getCoarsenTestValue(patch->getPatchLevelNumber(),
                                              patch->getPatchNumber(),
                                              soln_data_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > soln_data = 
      patch->getPatchData(soln_data_id);
   soln_data->fillAll(tbox::MathUtilities<double>::getSignalingNaN());
   soln_data->fillAll(soln_init_val, patch->getBox());
}

void LocallyActiveDataTester::initializeFunctionDataForVisitCheck()
{
   for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
         hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
            getLocallyActiveDataPatchLevelManager(level);

      for (int fn = 0; fn < d_functions.getSize(); fn++) {

         if (d_test_to_run != "COARSEN_TEST") {
            int f_id = d_functions[fn].d_func_data_index;
#if 0
            level_mgr->allocatePatchData(f_id);

            LocallyActiveDataPatchLevelManager::Iterator& ip =
               level_mgr->getIterator(f_id);
            for ( ; ip; ip++) {
#else
            level->allocatePatchData(f_id); 
            for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
#endif
               tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(ip());
               tbox::Pointer< pdat::CellData<NDIM,double> > func_data = patch->getPatchData(f_id);

               if ( level_mgr->getPatchDataActive( hier::PatchDataId(f_id), 
                                                   hier::PatchNumber(ip()) ) ) {
                  func_data->fillAll( static_cast<double>(fn) );
               } else {
                  func_data->fillAll( -1.0 );
               }
            }
         }

         if (d_test_to_run != "REFINE_TEST") {
            int s_id = d_functions[fn].d_soln_data_index;
#if 0
            level_mgr->allocatePatchData(s_id);

            LocallyActiveDataPatchLevelManager::Iterator& ip =
               level_mgr->getIterator(s_id);
            for ( ; ip; ip++) {
#else
            level->allocatePatchData(s_id); 
            for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
#endif
               tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(ip());
               tbox::Pointer< pdat::CellData<NDIM,double> > soln_data = patch->getPatchData(s_id);

               if ( level_mgr->getPatchDataActive( hier::PatchDataId(s_id), 
                                                   hier::PatchNumber(ip()) ) ) {
                  soln_data->fillAll( static_cast<double>(fn) );
               } else {
                  soln_data->fillAll( -1.0 );
               }
            }
         }

      } // iterate over functions

   } // iterate over levels

}

/*************************************************************************
 *
 * Setup communication schedules on MFHierarchy.  Note that this method
 * must be called only AFTER the hierarchy has been constructed and the
 * active data on each patch has been set. On each level, construct a 
 * refine schedule to fill ghost data, and a coarsen schedule to coarsen 
 * the finer level data to coarser level.
 *
 ************************************************************************/
void LocallyActiveDataTester::setupCommunication()
{
   tbox::Pointer< xfer::RefineOperator<NDIM> > refine_op;
   tbox::Pointer< xfer::CoarsenOperator<NDIM> > coarsen_op;

   if (d_test_to_run == "LAPLACIAN_TEST") {
      refine_op = new geom::CartesianCellDoubleLinearRefine<NDIM>();
      coarsen_op = new geom::CartesianCellDoubleWeightedAverage<NDIM>();
   }

   if (d_test_to_run == "REFINE_TEST") {
      refine_op = new pdat::CellDoubleConstantRefine<NDIM>();
   }

   if (d_test_to_run == "COARSEN_TEST") {
      coarsen_op = new geom::CartesianCellDoubleWeightedAverage<NDIM>();
   }

   for (int fn = 0; fn < d_functions.getSize(); fn++) {
      if (d_test_to_run == "LAPLACIAN_TEST") { 
         d_bdry_fill_alg.registerRefine(d_functions[fn].d_func_data_index,  // dst
                                        d_functions[fn].d_func_data_index,  // src
                                        d_functions[fn].d_func_data_index,  // scratch
                                        refine_op);
         d_coarsen_alg.registerCoarsen(d_functions[fn].d_soln_data_index,   // dst 
                                       d_functions[fn].d_soln_data_index,   // src
                                       coarsen_op);
      }
      if (d_test_to_run == "REFINE_TEST") {
         d_bdry_fill_alg.registerRefine(d_functions[fn].d_func_data_index,  // dst
                                        d_functions[fn].d_func_data_index,  // src
                                        d_functions[fn].d_func_data_index,  // scratch
                                        refine_op); 
      }
      if (d_test_to_run == "COARSEN_TEST") {
         d_coarsen_alg.registerCoarsen(d_functions[fn].d_soln_data_index,   // dst 
                                       d_functions[fn].d_soln_data_index,   // src
                                       coarsen_op);
      }
   }

   d_bdry_fill_sched.resizeArray(d_hierarchy->getNumberOfLevels());
   d_coarsen_sched.resizeArray(d_hierarchy->getNumberOfLevels());

   int finest_ln   = d_hierarchy->getFinestLevelNumber();
   int coarsest_ln = 0;

   for (int ln = finest_ln; ln >= coarsest_ln; ln--) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
         hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
            getLocallyActiveDataPatchLevelManager(level); 
      if (d_test_to_run != "COARSEN_TEST") {
         d_sum_refine_comm_setup_timer->start();
         d_refine_comm_setup_timers[ln]->start();
         d_bdry_fill_sched[ln] = 
            d_bdry_fill_alg.createSchedule(level, 
                                           level_mgr,
                                           ln-1,  
                                           d_hierarchy,
                                           this);
         d_refine_comm_setup_timers[ln]->stop();
         d_sum_refine_comm_setup_timer->stop();
      }
      if (ln > 0 && d_test_to_run != "REFINE_TEST") {
         tbox::Pointer< hier::PatchLevel<NDIM> > coarse_level = d_hierarchy->getPatchLevel(ln-1);
         tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > coarse_level_mgr =
            hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
               getLocallyActiveDataPatchLevelManager(coarse_level); 
         d_coarsen_comm_setup_timers[ln]->start();
         d_coarsen_sched[ln] = 
            d_coarsen_alg.createSchedule(coarse_level, coarse_level_mgr,
                                         level, level_mgr);  
         d_coarsen_comm_setup_timers[ln]->stop();
      }
   }
}

/*
************************************************************************
*
* Perform specified test.
*
************************************************************************
*/
void LocallyActiveDataTester::performTest()
{
   if (d_test_to_run == "LAPLACIAN_TEST") {

      computeLaplacians();

   } else {

      if (d_test_to_run == "REFINE_TEST") {

         d_sum_refine_comm_execute_timer->start();
         for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
            d_refine_comm_execute_timers[ln]->start();
            d_bdry_fill_sched[ln]->fillData(0.0);
            d_refine_comm_execute_timers[ln]->stop();
         }
         d_sum_refine_comm_execute_timer->stop();

      }

      if (d_test_to_run == "COARSEN_TEST") {

         for (int ln = d_hierarchy->getFinestLevelNumber(); ln > 0; ln--) {
            d_coarsen_comm_execute_timers[ln]->start();
            d_coarsen_sched[ln]->coarsenData();
            d_coarsen_comm_execute_timers[ln]->stop();
         }

      }

   }
}

/*
************************************************************************
*
* Compute laplacian of each function
*
************************************************************************
*/

void LocallyActiveDataTester::computeLaplacians()
{
   d_main_laplacian_timer->start();

   hier::LocallyActiveVariableDatabase<NDIM>* var_db =
      hier::LocallyActiveVariableDatabase<NDIM>::getDatabase();

   for (int lapcount = 0; lapcount < d_laplacian_iterations; lapcount++) {
      tbox::pout << "Laplacian iteration " << lapcount << "..." << endl;

      for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {

         tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

         tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
            var_db->getLocallyActiveDataPatchLevelManager(level);

         d_refine_comm_execute_timers[ln]->start();
         d_bdry_fill_sched[ln]->fillData(0.0);
         d_refine_comm_execute_timers[ln]->stop();

         d_laplacian_execute_timers[ln]->start();
         for (int fn = 0; fn < d_functions.getSize(); fn++) {

            int func_id = d_functions[fn].d_func_data_index;
            int soln_id = d_functions[fn].d_soln_data_index;

            hier::LocallyActiveDataPatchLevelManager<NDIM>::Iterator p =
               level_mgr->getIterator( hier::PatchDataId(func_id) );
            for ( ; p; p++) {
               tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());

               const hier::Index<NDIM> ifirst =patch->getBox().lower();
               const hier::Index<NDIM> ilast  =patch->getBox().upper();
               const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = 
                  patch->getPatchGeometry();
               const double* dx = patch_geom->getDx();

               tbox::Pointer< pdat::CellData<NDIM,double> > func_data = 
                  patch->getPatchData(func_id);
               tbox::Pointer< pdat::CellData<NDIM,double> > soln_data = 
                  patch->getPatchData(soln_id);

               laplacian_(ifirst(0),ilast(0),
                          ifirst(1),ilast(1),
#if (NDIM == 3)
                          ifirst(2),ilast(2),      
#endif
                          d_nghosts_function(0), d_nghosts_function(1),
#if (NDIM == 3)
                          d_nghosts_function(2),      
#endif
                          d_nghosts_solution(0), d_nghosts_solution(1),
#if (NDIM == 3)
                          d_nghosts_solution(2),      
#endif
                          dx,
                          func_data->getPointer(),
                          soln_data->getPointer());

            } // iterate over active patches for function
         
         } // loop over functions
         d_laplacian_execute_timers[ln]->stop();

         if (ln > 0) {
            d_coarsen_comm_execute_timers[ln]->start();
            d_coarsen_sched[ln]->coarsenData();
            d_coarsen_comm_execute_timers[ln]->stop();
         }

      } // loop over levels

   } // loop over lapcount
  
   d_main_laplacian_timer->stop(); 
}


/*************************************************************************
 *
 * Set boundary conditions for each function.  The different functions
 * may be distinguished by the "var_index" argument, which specifies
 * the variables descriptor index.
 *
 ************************************************************************/

void LocallyActiveDataTester::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const tbox::List<int>& scratch_data_ids,
   const double fill_time,
   const hier::IntVector<NDIM>& ghost_width_to_fill)
{

   const tbox::Pointer< geom::CartesianPatchGeometry<NDIM> > 
      patch_geom = patch.getPatchGeometry();

#if (NDIM == 3)
   const tbox::Array< hier::BoundaryBox<NDIM> > face_bdry =
      patch_geom->getCodimensionBoundaries(1);
   const int num_face_bdry_boxes = face_bdry.getSize();
#endif

   const tbox::Array< hier::BoundaryBox<NDIM> > edge_bdry =
      patch_geom->getCodimensionBoundaries(NDIM-1);
   const int num_edge_bdry_boxes = edge_bdry.getSize();

   const tbox::Array< hier::BoundaryBox<NDIM> > node_bdry =
      patch_geom->getCodimensionBoundaries(NDIM);
   const int num_node_bdry_boxes = node_bdry.getSize();

   const hier::Box<NDIM> interior(patch.getBox());

   for (tbox::List<int>::Iterator did(scratch_data_ids); did; did++) {

      int data_index = did();
     
      int function_id = -1;
      for (int fn = 0; fn < d_functions.getSize(); fn++) {
         if (d_functions[fn].d_func_data_index == data_index) {
            function_id = fn;
            break;
         }
      }
      if (function_id == -1) {
         TBOX_ERROR(d_object_name << "::setPhysicalBoundaryConditions():"
                    << "\nCatastrophic failure - no function ids match "
                    << "the supplied var_id." << endl);
      }

      hier::IntVector<NDIM> data_fill_gcw = 
         hier::IntVector<NDIM>::min(ghost_width_to_fill,
                        patch.getPatchData(data_index)->getGhostCellWidth());

#if (NDIM == 3)
      for (int ifbox = 0; ifbox < num_face_bdry_boxes; ifbox++ ) {
         hier::Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(face_bdry[ifbox],
                                                       interior,
                                                       data_fill_gcw);

         if (!fill_box.empty()) {
            setPatchBoundaryValues(patch,
                                   fill_box,
                                   function_id,
                                   data_index); 
         }

      } 
#endif

      for (int iebox = 0; iebox < num_edge_bdry_boxes; iebox++ ) {
         hier::Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(edge_bdry[iebox],
                                                       interior,
                                                       data_fill_gcw);

         if (!fill_box.empty()) {
            setPatchBoundaryValues(patch,
                                   fill_box,
                                   function_id,
                                   data_index);
         }

      }
   
      for (int inbox = 0; inbox < num_node_bdry_boxes; inbox++ ) {
         hier::Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(node_bdry[inbox],
                                                       interior,
                                                       data_fill_gcw);

         if (!fill_box.empty()) {
            setPatchBoundaryValues(patch,
                                   fill_box,
                                   function_id,
                                   data_index);
         }

      }

   }  // iterate over data ids to fill

}

void LocallyActiveDataTester::setPatchBoundaryValues(
   hier::Patch<NDIM>& patch,
   const hier::Box<NDIM>& fill_box,
   int function_id,
   int data_index)
{
   tbox::Pointer< pdat::CellData<NDIM,double> > var_data = patch.getPatchData(data_index);

   if (d_test_to_run == "LAPLACIAN_TEST") {

      tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

      const hier::Index<NDIM> ifirst = patch.getBox().lower();
      const hier::Index<NDIM> ilast  = patch.getBox().upper();
      const double* dx = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      const hier::Index<NDIM> bfirst =fill_box.lower();
      const hier::Index<NDIM> blast  =fill_box.upper();

      bdry_(ifirst(0),ilast(0),
            ifirst(1),ilast(1),
#if (NDIM == 3)
            ifirst(2),ilast(2),
#endif
            d_nghosts_function(0), d_nghosts_function(1),
#if (NDIM == 3)
            d_nghosts_function(2),
#endif
            bfirst(0),blast(0),
            bfirst(1),blast(1),
#if (NDIM == 3)
            bfirst(2),blast(2),
#endif
            xlo,dx,
            d_functions[function_id].d_centroid,
            d_functions[function_id].d_alpha,
            var_data->getPointer());

   } else {

      if (d_test_to_run == "REFINE_TEST") {
         double init_val = getRefineTestValue(patch.getPatchLevelNumber(),
                                              patch.getPatchNumber(),
                                              data_index);
         var_data->fill(init_val, fill_box);
      }

   }

}


/*************************************************************************
 *
 * Check values.
 *
 ************************************************************************/
int LocallyActiveDataTester::checkTestResult(ostream& os) const
{
   int fail_count = 0;
   if (d_check_results) {

      if (d_test_to_run == "REFINE_TEST") {
         fail_count += checkRefineTestResult(os);
      }
      if (d_test_to_run == "COARSEN_TEST") {
         fail_count += checkCoarsenTestResult(os);
      }
  
      if (d_test_to_run == "LAPLACIAN_TEST") {
         TBOX_ERROR("LocallyActiveDataTester::checkTestResult not implemented for "
                    << d_test_to_run << endl);
      }

   }
   return(fail_count); 
}

int LocallyActiveDataTester::checkCoarsenTestResult(ostream& os) const
{
   int fail_count = 0;

   bool all_correct_fine = true;
   bool all_correct_coarse = true;

   const hier::IntVector<NDIM>& refratio = d_ratio_to_coarser[1];
   const int coarse_ln = 0;
   const int fine_ln = 1;

   tbox::Pointer<hier::PatchLevel<NDIM> > coarse_level = d_hierarchy->getPatchLevel(coarse_ln);
   tbox::Pointer<hier::PatchLevel<NDIM> > fine_level = d_hierarchy->getPatchLevel(fine_ln);

   /*
    * Setup fine patch ids that intersect coarse patches.
    */

   hier::BoxArray<NDIM> refined_clev_boxes(coarse_level->getBoxes());
   refined_clev_boxes.refine(refratio);

   tbox::Pointer< hier::BoxTree<NDIM> > flev_boxtree = fine_level->getBoxTree();
   
   tbox::Array< tbox::Array<int> > fine_intersect_boxes(refined_clev_boxes.getNumberOfBoxes());

   for (int rcbi = 0; rcbi < refined_clev_boxes.getNumberOfBoxes(); rcbi++) {
      const hier::Box<NDIM>& rc_box = refined_clev_boxes[rcbi];
      flev_boxtree->findOverlapIndices(fine_intersect_boxes[rcbi], rc_box);
   }

   /*
    * Check patch data results.
    */

   if (tbox::SAMRAI_MPI::getRank() == 0) {
      tbox::pout << "\n\n***************************************************"
                 << "\nChecking results for fine level: fine_ln = " << fine_ln << endl;
   }

   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > fine_level_mgr =
      hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
         getLocallyActiveDataPatchLevelManager(fine_level);

   for (int fn = 0; fn < d_functions.getSize(); fn++) {
      int soln_id = d_functions[fn].d_soln_data_index;

      hier::LocallyActiveDataPatchLevelManager<NDIM>::Iterator p =
         fine_level_mgr->getIterator( hier::PatchDataId(soln_id) );
      for ( ; p; p++) {
         tbox::Pointer< hier::Patch<NDIM> > patch = fine_level->getPatch(p());
         const hier::Box<NDIM>& pbox = patch->getBox();
         tbox::Pointer< pdat::CellData<NDIM,double> > pdata = patch->getPatchData(soln_id);

         bool correct_on_patch = true;
         double correct_val = getCoarsenTestValue(fine_ln, p(), soln_id);
         for (pdat::CellIterator<NDIM> ic(pbox); ic; ic++) {
            pdat::CellIndex<NDIM> cell = ic();
            if ( !tbox::MathUtilities<double>::equalEps((*pdata)(cell), correct_val) ) {
               correct_on_patch = false;
            }
         }
         
         if (!correct_on_patch) {
            fail_count++; 
            all_correct_fine = false;
            tbox::perr << "\tCOARSEN_TEST FAILED: function " << fn
                       << " , patch " << p() << " , level " << fine_ln << endl; 
         }
         
      }  // iterate over fine patches on which function is defined
   }  // iterate over functions

   if (all_correct_fine) {
      tbox::pout << "\n\tAll results correct on level " << fine_ln << "..." << endl;
   }

   tbox::pout << "\n\n***************************************************"
        << "\nChecking results for coarse level: ln = " << coarse_ln << endl;

   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > coarse_level_mgr =
      hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
         getLocallyActiveDataPatchLevelManager(coarse_level);

   for (int fn = 0; fn < d_functions.getSize(); fn++) {
      int soln_id = d_functions[fn].d_soln_data_index;

      hier::LocallyActiveDataPatchLevelManager<NDIM>::Iterator p =
         coarse_level_mgr->getIterator( hier::PatchDataId(soln_id) );
      for ( ; p; p++) {
         tbox::Pointer< hier::Patch<NDIM> > patch = coarse_level->getPatch(p());
         const hier::Box<NDIM>& pbox = patch->getBox();
         hier::BoxList<NDIM> boxes_to_check(pbox);

         tbox::Pointer< pdat::CellData<NDIM,double> > pdata = patch->getPatchData(soln_id);

         for (int fib = 0; fib < fine_intersect_boxes[p()].size(); fib++) {
            int fine_patch_id = fine_intersect_boxes[p()][fib];
            if ( fine_level_mgr->getPatchDataActive( hier::PatchDataId(soln_id), 
                                                     hier::PatchNumber(fine_patch_id) ) ) {
               hier::Box<NDIM> cfbox(fine_level->getBoxForPatch(fine_patch_id));
               cfbox.coarsen(refratio);
               hier::Box<NDIM> check_box = pbox * cfbox; 
               double correct_val = getCoarsenTestValue(fine_ln, fine_patch_id, soln_id);

               bool correct_on_region = true;
               for (pdat::CellIterator<NDIM> ic(check_box); ic; ic++) {
                  pdat::CellIndex<NDIM> cell = ic();
                  if ( !tbox::MathUtilities<double>::equalEps((*pdata)(cell), correct_val) ) {
                     correct_on_region = false;
                  }
               }

               if (!correct_on_region) {
                  fail_count++;
                  all_correct_coarse = false;
                  tbox::perr << "COARSEN_TEST FAILED: function " << fn
                             << " , patch " << p() << " , level " << coarse_ln << endl; 
                  tbox::perr << "   coarse patch box = " << pbox << " , fine box region = "
                             << check_box << endl;
               }
         
               boxes_to_check.removeIntersections(check_box);
            }
         } // iterate over fine patch regions on coarse patch

         boxes_to_check.simplifyBoxes(); 

         if (boxes_to_check.getNumberOfBoxes() > 0) {
            double correct_val = getCoarsenTestValue(coarse_ln, p(), soln_id);
            for (hier::BoxList<NDIM>::Iterator lb(boxes_to_check); lb; lb++) {
               hier::Box<NDIM> check_box = pbox * lb(); 
               bool correct_on_region = true;
               for (pdat::CellIterator<NDIM> ic(check_box); ic; ic++) {
                  pdat::CellIndex<NDIM> cell = ic();
                  if ( !tbox::MathUtilities<double>::equalEps((*pdata)(cell), correct_val) ) {
                     correct_on_region = false;
                  }
               }

               if (!correct_on_region) {
                  fail_count++; 
                  all_correct_coarse = false;
                  tbox::perr << "COARSEN_TEST FAILED: function " << fn
                             << " , patch " << p() << " , level " << coarse_ln << endl; 
                  tbox::perr << "   coarse patch box = " << pbox << " , fine box region = "
                             << check_box << endl;
               }

            }
         }

      } // iterate over coarse patches on which function defined

   } // iterate over functions

   if (all_correct_coarse) {
      tbox::plog << "\n\tAll results correct on level " << coarse_ln << "..." << endl;
   }

   if (all_correct_coarse && all_correct_fine) {
      tbox::pout << "\n\n***************************************************" << endl;
      tbox::plog << "\n\tAll results correct..." << endl;
      tbox::pout << "\n\n***************************************************" << endl;
   } else {
      tbox::pout << "\n\n***************************************************" << endl;
   }

   return(fail_count);
}

int LocallyActiveDataTester::checkRefineTestResult(ostream& os) const
{
   int fail_count = 0;

   for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
      bool all_correct_on_physical_boundary = true;
      tbox::Array<bool> all_correct_from_level(ln+1);
      for (int j = 0; j < ln+1; j++) {
         all_correct_from_level[j] = true;
      }
      tbox::Array<int> global_check(tbox::SAMRAI_MPI::getNodes()); 
      int p, global_check_value;

      bool level_is_done = false;

      if (tbox::SAMRAI_MPI::getRank() == 0) {
         tbox::pout << "\n\n***************************************************"
              << "\n\tChecking refine results for level: ln = " << ln << endl;
      }

      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      const tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
         hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
         getLocallyActiveDataPatchLevelManager(level);
      tbox::Pointer< hier::BoxTree<NDIM> > level_boxtree = level->getBoxTree();
      const hier::BoxArray<NDIM>& level_boxes = level->getBoxes();
      const int nboxes = level_boxes.getNumberOfBoxes();

      /*
       *  Initialize arrays describing unchecked patches and patch regions.
       */

      tbox::Array< tbox::Array<hier::BoxList<NDIM> > > unchecked_boxes(nboxes);

      tbox::Array<bool> patch_is_done(nboxes);
      for (int j = 0; j < nboxes; j++) {
         patch_is_done[j] = true;
      }

      for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());

         hier::Box<NDIM> test_overlap_box = 
            hier::Box<NDIM>::grow(patch->getBox(), d_nghosts_function);
         for (int fn = 0; fn < d_functions.getSize(); fn++) {
            unchecked_boxes[ip()].resizeArray(d_functions.getSize());
            if ( level_mgr->getPatchDataActive( 
                            hier::PatchDataId(d_functions[fn].d_func_data_index), 
                            hier::PatchNumber(ip()) ) ) {
               unchecked_boxes[ip()][fn].unionBoxes(test_overlap_box);
               patch_is_done[ip()] = false;
            }
         }
      }

      /*
       *  Check results in physical boundary ghost regions.
       */

      if (tbox::SAMRAI_MPI::getRank() == 0) {
         tbox::pout << "\n\tChecking results in physical boundary ghost regions... " << endl;
      }

      for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(ip());
         if ( !patch_is_done[ip()] &&
              patch->getPatchGeometry()->getTouchesRegularBoundary() ) {
            for (int fn = 0; fn < d_functions.getSize(); fn++) {
               int data_id = d_functions[fn].d_func_data_index;
               if ( level_mgr->getPatchDataActive( hier::PatchDataId(data_id), 
                                                   hier::PatchNumber(ip()) ) ) {
                  double correct_val = getRefineTestValue(ln, ip(), data_id);
                  all_correct_on_physical_boundary =
                     ( all_correct_on_physical_boundary &&
                       checkRefineTestPhysicalBoundaries(patch, fn, correct_val,
                                                         unchecked_boxes[ip()][fn], os) );
               }
            }  // iterate over all functions
         }  // if patch is not done && patch touches physical boundary
      }  // iterate over all patches local to processor

      for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
         global_check[p] = 1;
      }
      global_check[tbox::SAMRAI_MPI::getRank()] = (all_correct_on_physical_boundary ? 1 : 0);
      tbox::SAMRAI_MPI::minReduction(global_check.getPointer(), tbox::SAMRAI_MPI::getNodes());

      if (tbox::SAMRAI_MPI::getRank() == 0) {
         global_check_value = 1;
         for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
            global_check_value = 
               tbox::MathUtilities<int>::Min(global_check_value, 
                                             global_check[p]);
         }
         if (global_check_value == 1) {
            tbox::pout << "\n\tAll results correct on physical boundary for level " 
                 << ln << " ..." << endl;
         } else {
            fail_count++;
            all_correct_on_physical_boundary = false;
            tbox::perr << "\n\tFAILED: - Incorrect results on physical boundary for level " 
                 << ln << " on processors: " << endl << "\t";
            for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
               if (global_check[p] == 0) {
                  tbox::perr << p << " , "; 
               }
            }
            tbox::perr << endl;
         }
      }

      level_is_done = checkLevelIsDone(level,
                                       unchecked_boxes,
                                       patch_is_done); 

      int check_ln = ln;

      while ( !level_is_done && (check_ln >= 0) ) {

         tbox::Pointer< hier::PatchLevel<NDIM> > src_level = d_hierarchy->getPatchLevel(check_ln);

         if (tbox::SAMRAI_MPI::getRank() == 0) { 
            tbox::pout << "\n\tChecking results on level " << ln
                 << " filled from level " << check_ln << "... " << endl;
         }

         all_correct_from_level[check_ln] = 
            checkRefineTestOnLevelInterior(level,  // check level (i.e., destination)
                                           src_level,  // source level
                                           unchecked_boxes,
                                           patch_is_done,
                                           os);

         for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
            global_check[p] = 1;
         }
         global_check[tbox::SAMRAI_MPI::getRank()] = (all_correct_from_level[check_ln] ? 1 : 0);
         tbox::SAMRAI_MPI::minReduction(global_check.getPointer(), tbox::SAMRAI_MPI::getNodes());

         if (tbox::SAMRAI_MPI::getRank() == 0) {
            global_check_value = 1;
            for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
               global_check_value = 
                  tbox::MathUtilities<int>::Min(global_check_value, 
                                                global_check[p]);
            }
            if (global_check_value == 1) {
               tbox::pout << "\n\tAll results correct on level " << ln 
                    << " filled from level " << check_ln << " ..." << endl;
            } else {
               fail_count++;
               all_correct_from_level[check_ln] = false;
               tbox::perr << "\n\tFAILED: - Incorrect results on level " << ln 
                    << " filled from level " << check_ln
                    << "  on processors (see log files): " << endl << "\t";
               for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
                  if (global_check[p] == 0) {
                     tbox::perr << p << " , ";
                  }
               }
               tbox::perr << endl;
            }
         }

         level_is_done = checkLevelIsDone(level,
                                          unchecked_boxes,
                                          patch_is_done); 

         check_ln--;

      }  // while !level_is_done && check_ln >= 0

      for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
         global_check[p] = 1;
      }
      global_check[tbox::SAMRAI_MPI::getRank()] = (level_is_done ? 1 : 0);
      tbox::SAMRAI_MPI::minReduction(global_check.getPointer(), tbox::SAMRAI_MPI::getNodes());

      if (tbox::SAMRAI_MPI::getRank() == 0) {
         global_check_value = 1;
         for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
            global_check_value = 
               tbox::MathUtilities<int>::Min(global_check_value, 
                                             global_check[p]);
         }
         if (global_check_value == 1) {
            tbox::pout << "\n\tLevel " << ln << " checking is done ..." << endl;
         } else {
            fail_count++;
            tbox::perr << "\n\tFAILED: - Result checking is done, but unchecked boxes"
                 << "  remain on processors (see log files): " << endl << "\t";
            for (p = 0; p < tbox::SAMRAI_MPI::getNodes(); p++) {
               if (global_check[p] == 0) {
                  tbox::perr << p << " , ";
               }
            }
            tbox::perr << endl;
         }
      }

      if (!level_is_done) {
         for (int ipatch = 0; ipatch < nboxes; ipatch++) {
            if (!patch_is_done[ipatch]) {
               tbox::plog << "\n Unchecked boxes for patch " << ipatch 
                    << " : " << level_boxes[ipatch] << endl;
               for (int fn = 0; fn < d_functions.getSize(); fn++) {
                  if ( !(unchecked_boxes[ipatch][fn].isEmpty()) ) {
                     tbox::plog << "   Function " << fn << ":" << endl;
                     for (hier::BoxList<NDIM>::Iterator lb(unchecked_boxes[ipatch][fn]); lb; lb++) {
                        tbox::plog << "      Box<NDIM> = " << lb() << endl;
                     }
                  }
               }
            } 
         }

      }

      if (tbox::SAMRAI_MPI::getRank() == 0) {
         bool level_test_successful = true;
      
         level_test_successful = level_test_successful && all_correct_on_physical_boundary;
         for (int j = 0; j < ln+1; j++) {
            level_test_successful = level_test_successful && all_correct_from_level[j];
         }

         if (level_test_successful) {
            tbox::pout << "\n\tAll results correct on level " << ln << "..." << endl;
            tbox::pout << "\n\n***************************************************" << endl;
         } else {
            tbox::perr << "\n\tFAILED: - Incorrect results encountered on level " 
                 << ln << " See log files..." << endl;
            tbox::pout << "\n\n***************************************************" << endl;
         }
      }

   } // iterate over levels

   return(fail_count);
}

bool LocallyActiveDataTester::checkRefineTestOnLevelInterior(
   tbox::Pointer<hier::PatchLevel<NDIM> > check_level,  
   tbox::Pointer<hier::PatchLevel<NDIM> > src_level,  
   tbox::Array< tbox::Array<hier::BoxList<NDIM> > >& unchecked_boxes,
   tbox::Array<bool>& patch_is_done,
   ostream& os) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!check_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(src_level->getLevelNumber() <= check_level->getLevelNumber());
   TBOX_ASSERT(unchecked_boxes.size() == check_level->getNumberOfPatches());
   TBOX_ASSERT(patch_is_done.size() == check_level->getNumberOfPatches());
#endif

   bool ret_all_correct = true;
  
   hier::IntVector<NDIM> refine_ratio = check_level->getRatio() / src_level->getRatio();

   const int src_ln   = src_level->getLevelNumber();

   hier::BoxArray<NDIM> refined_src_boxes(src_level->getBoxes());
   refined_src_boxes.refine(refine_ratio);
   tbox::Pointer< hier::BoxTree<NDIM> > src_level_boxtree = src_level->getBoxTree();

   const tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > src_level_mgr =
         hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->
         getLocallyActiveDataPatchLevelManager(src_level);

   for (hier::PatchLevel<NDIM>::Iterator ip(check_level); ip; ip++) {
      tbox::Pointer<hier::Patch<NDIM> > patch = check_level->getPatch(ip());

      if (!patch_is_done[ip()]) {

         hier::Box<NDIM> test_overlap_box(patch->getBox());
         test_overlap_box.coarsen(refine_ratio);
         test_overlap_box.grow(d_nghosts_function);

         tbox::Array<int> nabor_patches;
         src_level_boxtree->findOverlapIndices(nabor_patches, test_overlap_box);

         const int num_nabors = nabor_patches.getSize();
         for (int spp = 0; spp < num_nabors; spp++) {
            int sp = nabor_patches[spp];

            const hier::Box<NDIM>& nabor_box = refined_src_boxes[sp]; 
         
            tbox::List< hier::IntVector<NDIM> >::Iterator sh(src_level->getShiftsForPatch(sp));

            bool zero_shift = true;

            while (sh || zero_shift) {

               hier::IntVector<NDIM> shift(0);
               if (!zero_shift) {
                  shift = sh() * refine_ratio;
               }

               hier::Box<NDIM> shifted_nabor(hier::Box<NDIM>::shift(nabor_box, shift));

               for (int fn = 0; fn < d_functions.getSize(); fn++) {

                  hier::BoxList<NDIM> boxes_to_check = unchecked_boxes[ip()][fn];
                  if ( !boxes_to_check.isEmpty() ) {

                     const int data_id = d_functions[fn].d_func_data_index;

                     if ( src_level_mgr->getPatchDataActive( hier::PatchDataId(data_id), 
                                                             hier::PatchNumber(sp) ) ) {
                        double correct_val =
                               getRefineTestValue(src_ln, sp, data_id);
                        for (hier::BoxList<NDIM>::Iterator cbox(boxes_to_check); cbox; cbox++) {
                           hier::Box<NDIM> intersection = shifted_nabor * cbox();
                           ret_all_correct =
                               ( ret_all_correct &&
                                 checkRefineTestOnBox(patch, fn, intersection,
                                                      correct_val, os) );
                           unchecked_boxes[ip()][fn].removeIntersections(intersection);
                        }
                     }
                  }  // check values from shifted source patch on current patch 

               } // iterate over functions

               if (!zero_shift) {
                  sh++;
               } else {
                  zero_shift = false;
               }

            } // while loop (sh || zero_shift)

         } // loop over potential source patches

      }  // if patch is not done

   }  // iterate over all patches local to processor

   return(ret_all_correct);
}

bool LocallyActiveDataTester::checkLevelIsDone(
   tbox::Pointer<hier::PatchLevel<NDIM> > level,
   tbox::Array< tbox::Array<hier::BoxList<NDIM> > >& unchecked_boxes,
   tbox::Array<bool>& patch_is_done) const
{
   bool ret_level_done = true;
   for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
      if (!patch_is_done[ip()]) {
         bool now_patch_is_done = true;
         for (int fn = 0; fn < d_functions.getSize(); fn++) {
            now_patch_is_done = ( now_patch_is_done &&
                                  unchecked_boxes[ip()][fn].isEmpty() );
         }
         patch_is_done[ip()] = now_patch_is_done;
         ret_level_done = ret_level_done && now_patch_is_done;
      }
   }  // iterate over all patches local to processor

   return(ret_level_done);
}

bool LocallyActiveDataTester::checkRefineTestPhysicalBoundaries(
   tbox::Pointer<hier::Patch<NDIM> > patch,
   int function_id,
   double correct_value,
   hier::BoxList<NDIM>& check_boxes_to_modify,
   ostream& os) const
{
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
   const hier::Box<NDIM> interior(patch->getBox());
  
   bool all_good = true;

#if (NDIM == 3)
   const tbox::Array<hier::BoundaryBox<NDIM> > face_bdry =
      patch_geom->getCodimensionBoundaries(1);
   const int num_face_bdry_boxes = face_bdry.getSize();
   for (int ifbox = 0; ifbox < num_face_bdry_boxes; ifbox++ ) {
      hier::Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(face_bdry[ifbox],
                                                    interior,
                                                    d_nghosts_function);
      all_good = ( all_good && 
                   checkRefineTestOnBox(patch, function_id, fill_box, 
                                        correct_value, os) ); 
      check_boxes_to_modify.removeIntersections(fill_box);
   }
#endif

   const tbox::Array<hier::BoundaryBox<NDIM> > edge_bdry =
      patch_geom->getCodimensionBoundaries(NDIM-1);
   const int num_edge_bdry_boxes = edge_bdry.getSize();
   for (int iebox = 0; iebox < num_edge_bdry_boxes; iebox++ ) {
      hier::Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(edge_bdry[iebox],
                                                    interior,
                                                    d_nghosts_function);
      all_good = ( all_good && 
                   checkRefineTestOnBox(patch, function_id, fill_box,
                                        correct_value, os) ); 
      check_boxes_to_modify.removeIntersections(fill_box);
   }

   const tbox::Array<hier::BoundaryBox<NDIM> > node_bdry =
      patch_geom->getCodimensionBoundaries(NDIM);
   const int num_node_bdry_boxes = node_bdry.getSize();
   for (int inbox = 0; inbox < num_node_bdry_boxes; inbox++ ) {
      hier::Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(node_bdry[inbox],
                                                    interior,
                                                    d_nghosts_function);
      all_good = ( all_good && 
                   checkRefineTestOnBox(patch, function_id, fill_box,
                                        correct_value, os) ); 
      check_boxes_to_modify.removeIntersections(fill_box);
   }

   return(all_good);
}

bool LocallyActiveDataTester::checkRefineTestOnBox(
   tbox::Pointer<hier::Patch<NDIM> > patch,
   int function_id,
   const hier::Box<NDIM>& region,
   double correct_value,
   ostream& os) const
{
   tbox::Pointer< pdat::CellData<NDIM,double> > pdata = 
      patch->getPatchData(d_functions[function_id].d_func_data_index);
   bool all_good = true;
   for (pdat::CellIterator<NDIM> ic(region); ic; ic++) {
      pdat::CellIndex<NDIM> cell = ic();
      if ( !tbox::MathUtilities<double>::equalEps((*pdata)(cell), correct_value) ) {
         all_good = false;
      }
   }
   if (!all_good) {
      os << "REFINE_TEST FAILED: function " << function_id
         << " , patch " << patch->getPatchNumber() 
         << " , level " << patch->getPatchLevelNumber() << endl;
      os << "    patch box = " << patch->getBox() << " , box region = "
         << region << endl;
   }

   return(all_good);
}


/*************************************************************************
 *
 * Specify values to set on a patch, based on level number, function id,
 * and patch number.
 *
 ************************************************************************/

double LocallyActiveDataTester::getRefineTestValue(int level_number,
                                                   int patch_number,
                                                   int data_id) const
{
   double ret_value = (double) data_id;
   
   switch(d_function_specification) {

      case CONSTANT_TEST: {
         ret_value  = (double)level_number +
                      0.1*(double)data_id +
                      0.00001*(double)patch_number;
         break;
      }

      case CONSTANT_PERF: {
         ret_value = (double)data_id; 
         break;
      }

      default: {
         TBOX_ERROR("LocallyActiveDataTester::getRefineTestValue error ...\n"
                    << "function_spec case = " << d_function_specification
                    << " not implemented for refine test." << endl);

      }

   } //  switch(d_function_specification)

   return(ret_value);
}

double LocallyActiveDataTester::getCoarsenTestValue(int level_number,
                                            int patch_number,
                                            int data_id) const
{
   double ret_value = (double) data_id;

   switch(d_function_specification) {

      case CONSTANT_TEST: {
         ret_value  = (double)level_number +
                      0.1*(double)data_id +
                      0.00001*(double)patch_number;
         break;
      }

      case CONSTANT_PERF: {
         ret_value = (double)data_id;
         break;
      }

      default: {
         TBOX_ERROR("LocallyActiveDataTester::getCoarsenTestValue error ...\n"
                    << "function_spec case = " << d_function_specification
                    << " not implemented for coarsen test." << endl);
      }

   } //  switch(d_function_specification)

   return(ret_value);
}

void LocallyActiveDataTester::buildPatchHierarchy()
{
   tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
   const hier::BoxArray<NDIM> domain_boxes = grid_geom->getPhysicalDomain();

   tbox::Pointer<tbox::Database> lb_input_db = new tbox::InputDatabase("LoadBalancer");
   tbox::Array<int> num_proc(NDIM);
   for (int idim = 0; idim < NDIM; idim++) {
      num_proc[idim] = d_npatches_on_coarsest(idim);
   }
   lb_input_db->putIntegerArray("processor_layout", num_proc);
   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = new mesh::LoadBalancer<NDIM>(lb_input_db);

   hier::BoxList<NDIM> in0_boxes(domain_boxes);
   hier::BoxArray<NDIM> level0_boxes;
   hier::ProcessorMapping level0_mapping;
   int level_num = 0;
   hier::IntVector<NDIM> cut_factor(d_nghosts_function);
   hier::IntVector<NDIM> bad_interval(d_nghosts_function);
 
   hier::Box<NDIM> domain_box = in0_boxes.getBoundingBox();
   hier::IntVector<NDIM> box_size = domain_box.upper() - domain_box.lower() + hier::IntVector<NDIM>(1);
   hier::IntVector<NDIM> largest_patch_size = box_size/d_npatches_on_coarsest;
   hier::IntVector<NDIM> smallest_patch_size = largest_patch_size - hier::IntVector<NDIM>(1);

   tbox::plog << "domain_box = " << domain_box << endl;
   tbox::plog << "d_npatches_on_coarsest = " << d_npatches_on_coarsest << endl;
   tbox::plog << "processor_layout = ";
   for (int idim = 0; idim < NDIM; idim++) {
      tbox::plog << num_proc[idim] << " , ";
   }
   tbox::plog << endl;
   tbox::plog << "cut_factor = " << cut_factor << endl;
   tbox::plog << "bad_interval = " << bad_interval << endl;
   tbox::plog << "box_size = " << box_size << endl;
   tbox::plog << "largest_patch_size = " << largest_patch_size << endl;
   tbox::plog << "smallest_patch_size = " << smallest_patch_size << endl;
 
   load_balancer->loadBalanceBoxes(level0_boxes, level0_mapping,
                                   in0_boxes, d_hierarchy, level_num,
                                   grid_geom->getPhysicalDomain(),
                                   hier::IntVector<NDIM>(1),
                                   smallest_patch_size,
                                   largest_patch_size,
                                   cut_factor, bad_interval);
 
   tbox::plog << "\n\nPatch-Processor Mapping - Level 0" << endl;
   tbox::Array<int> mapping = level0_mapping.getProcessorMapping();
   for (int i = 0; i < mapping.getSize(); i++) {
      tbox::plog << "  Patch<NDIM>: " << i << " : " << level0_boxes[i] 
           << " : " << mapping[i] << endl;
   }

   hier::IntVector<NDIM> ratio_to_level_zero(1);
 
   d_hierarchy->makeNewPatchLevel(level_num, ratio_to_level_zero,
                                  level0_boxes, level0_mapping);
 
   for (int ln = 1; ln < d_num_levels; ln++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln-1);
      hier::BoxArray<NDIM> coarser_level_boxes = coarser_level->getBoxes();
 
      ratio_to_level_zero = coarser_level->getRatio() * d_ratio_to_coarser[ln];
 
      hier::BoxList<NDIM> new_level_domain(coarser_level_boxes);
      new_level_domain.refine(d_ratio_to_coarser[ln]);
 
      hier::BoxArray<NDIM> new_level_boxes;
      hier::ProcessorMapping new_level_mapping;

      hier::BoxArray<NDIM> physical_domain;
      grid_geom->computePhysicalDomain(physical_domain,
                                       ratio_to_level_zero);

      cut_factor = d_ratio_to_coarser[ln];
 
      load_balancer->loadBalanceBoxes(new_level_boxes, new_level_mapping,
                                      new_level_domain, d_hierarchy, ln,
                                      physical_domain,
                                      ratio_to_level_zero,
                                      smallest_patch_size,
                                      largest_patch_size,
                                      cut_factor, bad_interval);
 
      tbox::plog << "\n\nPatch-Processor Mapping - Level " << ln << endl;
      tbox::Array<int> mapping = new_level_mapping.getProcessorMapping();
      for (int i = 0; i < mapping.getSize(); i++) {
         tbox::plog << "  Patch<NDIM>: " << i << " : " << new_level_boxes[i] 
              << " : " << mapping[i] << endl;
      }
 
      d_hierarchy->makeNewPatchLevel(ln, ratio_to_level_zero,
                                     new_level_boxes, new_level_mapping);
   }

}

/*************************************************************************
 *
 * Read LocallyActiveDataTester data from input file.
 *
 ************************************************************************/
void LocallyActiveDataTester::getFromInput() 
{
   if (d_input_db->keyExists("test_to_run")) {
      d_test_to_run = d_input_db->getString("test_to_run");
      if ( (d_test_to_run != "REFINE_TEST") &&
           (d_test_to_run != "COARSEN_TEST") &&
           (d_test_to_run != "LAPLACIAN_TEST") ) {
         TBOX_ERROR("Input error in " << d_object_name << ": "
                    << "'test_to_run' input string must be either "
                    << "'REFINE_TEST', 'COARSEN_TEST', or 'LAPLACIAN_TEST'." << endl);
      }
   } else {
      TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nNo 'test_to_run' item found!" << endl);
   }

   if (d_test_to_run == "LAPLACIAN_TEST") {
      int default_laplacian_iterations = d_laplacian_iterations;
      d_laplacian_iterations = 
          d_input_db->getIntegerWithDefault("laplacian_iterations",
                                            default_laplacian_iterations);
      if (d_laplacian_iterations <= 0) {
         TBOX_WARNING("Input warning in " << d_object_name << ":"
                      << "\nUsing zero for 'laplacian_iterations'" << endl);
      }
   }

   bool default_check_results = d_check_results; 
   d_check_results =
      d_input_db->getBoolWithDefault("check_results", default_check_results);

   bool default_check_function_definitions = d_check_function_definitions;
   d_check_function_definitions =
      d_input_db->getBoolWithDefault("check_function_definitions",
                                     default_check_function_definitions);   

   if (d_input_db->keyExists("nghosts_function")) {
      int* temp_nghosts = d_nghosts_function;
      d_input_db->getIntegerArray("nghosts_function", temp_nghosts, NDIM);
      if (d_nghosts_function.min() < 1) {
         TBOX_WARNING("Input warning for " << d_object_name << ":  "
                      << "Key data 'nghosts_function' "
                      << "has entry < 1; using default value." << endl);
      }
   }

   if (d_input_db->keyExists("nghosts_solution")) {
      int* temp_nghosts = d_nghosts_solution;
      d_input_db->getIntegerArray("nghosts_solution", temp_nghosts, NDIM);
      if (d_nghosts_solution.min() < 0) {
         TBOX_WARNING("Input warning for " << d_object_name << ":  "
                      << "Key data 'nghosts_solution' "
                      << "has entry < 0; using default value." << endl);
         }
   }

   if (d_input_db->keyExists("GriddingParameters")) {
      getGriddingParametersFromInput(
         d_input_db->getDatabase("GriddingParameters")); 
   } else {
      TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nNo 'GriddingParameters' input found!" << endl);
   }

   if (d_input_db->keyExists("function_distribution")) {
      string func_dist_string = d_input_db->getString("function_distribution");
      if ( (func_dist_string != "EXPLICIT_FUNCTIONS") &&
           (func_dist_string != "UNIFORM_FUNCTIONS") &&
           (func_dist_string != "RANDOM_FUNCTIONS") &&
           (func_dist_string != "FINE_PATCH_ALL") &&
           (func_dist_string != "FINE_PATCH_MOD2") &&
           (func_dist_string != "FINE_PATCH_EVEN_ODD") ) {
         TBOX_ERROR("Input error in " << d_object_name << ": "
                    << "'function_distribution' input string given as " << func_dist_string << ".\n"
                    << " Input value must be either "
                    << "'EXPLICIT_FUNCTIONS', 'UNIFORM_FUNCTIONS', 'RANDOM_FUNCTIONS'\n"
                    << "'FINE_PATCH_ALL', 'FINE_PATCH_MOD2', or 'FINE_PATCH_EVEN_ODD'." << endl);
      }

      if (func_dist_string == "EXPLICIT_FUNCTIONS") {
         d_function_distribution = EXPLICIT_FUNCTIONS; 
         if (d_input_db->isDatabase("ExplicitFunctionData")) {
            getExplicitFunctionDataFromInput(
               d_input_db->getDatabase("ExplicitFunctionData"));
         } else {
            TBOX_ERROR("Input error in " << d_object_name
                       << "\nNo 'ExplicitFunctionData' input database found for test " 
                       << d_test_to_run << endl);
         }
      }

      if (func_dist_string == "UNIFORM_FUNCTIONS") {
         d_function_distribution = UNIFORM_FUNCTIONS;
         if (d_input_db->isDatabase("UniformFunctionData")) {
            getUniformFunctionDataFromInput(
               d_input_db->getDatabase("UniformFunctionData"));
         } else {
            TBOX_ERROR("Input error in " << d_object_name
                       << "\nNo 'UniformFunctionData' input input database found for test " 
                       << d_test_to_run << endl);
         }
      }

      if (func_dist_string == "RANDOM_FUNCTIONS") {
         d_function_distribution = RANDOM_FUNCTIONS;
         if (d_input_db->isDatabase("RandomFunctionData")) {
               getRandomFunctionDataFromInput(
                  d_input_db->getDatabase("RandomFunctionData"));
         } else {
            TBOX_ERROR("Input error in " << d_object_name
                       << "\nNo 'RandomFunctionData' input input database found for test "
                       << d_test_to_run << endl);
         }
      }

      if ( (func_dist_string == "FINE_PATCH_ALL") || 
           (func_dist_string == "FINE_PATCH_MOD2") ||
           (func_dist_string == "FINE_PATCH_EVEN_ODD") ) {
         if (func_dist_string == "FINE_PATCH_ALL") d_function_distribution = FINE_PATCH_ALL; 
         if (func_dist_string == "FINE_PATCH_MOD2") d_function_distribution = FINE_PATCH_MOD2; 
         if (func_dist_string == "FINE_PATCH_EVEN_ODD") d_function_distribution = FINE_PATCH_EVEN_ODD; 
         if (d_input_db->isDatabase("FinePatchFunctionData")) {
            getFinePatchFunctionDataFromInput(
               d_input_db->getDatabase("FinePatchFunctionData"));
         } else {
            TBOX_ERROR("Input error in " << d_object_name
                       << "\nNo 'FinePatchFunctionData' input database found for test "
                       << d_test_to_run << endl);
         }
      }
      
   } else {
      TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nNo 'function_distribution' item found!" << endl);
   }


}


/*************************************************************************
 *
 * Read gridding parameter data from input file.
 *
 ************************************************************************/

void LocallyActiveDataTester::getGriddingParametersFromInput(
   tbox::Pointer<tbox::Database> gridding_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!gridding_db.isNull());
#endif
 
   int ln;
 
   d_num_levels = 
      gridding_db->getIntegerWithDefault("num_levels", d_num_levels);
   if (d_num_levels < 1) {
      TBOX_ERROR("Input error for " << d_object_name << ":  "
                 << "Key data 'num_levels' found in input is < 1" << endl);
   }
   if ( (d_test_to_run == "COARSEN_TEST") && 
        d_num_levels > 2) {
      d_num_levels = 2;
      TBOX_WARNING("Setting 'num_levels' to 2 for COARSEN_TEST" << endl);
   }

   d_ratio_to_coarser.resizeArray(d_num_levels);
   d_ratio_to_coarser[0] = hier::IntVector<NDIM>(1);

   tbox::Pointer<tbox::Database> ratio_to_coarser_db;
   if ( d_num_levels > 1) {

      if (!gridding_db->keyExists("ratio_to_coarser")) {
         TBOX_ERROR("Input error for " << d_object_name << ":  "
                    << "Key data `ratio_to_coarser' not found in input.");
      } else {
         ratio_to_coarser_db = gridding_db->getDatabase("ratio_to_coarser");
         for (ln = 1; ln < d_num_levels; ln++) {

	    std::string level_name = "level_" + tbox::Utilities::intToString(ln);
 
            if (!ratio_to_coarser_db->keyExists(level_name)) {
                TBOX_ERROR("Input error for " << d_object_name << ":  "
                           <<"Key data `" << level_name
                           << "' not found in ratio_to_coarser input" << endl);
            }
 
            int* temp_ratio_to_coarser = d_ratio_to_coarser[ln];
            ratio_to_coarser_db->getIntegerArray(level_name,
                                                 temp_ratio_to_coarser, NDIM);
            if (d_ratio_to_coarser[ln].max() < 1) {
               TBOX_ERROR("Input error for " << d_object_name << ":  "
                          <<"Key data `" << level_name
                          << "' in ratio_to_coarser input has entry < 1" << endl);
            }
         }
      }

   }

   int* temp_npatches_on_coarsest = d_npatches_on_coarsest;
   gridding_db->getIntegerArray("npatches_on_coarsest",
                                 temp_npatches_on_coarsest, NDIM);
   if (d_npatches_on_coarsest.min() < 1) {
      TBOX_ERROR("Input error for " << d_object_name << ":  "
                 << "npatches_on_coarsest input has entries < 1" << endl);
   }

}

/*************************************************************************
 *
 * Functions to read locally-active data parameters from input file and
 * set up functions.
 *
 ************************************************************************/

void LocallyActiveDataTester::getFunctionSpecification(
   tbox::Pointer<tbox::Database> function_db)
{
   if (function_db->keyExists("function_spec")) {
      string fcn_spec_string = function_db->getString("function_spec");
      if ( (fcn_spec_string == "CONSTANT_TEST") ||
           (fcn_spec_string == "CONSTANT_PERF") ||
           (fcn_spec_string == "GAUSSIAN") ) { 
         if (fcn_spec_string == "CONSTANT_TEST") d_function_specification = CONSTANT_TEST;
         if (fcn_spec_string == "CONSTANT_PERF") d_function_specification = CONSTANT_PERF;
         if (fcn_spec_string == "GAUSSIAN") d_function_specification = GAUSSIAN;
      } else {
         TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nUnknown 'function_spec' string " << fcn_spec_string
                 << " found in database named " << function_db->getName() << endl); 
      }
   } else {
      TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nNo 'function_spec' item found in database named " 
                 << function_db->getName() << endl); 
   }
}

void LocallyActiveDataTester::getRandomFunctionDataFromInput(
   tbox::Pointer<tbox::Database> random_db)
{
   getFunctionSpecification(random_db);

   int default_num_test_functions = d_num_test_functions;
   d_num_test_functions = 
      d_input_db->getIntegerWithDefault("num_test_functions",
                                        default_num_test_functions);
   if (d_num_test_functions < 0) {
      TBOX_ERROR("Input error for " << d_object_name
         << "\n Negative 'num_test_functions' value found in 'RandomFunctionData' input!" 
         << endl);
   }

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
      d_hierarchy->getGridGeometry();
   const double* xlo = grid_geometry->getXLower();
   const double* xhi = grid_geometry->getXUpper();

   d_functions.resizeArray(d_num_test_functions);

   for (int fn = 0; fn < d_num_test_functions; fn++) {
      string fn_string = tbox::Utilities::intToString(fn) + "\0";

      d_functions[fn].d_name = "Function:" + fn_string;

      int id; 
      for (id = 0; id < NDIM; id++) {
         d_functions[fn].d_centroid[id] = drand48() * (xhi[id]-xlo[id]);
      }

      d_functions[fn].d_alpha = drand48();

      double min_len = xhi[0]-xlo[0];
      for (id = 1; id < NDIM; id++) {
         min_len = tbox::MathUtilities<double>::Min(min_len, xhi[id]-xlo[id]);
      }
      
      double fun1d = pow( (double)d_num_test_functions, (1.0/double(NDIM)) );
      d_functions[fn].d_radius = min_len/fun1d;

      d_functions[fn].d_func_data_index = -1;
      d_functions[fn].d_soln_data_index = -1;
   }
}

void LocallyActiveDataTester::getExplicitFunctionDataFromInput(
   tbox::Pointer<tbox::Database> explicit_data_db) 
{
   getFunctionSpecification(explicit_data_db);

   tbox::Array<string> db_names = explicit_data_db->getAllKeys();
   int d_num_test_functions = db_names.getSize() - 1;
   if (d_num_test_functions < 1) {
      TBOX_ERROR("Input error for " << d_object_name 
         << "\n No function data found in 'ExplicitFunctionData' input!" << endl);
   }
   d_functions.resizeArray(d_num_test_functions);

   int fn_id = -1;
   int name_count = 0;
   while ( (name_count < db_names.getSize()) && 
           (fn_id < d_num_test_functions) ) {
      if (db_names[name_count] != "function_spec") {
         fn_id++;
         d_functions[fn_id].d_name = db_names[name_count];            

         tbox::Pointer<tbox::Database> single_function_db = 
            explicit_data_db->getDatabase(d_functions[fn_id].d_name);

         double* tmp_centroid = d_functions[fn_id].d_centroid;
         single_function_db->getDoubleArray("centroid", tmp_centroid, NDIM);
         
         d_functions[fn_id].d_alpha = single_function_db->getDouble("alpha");
         d_functions[fn_id].d_radius = single_function_db->getDouble("radius");
         d_functions[fn_id].d_func_data_index = -1;
         d_functions[fn_id].d_soln_data_index = -1;
      }
      name_count++;
   }
}

void LocallyActiveDataTester::getUniformFunctionDataFromInput(
   tbox::Pointer<tbox::Database> uniform_db) 
{
   getFunctionSpecification(uniform_db);

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
   const double* xlo = grid_geometry->getXLower();
   const double* xhi = grid_geometry->getXUpper();

   double radius = 0.;
   double alpha = 0.;
   tbox::Array<int> nfunctions;
   tbox::Array<double> lower_bound(NDIM);
   tbox::Array<double> upper_bound(NDIM);

   radius = uniform_db->getDouble("radius");

   alpha = uniform_db->getDouble("alpha");
      
   nfunctions = uniform_db->getIntegerArray("nfunctions");
   if (nfunctions.getSize() < NDIM) {
      TBOX_ERROR("Input error for " << d_object_name << ":"
         << "\n   Size of 'nfunctions' array must = NDIM" << endl);
   }
 
   int id;
   for (id = 0; id < NDIM; id++) {
      lower_bound[id] = xlo[id];
      upper_bound[id] = xhi[id];
   }
   if (uniform_db->keyExists("lower_bound")) {
      lower_bound = uniform_db->getDoubleArray("lower_bound");
   }
   if (uniform_db->keyExists("upper_bound")) {
      upper_bound = uniform_db->getDoubleArray("upper_bound");
   }
   for (id = 0; id < NDIM; id++) {
      if (lower_bound[id] < xlo[id] || upper_bound[id] > xhi[id]) { 
         TBOX_ERROR("Input error for " << d_object_name << ":"
                    << "\n lower_bound or upper_bound out of range for domain" << endl);
      }
   }
      
   int num_functions = 1;
   for (id = 0; id < NDIM; id++) {
      num_functions *= nfunctions[id];
   }
   d_functions.resizeArray(num_functions);

   int n[3] = {0,0,0};
   int fn_id = -1;
      
#if (NDIM == 3)
   for (int nz = 0; nz < nfunctions[2]; nz++) {
      n[2] = nz;
#endif
      for (int ny = 0; ny < nfunctions[1]; ny++) {
         n[1] = ny;
         for (int nx = 0; nx < nfunctions[0]; nx++) {
            n[0] = nx;

            fn_id++; 
 
	    std::string fn_string = tbox::Utilities::intToString(fn_id) + "\0";
            d_functions[fn_id].d_name = "Function:" + fn_string;
            
            for (id = 0; id < NDIM; id++) {
               double incx = (upper_bound[id] - lower_bound[id]) /
                             (double)nfunctions[id];
               d_functions[fn_id].d_centroid[id] = 
                  lower_bound[id] + incx*n[id];
            }
               
            d_functions[fn_id].d_radius = radius;
            d_functions[fn_id].d_alpha = alpha;
               
            d_functions[fn_id].d_func_data_index = -1;
            d_functions[fn_id].d_soln_data_index = -1;
               
         }
      }
#if (NDIM == 3)
   }
#endif

}

void LocallyActiveDataTester::getFinePatchFunctionDataFromInput(
   tbox::Pointer<tbox::Database> fine_patch_data_db) 
{
   getFunctionSpecification(fine_patch_data_db);

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
   const double* xlo = grid_geometry->getXLower();
   const double* xhi = grid_geometry->getXUpper();

   double radius = fine_patch_data_db->getDouble("radius");
   double alpha = fine_patch_data_db->getDouble("alpha");

   int idim;

   int max_nfunctions[3] = {0, 0, 0};
   double function_spacing[3] = {0.0, 0.0, 0.0};
   int prod_max_nfunctions = 1;
   int actual_nfunctions = 0;

   for (idim = 0; idim < NDIM; idim++) {
      max_nfunctions[idim] = d_npatches_on_coarsest[idim];
      for (int ln = 1; ln < d_num_levels; ln++) {
         max_nfunctions[idim] *= d_ratio_to_coarser[ln](idim);
      }
      prod_max_nfunctions *= max_nfunctions[idim];
      function_spacing[idim] = (xhi[idim] - xlo[idim]) / (double)max_nfunctions[idim];
   }

   tbox::Array<LocallyActiveDataTester::Function> tmp_functions(prod_max_nfunctions);

   int n[3] = {0,0,0};
   int fn_id = -1;
   
   int loop_count = 0; 
#if (NDIM == 3)
   for (int nz = 0; nz < max_nfunctions[2]; nz++) {
      n[2] = nz;
#endif
      for (int ny = 0; ny < max_nfunctions[1]; ny++) {
         n[1] = ny;
         for (int nx = 0; nx < max_nfunctions[0]; nx++) {
            n[0] = nx;

            bool make_function = false;

            switch(d_function_distribution) {
               
               case FINE_PATCH_ALL: {
                  make_function = true;
                  break;
               }

               case FINE_PATCH_MOD2: {
                  int sum = 0;
                  for (idim = 0; idim < NDIM; idim++) {
                     sum += n[idim];
                  }
                  if (sum % 2 == 0) {
                     make_function = true;
                  }
                  break; 
               }

               case FINE_PATCH_EVEN_ODD: {
#if (NDIM == 2)
                  TBOX_ERROR("LocallyActiveDataTester::getFinePatchFunctionDataFromInput error :  "
                             << "FCN_SPEC_CASE = FINE_PATCH_EVEN_ODD only applies in 3d!" << endl);
#endif
                  if (n[2] % 2 == 0) {
                     if ( (n[0] % 2 == 0) && (n[1] % 2 == 0) ) {
                        make_function = true;
                     } 
                  } else {
                     if ( (n[0] % 2 != 0) && (n[1] % 2 != 0) ) {
                        make_function = true;
                     } 
                  }
                  break;
               }

               default: {

                  TBOX_ERROR("LocallyActiveDataTester::getFinePatchFunctionDataFromInput error :  "
                             << "unknown FCN_DIST_CASE = " << d_function_distribution << endl);

               }

            } //  switch(fcn_case)

            if (make_function) {

               fn_id++;

	       std::string fn_string = tbox::Utilities::intToString(loop_count) + "\0";
               tmp_functions[fn_id].d_name = "Function:" + fn_string;
            
               for (idim = 0; idim < NDIM; idim++) {
                  tmp_functions[fn_id].d_centroid[idim] = 
                     xlo[idim] + function_spacing[idim] * n[idim];
               }

               actual_nfunctions++;

            }

            loop_count++;
               
         }  // for nx = ...
      }  // for ny = ...
#if (NDIM == 3)
   } // for nz = ...
#endif

    d_functions.resizeArray(actual_nfunctions);
    for (fn_id = 0; fn_id < actual_nfunctions; fn_id++) {
       d_functions[fn_id].d_name = tmp_functions[fn_id].d_name; 
       d_functions[fn_id].d_radius = radius;
       d_functions[fn_id].d_alpha = alpha;
       for (idim = 0; idim < NDIM; idim++) {
          d_functions[fn_id].d_centroid[idim] = 
             tmp_functions[fn_id].d_centroid[idim];
       }
           
       d_functions[fn_id].d_func_data_index = -1;
       d_functions[fn_id].d_soln_data_index = -1;
    }
               
}

void LocallyActiveDataTester::printHierarchyData(
   ostream &os,
   bool dump_function_data) const
{

   hier::LocallyActiveVariableDatabase<NDIM>* var_db =
      hier::LocallyActiveVariableDatabase<NDIM>::getDatabase();

   for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<NDIM> > level_mgr =
          var_db->getLocallyActiveDataPatchLevelManager(level);

      os << "\n Patches on level: " << level->getLevelNumber() << endl;

      const hier::BoxArray<NDIM>& level_boxes = level->getBoxes();
      for (int ib = 0; ib < level_boxes.getNumberOfBoxes(); ib++) {
         os << "   Active functions on Patch<NDIM> " << ib << " : " 
            << level_boxes[ib] << endl << "      ";
         for (int fn = 0; fn < d_functions.getSize(); fn++) {
            int data_id = d_functions[fn].d_func_data_index;
            if (d_test_to_run == "COARSEN_TEST") {
               data_id = d_functions[fn].d_soln_data_index;
            }
            if ( level_mgr->getPatchDataActive( hier::PatchDataId(data_id), 
                                                hier::PatchNumber(ib) ) ) {
               os << fn << " , "; 
            }
         }
         os << endl;
      }
   }

   if (dump_function_data) {      
      for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ln++) {
         tbox::Pointer< hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

         os << "\n\nDumping function data on level: " << endl;  
      
         for (int fn = 0; fn < d_functions.getSize(); fn++) {

            int func_id = d_functions[fn].d_func_data_index;
            int soln_id = d_functions[fn].d_soln_data_index;

            os << "\n\nVariable: " << d_functions[fn].d_name << endl;

            hier::LocallyActiveDataPatchLevelManager<NDIM>::Iterator p =
               var_db->getLocallyActiveDataPatchLevelManager(level)->
               getIterator( hier::PatchDataId(func_id) );
            for ( ; p; p++) {
               tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());

               const hier::Index<NDIM> ifirst =patch->getBox().lower();
               const hier::Index<NDIM> ilast  =patch->getBox().upper();
               const tbox::Pointer< geom::CartesianPatchGeometry<NDIM> > patch_geom = 
                  patch->getPatchGeometry();

               if (d_test_to_run != "COARSEN_TEST") {
                  tbox::Pointer< pdat::CellData<NDIM,double> > func_data = 
                     patch->getPatchData(func_id);
                  hier::Box<NDIM> ghostbox = func_data->getGhostBox();

                  os << "\n  Function data on patch: " 
                     << p() << " : " << patch->getBox() << endl;
                  func_data->print(ghostbox, os);
               }

               if (d_test_to_run != "REFINE_TEST") { 
                  tbox::Pointer< pdat::CellData<NDIM,double> > soln_data = 
                     patch->getPatchData(soln_id);

                  os << "\n\n  Solution data on patch: " 
                     << p() << " : " << patch->getBox() << endl;
                  soln_data->print(patch->getBox(), os);
               }
            
            } // patch loop
         } // function loop
      } // level loop

   } // if dump data
   
}
