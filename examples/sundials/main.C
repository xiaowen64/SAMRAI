//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/sundials/main.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2141 $
// Modified:    $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Main program for testing Sundials/SAMRAI interface.
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

using namespace std;

using namespace std;

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include "tbox/SAMRAIManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

#include "tbox/Array.h"
#include "BergerRigoutsos.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "CellVariable.h"
#include "CellData.h"
#include "tbox/Database.h"
#include "StandardTagAndInitialize.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchLevel.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Utilities.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "LoadBalancer.h"
#include "VariableContext.h"
#include "VariableDatabase.h"

#include "SundialsAbstractVector.h"
#include "CVODESolver.h"
#include "Sundials_SAMRAIVector.h"
#include "CVODEModel.h"

// CVODE includes
#ifdef HAVE_SUNDIALS
#ifndef included_cvspgmr_h
#define included_cvspgmr_h
#include "cvode/cvode_spgmr.h"
#endif
#endif


using namespace SAMRAI;

/*
 * The cvode_test program is a general skeleton for using the 
 * CVODESolver interface.  The main stages of this program are:
 * 
 * (1)  Retrieving integration parameters from the input database.
 * (2)  Creating hierarchy, geometry, gridding, and CVODEModel
 *      objects.
 * (3)  Setting up the hierarchy configuration (grid configuration).
 * (4)  Setting the initial condition vector.
 * (5)  Creating a CVODESolver object.
 * (6)  Setting the integration parameters for CVODESolver.
 * (7)  Printing out to the log file the initial condition vector 
 *      and computing some norm for checking purposes.
 * (8)  Solving the ODE system.
 * (9)  Printing out to the log file the solution vector produced
 *      by CVODE and computing some norms.
 * (10) Printing out the CVODE statistics.
 * (11) Cleaning up the memory allocated for the program.
 */

int main( int argc, char *argv[] )
{

   /*
    * Initialize tbox::MPI and SAMRAI.  Enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   
   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("cvode_test.log");

#if !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)
      tbox::pout << "Library compiled WITHOUT CVODE -and- HYPRE...\n"
		 << "SAMRAI was not configured with one, or both, of "
		 << "these packages.  Cannot run this example."<< endl;
#else 

      /*
       * Process command line arguments. 
       */
      string input_filename;

      if (argc != 2) {
	 tbox::pout << "USAGE:  " << argv[0] << " <input filename> " << endl;
	 exit(-1);
      } else {
	 input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.  
       */
      tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename,input_db);

      /**************************************************************************
       * Read input data and setup objects.
       **************************************************************************/

      /*
       * Retreive "Main" section of input db.
       */
      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");
      int        max_order = main_db->getInteger("max_order");
      int        max_internal_steps = main_db->getInteger("max_internal_steps");
      double     init_time = main_db->getDouble("init_time");
      double     print_interval       = main_db->getDouble("print_interval");  
      int        num_print_intervals  = main_db->getInteger("num_print_intervals");


      double     relative_tolerance = main_db->getDouble("relative_tolerance");
      double     absolute_tolerance = main_db->getDouble("absolute_tolerance");
      bool       uses_newton = main_db->getBool("uses_newton");
      int        stepping_method = main_db->getInteger("stepping_method");
      bool       uses_preconditioning = 
	 main_db->getBoolWithDefault("uses_preconditioning", false);
      int        viz_dump_interval = 
	 main_db->getIntegerWithDefault("viz_dump_interval", 0);
      bool       solution_logging = 
	 main_db->getBoolWithDefault("solution_logging", false);
   
      string viz_dump_filename;
      string viz_dump_dirname;
      if ( viz_dump_interval > 0 ) {
	 if (main_db->keyExists("viz_dump_filename")) {
	    viz_dump_filename = main_db->getString("viz_dump_filename");
	 }
	 if (main_db->keyExists("viz_dump_dirname")) {
	    viz_dump_dirname = main_db->getString("viz_dump_dirname");
	 }
      }
      const bool viz_dump_data = (viz_dump_interval > 0)
	 && !(viz_dump_filename.empty());

      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > viz_data_writer =
	 new appu::CartesianVizamraiDataWriter<NDIM>("TwoTemperature Viz Writer");
   
      /*
       * Create geometry and hierarchy objects.
       */
      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > geometry = new geom::CartesianGridGeometry<NDIM>(
	 "Geometry",
	 input_db->getDatabase("Geometry"));

      tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = new hier::PatchHierarchy<NDIM>(
	 "Hierarchy",
	 geometry);

      /*
       * Create gridding algorithm objects that will handle construction of
       * of the patch levels in the hierarchy.
       */
      tbox::Pointer<CVODEModel> cvode_model = 
	 new CVODEModel("CVODEModel",
			input_db->getDatabase("CVODEModel"),
			geometry);
   
      tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_est =
	 new mesh::StandardTagAndInitialize<NDIM>(
	    "StandardTagAndInitialize",
	    cvode_model,
	    input_db->getDatabase("StandardTagAndInitialize"));

      tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();
  
      tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = new mesh::LoadBalancer<NDIM>(
	 "LoadBalancer",
	 input_db->getDatabase("LoadBalancer"));

      tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
	 new mesh::GriddingAlgorithm<NDIM>(
	    "GriddingAlgorithm",
	    input_db->getDatabase("GriddingAlgorithm"),
	    error_est,
	    box_generator,
	    load_balancer);

      /*
       * Set up Vizamrai plot file writer.
       */
      if (viz_dump_data) {
	 viz_data_writer->setDirectoryName(viz_dump_dirname);
	 viz_data_writer->setPlotDataToFloat();
	 viz_data_writer->
	    setFinestLevelToPlot(gridding_algorithm->getMaxLevels()-1);
	 for (int ln = 1; ln < gridding_algorithm->getMaxLevels(); ln++)
	 {
	    const hier::IntVector<NDIM>& lratio =
	       gridding_algorithm->getRatioToCoarserLevel(ln);
	    viz_data_writer->setRatioToCoarserLevel(ln, lratio);
	 }
      }
   
      /*
       * Setup hierarchy.  
       */
      gridding_algorithm->makeCoarsestLevel(hierarchy,init_time);
   
      tbox::Array<int> tag_buffer_array(gridding_algorithm->getMaxLevels());
      for (int il = 0; il < gridding_algorithm->getMaxLevels(); il++) {
	 tag_buffer_array[il] = 1;
      }

      bool done = false;
      bool initial_time = true;
      for (int ln = 0; gridding_algorithm->levelCanBeRefined(ln) && !done; ln++) {
	 gridding_algorithm->makeFinerLevel(hierarchy,
					    init_time,
					    initial_time,
					    tag_buffer_array[ln]);
	 done = !(hierarchy->finerLevelExists(ln));
      }

      /*
       * Setup timer manager for profiling code.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));  
      tbox::Pointer<tbox::Timer> t_cvode_solve = tbox::TimerManager::getManager()->
	 getTimer("apps::main::cvode_solver");
      tbox::Pointer<tbox::Timer> t_viz_dump = tbox::TimerManager::getManager()->
	 getTimer("apps::main::Vizamrai dump");
      tbox::Pointer<tbox::Timer> t_log_dump = tbox::TimerManager::getManager()->
	 getTimer("apps::main::Solution log dump");
      /*
       * Setup solution vector. 
       */
      cvode_model->setupSolutionVector(hierarchy);
      SundialsAbstractVector* solution_vector = 
	 cvode_model->getSolutionVector();

      /*
       * Set initial conditions vector.
       */
      cvode_model->setInitialConditions(solution_vector);
   
      /**************************************************************************
       * Setup CVODESolver object.
       **************************************************************************/
      solv::CVODESolver* cvode_solver = 
	 new solv::CVODESolver("cvode_solver",
			       cvode_model,
			       uses_preconditioning);

      int neq = 0;
      tbox::Pointer<hier::PatchLevel<NDIM> > level_zero 
	 = hierarchy->getPatchLevel(0);
      hier::BoxArray<NDIM> level_0_boxes = level_zero->getBoxes();
      for (int i = 0; i < level_0_boxes.getNumberOfBoxes(); i++) {
	 neq += level_0_boxes[i].size();
      }
      cvode_solver->setIterationType(uses_newton ? CV_NEWTON : CV_FUNCTIONAL);
      //cvode_solver->setToleranceType(SV); // this is in craig's code, but
      // causes mine to bomb.  Why??
      cvode_solver->setRelativeTolerance(relative_tolerance);
      cvode_solver->setAbsoluteTolerance(absolute_tolerance);
      cvode_solver->setMaximumNumberOfInternalSteps(max_internal_steps);
      cvode_solver->setSteppingMethod(stepping_method);
      cvode_solver->setMaximumLinearMultistepMethodOrder(max_order);
      if (uses_preconditioning) {
	 cvode_solver->setPreconditioningType(PREC_LEFT);
      }


      cvode_solver->setInitialValueOfIndependentVariable(init_time);
      cvode_solver->setInitialConditionVector(solution_vector);
      cvode_solver->initialize(solution_vector);

      /*
       * Print initial vector (if solution logging is enabled)
       */
      tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > y_init = 
	 solv::Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(solution_vector);

      if (solution_logging) {

	 tbox::Pointer<hier::PatchHierarchy<NDIM> > init_hierarchy =
	    y_init->getPatchHierarchy();
      
	 tbox::pout << "Initial solution vector y() at initial time: " << endl;
	 int ln;
	 tbox::pout << "y(" << init_time << "): "<< endl;
	 for (ln = 0; ln < init_hierarchy->getNumberOfLevels(); ln++) {
	    tbox::Pointer<hier::PatchLevel<NDIM> > level;
         
	    level = init_hierarchy->getPatchLevel(ln);
	    tbox::pout << "level = " << ln << endl;
         
	    for (hier::PatchLevel<NDIM>::Iterator  p(level); p; p++) {
	       tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
            
	       tbox::Pointer< CellData<NDIM,double> > y_data =
		  y_init->getComponentPatchData(0,*patch);
//         y_data->print(y_data->getBox(),tbox::pout); 
	    }
	 }
      }

      /*
       * Compute maxNorm and L1Norm of initial vector
       */
      //pout.precision(12);
      if (solution_logging) {
	 tbox::pout << "\n\nBefore solve..." << endl;
	 tbox::pout << "Max Norm of y()= " << y_init->maxNorm() << endl;
	 tbox::pout << "L1 Norm of y()= " << y_init->L1Norm() << endl;
	 tbox::pout << "L2 Norm of y()= " << y_init->L2Norm() << endl;
      }

      /*
       * Write initial data to vizamrai file
       */
      if ( viz_dump_data ) {
   
	 int y_id = y_init->getComponentDescriptorIndex(0);

	 viz_data_writer->registerPlotScalar("y", y_id, 0);
 
	 tbox::Pointer<hier::PatchHierarchy<NDIM> > init_hierarchy =
	    y_init->getPatchHierarchy();

	 viz_data_writer->writePlotData(init_hierarchy,
					viz_dump_filename,
					0,
					init_time);
      }


      /**************************************************************************
       * Start time-stepping.
       **************************************************************************/

      tbox::Array<double> time(num_print_intervals);
      tbox::Array<double> maxnorm(num_print_intervals);
      tbox::Array<double> l1norm(num_print_intervals);
      tbox::Array<double> l2norm(num_print_intervals);

      double final_time = init_time;
      int interval;
      for (interval = 1; interval <= num_print_intervals; interval++) {

	 /*
	  * Set time interval
	  */
	 final_time += print_interval;
	 cvode_solver->setFinalValueOfIndependentVariable(final_time, false);
      
	 /*
	  * Perform CVODE solve to the requested interval time.
	  */
	 t_cvode_solve->start();
	 tbox::plog << "return code = " << cvode_solver->solve() << endl;
	 t_cvode_solve->stop();
	 double actual_time = 
	    cvode_solver->getActualFinalValueOfIndependentVariable();
      
	 /*
	  * Print statistics
	  * Format:  time  max norm   l1 norm   l2 norm
	  */
	 tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > y_result = 
	    solv::Sundials_SAMRAIVector<NDIM>::getSAMRAIVector(solution_vector);
	 tbox::Pointer<hier::PatchHierarchy<NDIM> > result_hierarchy =
	    y_result->getPatchHierarchy();

 
	 time[interval-1] = actual_time;
	 maxnorm[interval-1] = y_result->maxNorm();
	 l1norm[interval-1] = y_result->L1Norm();
	 l2norm[interval-1] = y_result->L2Norm();
      
	 if ( solution_logging ) cvode_solver->printStatistics(tbox::pout);
      
	 /*
	  * Dump vizamrai data.
	  */
	 if ( viz_dump_data ){
	    if ( (interval % viz_dump_interval) == 0 ){
	       t_viz_dump->start();
	       viz_data_writer->writePlotData(result_hierarchy,
					      viz_dump_filename,
					      interval,
					      final_time);
	       t_viz_dump->stop();
	    }
	 }
      
	 /*
	  * Write solution (if desired).
	  */
	 if ( solution_logging ) {
	    tbox::plog << "y(" << final_time << "): "<< endl << endl;
	    t_log_dump->start();
	    for (int ln = 0; ln < result_hierarchy->getNumberOfLevels(); ln++)
	    {
	       tbox::Pointer<hier::PatchLevel<NDIM> > level;
            
	       level = result_hierarchy->getPatchLevel(ln);
	       tbox::plog << "level = " << ln << endl;
            
	       for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
	       {
		  tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
               
		  tbox::Pointer< CellData<NDIM,double> > y_data =
		     y_result->getComponentPatchData(0,*patch);
		  y_data->print(y_data->getBox());
	       }
	    }
	    t_log_dump->stop();
	 }
      } // end of timestep loop
   
      /*************************************************************************
       * Write summary information
       ************************************************************************/
      /*
       * Write CVODEModel stats
       */
      tbox::Array<int> counters;
      cvode_model->getCounters(counters);

#if (TESTING==1)
      int correct_rhs_evals      = main_db->getInteger("correct_rhs_evals");
      int correct_precond_setups = main_db->getInteger("correct_precond_setups");
      int correct_precond_solves = main_db->getInteger("correct_precond_solves");
   
      if (counters[0] == correct_rhs_evals) {
	 tbox::pout << "Test 0: Number RHS evals CORRECT" << endl;
      } else { 
	 tbox::perr << "Test 0 FAILED: Number RHS evals INCORRECT" << endl;
	 tbox::perr << "Correct Number RHS evals:  " << correct_rhs_evals << endl;
	 tbox::perr << "Number RHS evals computed: " << counters[0] << endl;
      }

      if (counters[1] == correct_precond_setups) {
	 tbox::pout << "Test 1: Number precond setups CORRECT" << endl;
      } else { 
	 tbox::perr << "Test 1 FAILED: Number precond setups INCORRECT" << endl;
	 tbox::perr << "Correct number precond setups:  " << correct_precond_setups << endl;
	 tbox::perr << "Number precond setups computed: " << counters[1] << endl;
      }

      if (counters[2] == correct_precond_solves) {
	 tbox::pout << "Test 2: Number precond solves CORRECT" << endl;
      } else { 
	 tbox::perr << "Test 2 FAILED: Number precond solves INCORRECT" << endl;
	 tbox::perr << "Correct number precond solves:  " << correct_precond_solves << endl;
	 tbox::perr << "Number precond solves computed: " << counters[2] << endl;
      }
#endif


      if (solution_logging) {
	 tbox::pout << "\n\nEnd Timesteps - final time = " << final_time
		    << "\n\tTotal number of RHS evaluations = " << counters[0]
		    << "\n\tTotal number of precond setups = " << counters[1]
		    << "\n\tTotal number of precond solves = " << counters[2]
		    << endl;
   
	 /*
	  * Write out timestep sequence information
	  */   
	 tbox::pout << "\n\nTimestep Summary of solution vector y()\n"
		    << "  time                   \t"
		    << "  Max Norm  \t"
		    << "  L1 Norm  \t"
		    << "  L2 Norm  " << endl;

	 for (interval = 0; interval < num_print_intervals; interval++) {
	    tbox::pout.precision(18);
	    tbox::pout << "  "  << time[interval] << "  \t";
	    tbox::pout.precision(6);
	    tbox::pout << "  " << maxnorm[interval] << "  \t"
		       << "  " << l1norm[interval] << "  \t"
		       << "  " << l2norm[interval] << endl;
	 }
      }
   
      /*
       * Write out timings
       */   
#if (TESTING != 1)
      tbox::TimerManager::getManager()->print(tbox::pout);
#endif

      /*
       * Memory cleanup.
       */
      if (cvode_solver) delete cvode_solver;
   
      cvode_model.setNull();
      gridding_algorithm.setNull();
      error_est.setNull();
      load_balancer.setNull();
      box_generator.setNull();
      hierarchy.setNull();
      geometry.setNull();
      
#endif // HAVE_SUNDIALS

      tbox::pout << "\nPASSED:  cvode" << endl;

   }

   /*
    * Shutdown SAMRAI and tbox::MPI.
    */
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(0); 
}

