//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/LinAdv/main.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Main program for SAMRAI Linear Advection example problem.
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>

// Headers for basic SAMRAI objects

#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h" 
#include "tbox/InputDatabase.h" 
#include "tbox/InputManager.h" 
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h" 
#include "tbox/Utilities.h"
#include "VariableDatabase.h"

// Headers for major algorithm/data structure objects

#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "GriddingAlgorithm.h"
#include "HyperbolicLevelIntegrator.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "StandardTagAndInitialize.h"
#include "TimeRefinementIntegrator.h"
#include "TimeRefinementLevelStrategy.h"
#include "VisItDataWriter.h"

// Header for application-specific algorithm/data structure object

#include "LinAdv.h"

// Classes for run-time plotting and autotesting.

#if (TESTING == 1)
#include "AutoTester.h"
#endif

using namespace SAMRAI;


/************************************************************************
 *                                                                      *
 * This is the main program for an AMR solution of the linear advection *
 * equation: du/dt + div(a*u) = 0, where  "u" is a scalar-valued        *
 * function and "a" is a constant vector.  This application program is  *
 * constructed by composing several algorithm objects found in the      *
 * SAMRAI library with a few that are specific to this application.     *
 * A brief description of these object follows.                         *
 *                                                                      *
 * There are two main data containment objects.  These are:             *
 *                                                                      *
 *    hier::PatchHierarchy - A container for the AMR patch hierarchy and      *
 *       the data on the grid.                                          *
 *                                                                      *
 *    geom::CartesianGridGeometry<NDIM> - Defines and maintains the Cartesian       *
 *       coordinate system on the grid.  The hier::PatchHierarchy             *
 *       maintains a reference to this object.                          *
 *                                                                      *
 * A single overarching algorithm object drives the time integration    *
 * and adaptive gridding processes:                                     *
 *                                                                      *
 *    algs::TimeRefinementIntegrator<NDIM> - Coordinates time integration and       *
 *       adaptive gridding procedures for the various levels            *
 *       in the AMR patch hierarchy.  Local time refinement is          *
 *       employed during hierarchy integration; i.e., finer             *
 *       levels are advanced using smaller time increments than         *
 *       coarser level.  Thus, this object also invokes data            *
 *       synchronization procedures which couple the solution on        *
 *       different patch hierarchy levels.                              *
 *                                                                      *
 * The time refinement integrator is not specific to the numerical      *
 * methods used and the problem being solved.   It maintains references *
 * to two other finer grain algorithmic objects that are more specific  *
 * to the problem at hand and with which it is configured when they are *
 * passed into its constructor.   These finer grain algorithm objects   * 
 * are:                                                                 *
 *                                                                      *
 *    algs::HyperbolicLevelIntegrator<NDIM> - Defines data management procedures    *
 *       for level integration, data synchronization between levels,    *
 *       and tagging cells for refinement.  These operations are        *
 *       tailored to explicit time integration algorithms used for      *
 *       hyperbolic systems of conservation laws, such as the Euler     *
 *       equations.  This integrator manages data for numerical         *
 *       routines that treat individual patches in the AMR patch        *
 *       hierarchy.  In this particular application, it maintains a     *
 *       pointer to the LinAdv object that defines variables and        *
 *       provides numerical routines for the linear advection problem.  *
 *                                                                      *
 *       LinAdv - Defines variables and numerical routines for the      *
 *          discrete linear advection equation on each patch in the     * 
 *          AMR hierarchy.                                              *
 *                                                                      *
 *    mesh::GriddingAlgorithm<NDIM> - Drives the AMR patch hierarchy generation     *
 *       and regridding procedures.  This object maintains              *
 *       references to three other algorithmic objects with             *
 *       which it is configured when they are passed into its           *
 *       constructor.   They are:                                       *
 *                                                                      *
 *       mesh::BergerRigoutsos<NDIM> - Clusters cells tagged for refinement on a    *
 *          patch level into a collection of logically-rectangular      *
 *          box domains.                                                *
 *                                                                      *
 *       mesh::LoadBalancer<NDIM> - Processes the boxes generated by the            *
 *          mesh::BergerRigoutsos<NDIM> algorithm into a configuration from         *
 *          which patches are contructed.  The algorithm used in this   *
 *          class assumes a spatially-uniform workload distribution;    *
 *          thus, it attempts to produce a collection of boxes          *
 *          each of which contains the same number of cells.  The       *
 *          load balancer also assigns patches to processors.           *
 *                                                                      *
 *       mesh::StandardTagAndInitialize<NDIM> - Couples the gridding algorithm      *
 *          to the HyperbolicIntegrator. Selects cells for              *
 *          refinement based on either Gradient detection, Richardson   *
 *          extrapolation, or pre-defined Refine box region.  The       *
 *          object maintains a pointer to the algs::HyperbolicLevelIntegrator<NDIM>,*
 *          which is passed into its constructor, for this purpose.     *
 *                                                                      *
 ************************************************************************
 */

/*
 *******************************************************************
 *                                                                 *
 * For each run, the input filename and restart information        *
 * (if needed) must be given on the command line.                  *
 *                                                                 *
 *      For non-restarted case, command line is:                   *
 *                                                                 *
 *          executable <input file name>                           *
 *                                                                 *
 *      For restarted run, command line is:                        *
 *                                                                 *
 *          executable <input file name> <restart directory> \     *
 *                     <restart number>                            *
 *                                                                 *
 *******************************************************************
 */

int main( int argc, char *argv[])
{

   /*
    * Initialize tbox::MPI and SAMRAI, enable logging, and process command line.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   int num_failures = 0;

   {

      string input_filename;
      string restart_read_dirname;
      int restore_num = 0;

      bool is_from_restart = false;

      if ( (argc != 2) && (argc != 4) ) {
	 tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
	      << "<restart dir> <restore number> [options]\n"
	      << "  options:\n"
	      << "  none at this time"
	      << endl;
	 tbox::SAMRAI_MPI::abort();
	 return (-1);
      } else {
	 input_filename = argv[1];
	 if (argc == 4) {
	    restart_read_dirname = argv[2];
	    restore_num = atoi(argv[3]);

	    is_from_restart = true;
	 }
      }

      tbox::plog << "input_filename = " << input_filename << endl;
      tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
      tbox::plog << "restore_num = " << restore_num << endl;

      /*
       * Create input database and parse all data in input file.  
       */
   
      tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "GlobalInputs" section of the input database and set
       * values accordingly.
       */

      if (input_db->keyExists("GlobalInputs")) {
         tbox::Pointer<tbox::Database> global_db =
            input_db->getDatabase("GlobalInputs");
         if (global_db->keyExists("tag_clustering_method")) {
            string tag_clustering_method =
               global_db->getString("tag_clustering_method");
            mesh::BergerRigoutsos<NDIM>::setClusteringOption(tag_clustering_method);
         }
         if (global_db->keyExists("refine_schedule_generation_method")) {
            string refine_schedule_generation_method =
               global_db->getString("refine_schedule_generation_method");
            xfer::RefineSchedule<NDIM>::setScheduleGenerationMethod(
                                        refine_schedule_generation_method);
         }
         if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
            bool flag = global_db->
               getBool("call_abort_in_serial_instead_of_exit");
            tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
         }
      }

      /*
       * Retrieve "Main" section of the input database.  First, read dump
       * information, which is used for writing vizamrai plot files.  Second,
       * if proper restart information was given on command line, and the restart
       * interval is non-zero, create a restart database.
       */

      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

      string log_file_name = "linadv.log";
      if (main_db->keyExists("log_file_name")) {
	 log_file_name = main_db->getString("log_file_name");
      }
      bool log_all_nodes = false;
      if (main_db->keyExists("log_all_nodes")) {
	 log_all_nodes = main_db->getBool("log_all_nodes");
      }
      if (log_all_nodes) {
	 tbox::PIO::logAllNodes(log_file_name);
      } else {
	 tbox::PIO::logOnlyNodeZero(log_file_name);
      }

      int viz_dump_interval = 0;
      if (main_db->keyExists("viz_dump_interval")) {
	 viz_dump_interval = main_db->getInteger("viz_dump_interval");
      }

      tbox::Array<string> viz_writer(1);
      viz_writer[0] = "Vizamrai";
      string viz_dump_filename;
      string vizamrai_dump_dirname;
      string visit_dump_dirname;
      bool uses_vizamrai = false;
      bool uses_visit    = false;
      int visit_number_procs_per_file = 1;
      if ( viz_dump_interval > 0 ) {
	 if (main_db->keyExists("viz_writer")) {
	    viz_writer = main_db->getStringArray("viz_writer");
	 }
	 if (main_db->keyExists("viz_dump_filename")) {
	    viz_dump_filename = main_db->getString("viz_dump_filename");
	 }
	 string viz_dump_dirname;
	 if (main_db->keyExists("viz_dump_dirname")) {
	    viz_dump_dirname = main_db->getString("viz_dump_dirname");
	 }
	 for (int i = 0; i < viz_writer.getSize(); i++) {
	    if (viz_writer[i] == "Vizamrai") uses_vizamrai = true;
	    if (viz_writer[i] == "VisIt") uses_visit = true;
	 }
	 if (uses_vizamrai && !uses_visit) {
	    vizamrai_dump_dirname = viz_dump_dirname;
	 } else if (!uses_vizamrai && uses_visit) {
	    visit_dump_dirname = viz_dump_dirname;         
	 } else if (uses_vizamrai && uses_visit) {
	    vizamrai_dump_dirname = viz_dump_dirname + "_Vizamrai";
	    visit_dump_dirname = viz_dump_dirname + "_VisIt";
	 } else {
	    TBOX_ERROR("main(): "
		       << "\nUnrecognized 'viz_writer' entry..."
		       << "\nOptions are 'Vizamrai' and/or 'VisIt'" 
		       << endl);
	 }
	 if (uses_visit) {
	    if (viz_dump_dirname.empty()) {
	       TBOX_ERROR("main(): "
			  << "\nviz_dump_dirname is null ... "
			  << "\nThis must be specified for use with VisIt"
			  << endl);
	    }
	    if (main_db->keyExists("visit_number_procs_per_file")) {
	       visit_number_procs_per_file = 
		  main_db->getInteger("visit_number_procs_per_file");
	    }
	 }
	 if (uses_vizamrai) {
	    if (viz_dump_filename.empty()) {
	       TBOX_ERROR("main(): "
			  << "\nviz_dump_filename is null ... "
			  << "\nThis must be specified for use with Vizamrai"
			  << endl);
	    }
	 }
      }

      const bool viz_dump_data = (viz_dump_interval > 0);

      int restart_interval = 0;
      if (main_db->keyExists("restart_interval")) {
	 restart_interval = main_db->getInteger("restart_interval");
      } 

      string restart_write_dirname;
      if ( restart_interval > 0 ) {
	 if (main_db->keyExists("restart_write_dirname")) {
	    restart_write_dirname = main_db->getString("restart_write_dirname");
	 } else {
	    TBOX_ERROR("restart_interval > 0, but key `restart_write_dirname'"
		       << " not specifed in input file");
	 }
      }

      bool use_refined_timestepping = true;
      if ( main_db->keyExists("timestepping")) {
	 string timestepping_method = main_db->getString("timestepping");
	 if (timestepping_method == "SYNCHRONIZED") {
	    use_refined_timestepping = false;
	 }
      }

#if (TESTING == 1) && !(HAVE_HDF5)
      /* 
       * If we are autotesting on a system w/o HDF5, the read from
       * restart will result in an error.  We want this to happen
       * for users, so they know there is a problem with the restart,
       * but we don't want it to happen when autotesting.
       */
      is_from_restart  = false;
      restart_interval = 0;
#endif

      const bool write_restart = (restart_interval > 0)
	 && !(restart_write_dirname.empty());

      /*
       * Get the restart manager and root restart database.  If run is from 
       * restart, open the restart file.
       */
 
      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

      if (is_from_restart) {
	 restart_manager->
	    openRestartFile(restart_read_dirname, restore_num, 
			    tbox::SAMRAI_MPI::getNodes() );
      }

      /*
       * Create major algorithm and data objects which comprise application.
       * Each object will be initialized either from input data or restart
       * files, or a combination of both.  Refer to each class constructor
       * for details.  For more information on the composition of objects
       * for this application, see comments at top of file.
       */

      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
	 new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
				   input_db->getDatabase("CartesianGeometry"));

      tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy = 
	 new hier::PatchHierarchy<NDIM>("PatchHierarchy",grid_geometry);

      LinAdv* linear_advection_model = new LinAdv("LinAdv",
						  input_db->getDatabase("LinAdv"),
						  grid_geometry);

      tbox::Pointer<algs::HyperbolicLevelIntegrator<NDIM> > hyp_level_integrator =
	 new algs::HyperbolicLevelIntegrator<NDIM>(
	    "HyperbolicLevelIntegrator",
	    input_db->getDatabase("HyperbolicLevelIntegrator"),
	    linear_advection_model, true, use_refined_timestepping);

      tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector = 
	 new mesh::StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
				      hyp_level_integrator,
				      input_db->getDatabase("StandardTagAndInitialize"));

      tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();

      tbox::Pointer<mesh::LoadBalancer<NDIM> >  load_balancer = 
	 new mesh::LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));

      tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm = 
	 new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
			       input_db->getDatabase("GriddingAlgorithm"),
			       error_detector,
			       box_generator,
			       load_balancer);

      tbox::Pointer<algs::TimeRefinementIntegrator<NDIM> > time_integrator =
	 new algs::TimeRefinementIntegrator<NDIM>("TimeRefinementIntegrator",
				      input_db->getDatabase(
					 "TimeRefinementIntegrator"),
				      patch_hierarchy,
				      hyp_level_integrator,
				      gridding_algorithm);

      /*
       * Set up Visualization plot file writer(s).
       */
      tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > vizamrai_data_writer = 
	 new appu::CartesianVizamraiDataWriter<NDIM>("LinAdv Vizamrai Writer");

      // VisitDataWriter is only present if HDF is available
#ifdef HAVE_HDF5
      tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer = 
	 new appu::VisItDataWriter<NDIM>("LinAdv VisIt Writer",
			     visit_dump_dirname,
			     visit_number_procs_per_file);
      if (uses_visit) {
	 linear_advection_model->
	    registerVisItDataWriter(visit_data_writer);
      }
#endif

      if (uses_vizamrai) {
	 linear_advection_model->
	    registerVizamraiDataWriter(vizamrai_data_writer);
      }
   
      if (viz_dump_data) {
	 if (uses_vizamrai) {
	    vizamrai_data_writer->setDirectoryName(vizamrai_dump_dirname);
	    vizamrai_data_writer->setPlotDataToFloat();
	    vizamrai_data_writer->setFinestLevelToPlot(
	       gridding_algorithm->getMaxLevels()-1);
	    for (int ln = 1; ln < gridding_algorithm->getMaxLevels(); ln++) {
	       const hier::IntVector<NDIM>& lratio =
		  gridding_algorithm->getRatioToCoarserLevel(ln);
	       vizamrai_data_writer->setRatioToCoarserLevel(ln, lratio);
	    }
	 }
      }

      /*
       * Initialize hierarchy configuration and data on all patches.
       * Then, close restart file and write initial state for visualization.
       */

      double dt_now = time_integrator->initializeHierarchy();

      tbox::RestartManager::getManager()->closeRestartFile();

#if (TESTING == 1)
      /*
       * Create the autotesting object which will verify correctness
       * of the problem. If no automated testing is done, the object does 
       * not get used.
       */
      AutoTester autotester("AutoTester",input_db);
#endif

      /*
       * After creating all objects and initializing their state, we
       * print the input database and variable database contents
       * to the log file.
       */

      tbox::plog << "\nCheck input data and variables before simulation:" << endl;
      tbox::plog << "Input database..." << endl;
      input_db->printClassData(tbox::plog);
      tbox::plog << "\nVariable database..." << endl;
      hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

      tbox::plog << "\nCheck Linear Advection data... " << endl;
      linear_advection_model->printClassData(tbox::plog);

      if ( viz_dump_data ) {
	 if (uses_vizamrai) {
	    vizamrai_data_writer->writePlotData(
	       patch_hierarchy,
	       viz_dump_filename,
	       time_integrator->getIntegratorStep(),
	       time_integrator->getIntegratorTime()); 
	 }
#ifdef HAVE_HDF5
	 if (uses_visit) {
	    visit_data_writer->writePlotData(
	       patch_hierarchy,
	       time_integrator->getIntegratorStep(),
	       time_integrator->getIntegratorTime()); 
	 }
#endif
      }

      /*
       * Time step loop.  Note that the step count and integration
       * time are maintained by algs::TimeRefinementIntegrator<NDIM>.
       */

      double loop_time = time_integrator->getIntegratorTime();
      double loop_time_end = time_integrator->getEndTime();

      int iteration_num = time_integrator->getIntegratorStep();

#if (TESTING == 1)
      /*
       * If we are doing autotests, check result...
       */
      num_failures += autotester.evalTestData(iteration_num,
                                              patch_hierarchy,
                                              time_integrator,
                                              hyp_level_integrator,
                                              gridding_algorithm);
#endif

      while ( (loop_time < loop_time_end) && 
	      time_integrator->stepsRemaining() ) {

	 iteration_num = time_integrator->getIntegratorStep() + 1;

	 tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	 tbox::pout << "At begining of timestep # " <<  iteration_num - 1 << endl;
	 tbox::pout << "Simulation time is " << loop_time << endl;
	 tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
 
	 double dt_new = time_integrator->advanceHierarchy(dt_now); 

	 loop_time += dt_now;
	 dt_now = dt_new;

	 tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
	 tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
	 tbox::pout << "Simulation time is " << loop_time << endl;
	 tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
   
	 /* 
	  * At specified intervals, write restart and visualization files.
	  */
	 if ( write_restart ) {

	    if ( (iteration_num % restart_interval) == 0 ) {
	       tbox::RestartManager::getManager()->
		  writeRestartFile(restart_write_dirname,
				   iteration_num);
	    }
	 }

	 /*
	  * At specified intervals, write out data files for plotting.
	  */

	 if ( viz_dump_data ) {
	    if ( (iteration_num % viz_dump_interval) == 0 ) {
	       if (uses_vizamrai) {
		  vizamrai_data_writer->writePlotData(patch_hierarchy,
						      viz_dump_filename,
						      iteration_num,
						      loop_time);
	       }
#ifdef HAVE_HDF5
	       if (uses_visit) {
		  visit_data_writer->writePlotData(patch_hierarchy,
						   iteration_num,
						   loop_time);
	       }
#endif
	    }
	 }

#if (TESTING == 1)
	 /*
	  * If we are doing autotests, check result...
	  */
	 num_failures += autotester.evalTestData(iteration_num,
                                                 patch_hierarchy,
                                                 time_integrator,
                                                 hyp_level_integrator,
                                                 gridding_algorithm);
#endif

      }   

      /*
       * At conclusion of simulation, deallocate objects.
       */

#ifdef HAVE_HDF5
      visit_data_writer.setNull();
#endif

      vizamrai_data_writer.setNull();
      time_integrator.setNull();
      gridding_algorithm.setNull();
      load_balancer.setNull();
      box_generator.setNull(); 
      error_detector.setNull();
      hyp_level_integrator.setNull();

      if (linear_advection_model) delete linear_advection_model;

      patch_hierarchy.setNull();
      grid_geometry.setNull();

      input_db.setNull();
      main_db.setNull();

   }

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  LinAdv" << endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(num_failures);
}
