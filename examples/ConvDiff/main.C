//
// File:        main.C
// Package:     SAMRAI application
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Main program for SAMRAI convection-diffusion ex. problem.
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
#include "tbox/Array.h"
#include "CellData.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "MainRestartData.h"
#include "tbox/MPI.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"

// Headers for major algorithm/data structure objects
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CartesianVizamraiDataWriter.h"
#include "GriddingAlgorithm.h"
#include "StandardTagAndInitialize.h"
#include "MethodOfLinesIntegrator.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "LoadBalancer.h"

// Header for application-specific algorithm/data structure object
#include "ConvDiff.h"

// Classes for autotesting.
#if (TESTING == 1)
#include "AutoTester.h"
#endif

using namespace SAMRAI;
using namespace algs;

/************************************************************************
 *                                                                      *
 * This is the main program for an AMR application of the method        *
 * of lines algorithm.  Specifically, it uses a Runge-Kutta method      *
 * to integrate the convection-diffusion in time with spatial           *
 * refinement.  The timestep loop in this file controls gridding        *
 * operations.                                                          *
 *                                                                      *
 * The gridding objects used in this application are as follows:        *
 *                                                                      *
 *    hier::PatchHierarchy<NDIM> - A container for the AMR patch hierarchy and      *
 *       the data on the grid.                                          *
 *                                                                      *
 *    geom::CartesianGridGeometry<NDIM> - Defines and maintains the Cartesian       *
 *       coordinate system on the grid.  The hier::PatchHierarchy<NDIM>             *
 *       maintains a reference to this object.                          *
 *                                                                      *
 * The Method-of-Lines algorithm uses these two components to integrate *
 * data one timestep.  Time-stepping is done in the main routine.       *
 * This case uses the same timestep across all patches, regardless      *
 * of the refinement.  See either the Euler or Linear Advection         *
 * cases for examples of time-refinement cases.                         *
 *                                                                      *
 *    algs::MethodOfLinesIntegrator<NDIM> - advances solution on patches.           *
 *             Controls integration of quantities over entire patch     *
 *             hierarchy. Takes the user-defined class that defines     *
 *             patch operations (ConvDiff in this case) as an argument  *
 *                                                                      *
 * The user-supplied object defines characteristics and operations to   *
 * be carried out on each patch.  We pass this object to the integrator *
 * which will orchestrate the integration on all patches.               *
 *                                                                      *
 *    ConvDiff - Defines variables and numerical routines on a single   *
 *               patch.                                                 *
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
 *          which patches are contructed.  The algorithm we use in this *
 *          class assumes a spatially-uniform workload distribution;    *
 *          thus, it attempt to produce a collection of boxes           *
 *          each of which contains the same number of cells.  The       *
 *          load balancer also assigns patches to processors.           *
 *                                                                      *
 *       mesh::StandardTagAndInitialize<NDIM> - Couples the gridding algorithm      *
 *          to the HyperbolicIntegrator. Selects cells for              *
 *          refinement based on either Gradient detection, Richardson   *
 *          extrapolation, or pre-defined Refine box region.  The       *
 *          object maintains a pointer to the algs::MethodOfLinesIntegrator<NDIM>,  *
 *          which is passed into its constructor, for this purpose.     *
 *                                                                      *
 ************************************************************************
 */


int main( int argc, char *argv[] )
{
   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();


   /* This extra code block is used to scope some temporaries that are
    * created, it forces the destruction before the manager is
    * shutdown.
    */
   {

   /*
    * Process command line arguments and dump to log file.
    * For non-restarted case, command line is:
    *
    *    executable <input file name>
    *
   * For restarted run, command line is:
    *
    *    executable <input file name> <restart directory> \
    *               <restart number>
    */

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
        tbox::MPI::abort();
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
   tbox::InputManager::getManager()->parseInputFile(input_filename,input_db);

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
         tbox::MPI::setCallAbortInSerialInsteadOfExit(flag);
      }
   }

   /*
    * Retrieve "Main" section of the input database.  First, read dump
    * information, which is used for writing plot files.  Second, 
    * if proper restart information was given on command line, and the restart
    * interval is non-zero, create a restart database.
    */

   tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

   string log_file_name = "convdiff.log";
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
   tbox::Array<string> *viz_writer = new tbox::Array<string>(1);
   string viz_dump_filename;
   string vizamrai_dump_dirname;
   string visit_dump_dirname;
   bool uses_vizamrai = false;
   bool uses_visit    = false;
   int visit_number_procs_per_file = 1;
   if ( viz_dump_interval > 0 ) {
      if (main_db->keyExists("viz_writer")) {
         *viz_writer = main_db->getStringArray("viz_writer");
      } else {
	 (*viz_writer)[0] = "Vizamrai";
      }

      if (main_db->keyExists("viz_dump_filename")) {
         viz_dump_filename = main_db->getString("viz_dump_filename");
      }
      string viz_dump_dirname;
      if (main_db->keyExists("viz_dump_dirname")) {
         viz_dump_dirname = main_db->getString("viz_dump_dirname");
      }
      for (int i = 0; i < viz_writer -> getSize(); i++) {
         if ((*viz_writer)[i] == "Vizamrai") uses_vizamrai = true;
         if ((*viz_writer)[i] == "VisIt") uses_visit = true;
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
         TBOX_ERROR("`restart_interval' > 0, but key `restart_write_dirname'"
                    << "not found in input file.");
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
    * Get restart manager and root restart database.  If run is from
    * restart, open the restart file.
    */

   tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

   if (is_from_restart) {
     restart_manager->
         openRestartFile(restart_read_dirname, restore_num, 
                         tbox::MPI::getNodes() );
   }
   

   /* 
    * Initialize the MainRestartData object which stores the state of the
    * main program for restart.
    */
   MainRestartData* main_restart_data = new MainRestartData(
      "MainRestartData",
      input_db->getDatabase("MainRestartData"));

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
      new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   ConvDiff* convdiff_model = new ConvDiff("ConvDiff",
                                           input_db->getDatabase("ConvDiff"),
                                           grid_geometry);

   tbox::Pointer<algs::MethodOfLinesIntegrator<NDIM> > mol_integrator =
      new algs::MethodOfLinesIntegrator<NDIM>(
             "MethodOfLinesIntegrator",
             input_db->getDatabase("MethodOfLinesIntegrator"),
             convdiff_model);

   tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
      new mesh::StandardTagAndInitialize<NDIM>(
             "StandardTagAndInitialize",
             mol_integrator,
             input_db->getDatabase("StandardTagAndInitialize"));

   tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();

   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = 
      new mesh::LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));

   tbox::Pointer< mesh::GriddingAlgorithm<NDIM> > gridding_algorithm = 
      new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            input_db->getDatabase("GriddingAlgorithm"),
                            error_detector,
                            box_generator,
                            load_balancer);
   
   /*
    * Set up Visualization plot file writer(s).
    */
   tbox::Pointer<appu::CartesianVizamraiDataWriter<NDIM> > vizamrai_data_writer = 
      new appu::CartesianVizamraiDataWriter<NDIM>("ConvDiff Vizamrai Writer");

#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer = 
      new appu::VisItDataWriter<NDIM>("ConvDiff VisIt Writer",
                          visit_dump_dirname,
                          visit_number_procs_per_file);
   if (uses_visit) {
      convdiff_model->registerVisItDataWriter(visit_data_writer);
   }
#endif
   if (uses_vizamrai) {
      convdiff_model->registerVizamraiDataWriter(vizamrai_data_writer);
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
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents to
    * the log file.
    */

   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << endl;
   hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);
   mol_integrator->initializeIntegrator(gridding_algorithm);

   /****************************************************************
    *                                                              *
    *  INITIALIZE DATA ON PATCHES                                  *
    *                                                              *
    ****************************************************************
    *                                                              *
    *  Build patch hierarchy and initialize the data on the patches*
    *  in the hierarchy. Note: this step is performed by the       *
    *  algs::TimeRefinementIntegrator<NDIM> in Euler/LinAdv example cases.     *
    *  1) Create a "tag_buffer" for each level in the Hierarchy.   *
    *  2) Create the coarse (i.e. level 0) grid.                   *
    *  3) Cycle through levels 1-max_levels, initializing data     *
    *     on each.  The makeFinerLevel method calls the error      *
    *     estimator (remember, it was registered with the          *
    *     gridding algorithm object) and tags cells for refinement *
    *     as it generates patches on the finer levels.             *
    *  4) Dump this initial data to viz file.                      *
    *                                                              *
    ****************************************************************/


   tbox::Array<int> *tag_buffer_array = new tbox::Array<int>(gridding_algorithm->getMaxLevels());
   for (int il = 0; il < gridding_algorithm->getMaxLevels(); il++) {
      (*tag_buffer_array)[il] = main_restart_data->getTagBuffer();
      tbox::pout << "il = " << il << " tag_buffer = " << (*tag_buffer_array)[il] 
           << endl;
   }

   double loop_time = main_restart_data->getLoopTime();

   if (tbox::RestartManager::getManager()->isFromRestart()) {

      patch_hierarchy->getFromRestart(gridding_algorithm->getMaxLevels());
  
      gridding_algorithm->getTagAndInitializeStrategy()->
         resetHierarchyConfiguration(patch_hierarchy, 
                                     0, 
                                     patch_hierarchy->getFinestLevelNumber()); 

   } else {

      gridding_algorithm->makeCoarsestLevel(patch_hierarchy,loop_time);

      bool done = false;
      bool initial_time = true;
      for (int ln = 0; 
           gridding_algorithm->levelCanBeRefined(ln) && !done; 
           ln++) {
          gridding_algorithm->makeFinerLevel(patch_hierarchy,
                               loop_time,
                               initial_time,
                               (*tag_buffer_array)[ln]);
          done = !(patch_hierarchy->finerLevelExists(ln));
      }
   }

   tbox::RestartManager::getManager()->closeRestartFile();

#if (TESTING == 1)
   /*
    * Create the autotesting component which will verify correctness
    * of the problem. If no automated testing is done, the object does
    * not get used.
    */
   AutoTester *autotester = new AutoTester("AutoTester",input_db);
#endif

   /*******************************************************************
    *                                                                 *
    *  MAIN TIME ADVANCE LOOP                                         *
    *                                                                 *
    *******************************************************************
    *                                                                 *
    *  1) Set start and end time.                                     *
    *  2) Start integration timesteps.                                *
    *     While (loop_time < end_time) {                              *
    *        2a) Write restart and vizamrai data.                     *
    *        2b) Advance all levels in the hierarchy by time          *
    *            dt by calling algs::MethodOfLinesIntegrator<NDIM>'s              *
    *            advanceHierarchy method.                             *
    *        2c) Check if it is time to do a regrid step.  If so,     *
    *            have the mesh::GriddingAlgorithm<NDIM> call its regridAllFiner   *
    *            Levels method, passing in the coarsest (0) level.    *
    *            This method will invoke the Gradient detector,       *
    *            Berger Rigoutsos algorithm, and load balancer        *
    *            while generating finer grid levels.                  *
    *     }                                                           *
    *                                                                 *
    ******************************************************************/
   int iteration_num = main_restart_data->getIterationNumber();

#if (TESTING == 1)
      /*
       * If we are doing autotests, check result...
       */
      autotester -> evalTestData( iteration_num,
                               patch_hierarchy,
                               loop_time,
                               mol_integrator,
                               gridding_algorithm);
#endif

   if ( viz_dump_data ) {
      if (uses_vizamrai) {
         vizamrai_data_writer->writePlotData(
            patch_hierarchy,
            viz_dump_filename,
            iteration_num,
            loop_time); 
      }

#ifdef HAVE_HDF5
      if (uses_visit) {
         visit_data_writer->writePlotData(
            patch_hierarchy,
            iteration_num,
            loop_time); 
      }
#endif

   }
   
   while ( (loop_time < main_restart_data->getEndTime()) && 
           (iteration_num < main_restart_data->getMaxTimesteps()) ) {
  
      iteration_num = main_restart_data->getIterationNumber();
      iteration_num++;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "\n++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      double dt = mol_integrator->getTimestep(patch_hierarchy, loop_time);

      mol_integrator->advanceHierarchy(patch_hierarchy, loop_time, dt);

      loop_time += dt;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      /*
       * Write restart file at specified intervals.  Set current state 
       * of "main" in the main_restart_data object, before writing restart.
       */
      main_restart_data->setLoopTime(loop_time);
      main_restart_data->setIterationNumber(iteration_num);
      if ( write_restart ) {
         if ( (iteration_num % restart_interval) == 0 ) {
            tbox::RestartManager::getManager()->
               writeRestartFile(restart_write_dirname, 
				iteration_num);
         }
      }

      /*
       * At specified intervals, write out data files for plotting.
       * The viz_data_writer dumps data in a format that can
       * be processed by the Vizamrai tool. 
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

      /* 
       *  At specified intervals, regrid.
       */
      if ((iteration_num % main_restart_data->getRegridStep()) == 0 && 
           gridding_algorithm->getMaxLevels() > 1 ) {
         tbox::pout << "\n\n############################################" << endl;
         tbox::pout << "                 REGRIDDING" << endl;
         tbox::pout << "Finest level before regrid: " 
              << patch_hierarchy->getFinestLevelNumber() << endl;

         gridding_algorithm->regridAllFinerLevels(patch_hierarchy,0,
                                          loop_time,
                                          *tag_buffer_array);

         tbox::pout << "Finest level after regrid: " 
              << patch_hierarchy->getFinestLevelNumber() << endl;
         tbox::pout << "############################################\n\n" << endl;
      }

#if (TESTING == 1)
      /*
       * If we are doing autotests, check result...
       */
      autotester -> evalTestData( iteration_num,
                               patch_hierarchy,
                               loop_time,
                               mol_integrator,
                               gridding_algorithm);
#endif

   }




   /*
    * At conclusion of simulation, deallocate objects.
    */

   box_generator.setNull();


   load_balancer.setNull();
   gridding_algorithm.setNull();
   vizamrai_data_writer.setNull();

#ifdef HAVE_HDF5
   visit_data_writer.setNull();
#endif


   delete viz_writer;
   delete tag_buffer_array;


   error_detector.setNull();
   mol_integrator.setNull();

   if (convdiff_model) delete convdiff_model;

   patch_hierarchy.setNull();



   grid_geometry.setNull();

   if (main_restart_data) delete main_restart_data;

   main_db.setNull();
   input_db.setNull();

#if (TESTING == 1)
   delete autotester;
#endif

   } 

   tbox::pout << "\nPASSED:  ConvDiff" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();
 
   return(0); 
}



