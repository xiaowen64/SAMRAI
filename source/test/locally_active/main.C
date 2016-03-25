//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/locally_active/main.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Tests multi-function operations
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>


/*
 * Header files for classes particular to this application
 */
#include "LocallyActiveDataTester.h"

/*
 * Header files for SAMRAI library classes
 */
#include "LocallyActiveVariableDatabase.h"
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/IOStream.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/SAMRAIManager.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

int main( int argc, char *argv[] ) {

   int fail_count;

   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {


      /*
       * Process command line arguments and dump to log file.
       * To run:
       *
       *    executable <input file name>
       */

      string input_filename;

      if (argc != 2) {
	 tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
		    << "  options:\n"
		    << "  none at this time"
		    << endl;
	 tbox::SAMRAI_MPI::abort();
	 return (-1);
      }
 
      input_filename = argv[1];
      tbox::plog << "input_filename = " << input_filename << endl;


      /**********************************************************************
       * Read input
       *********************************************************************/

      /*
       * Create input database and parse all data in input file.
       */
      tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename,input_db);

      /*
       * Retrieve "Main" section of the input database.  First, read dump
       * information, which is used for writing plot files.  Second,
       * if proper restart information was given on command line, and the restart
       * interval is non-zero, create a restart database.
       */
      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

      string log_file_name = 
	 main_db->getStringWithDefault("log_file_name", "locally-active-test.log");
      bool log_all_nodes = 
	 main_db->getBoolWithDefault("log_all_nodes", false);
      if (log_all_nodes) {
	 tbox::PIO::logAllNodes(log_file_name);
      } else {
	 tbox::PIO::logOnlyNodeZero(log_file_name);
      }

      bool dump_databases_to_log = 
	 main_db->getBoolWithDefault("dump_databases_to_log", false);
      bool dump_hierarchy_to_log = 
	 main_db->getBoolWithDefault("dump_hierarchy_to_log", false);
      bool dump_patchdata_to_log = 
	 main_db->getBoolWithDefault("dump_patchdata_to_log", false);

      bool dump_viz_data = 
	 main_db->getBoolWithDefault("dump_viz_data", false);
      string visit_dump_dirname = 
	 main_db->getStringWithDefault("visit_dump_dirname", "visit_data");
      int visit_number_procs_per_file = 
	 main_db->getIntegerWithDefault("visit_number_procs_per_file", 1);

      bool get_total_work_units =
	 main_db->getBoolWithDefault("get_total_work_units", true);

      /*
       * Setup timer manager for profiling code.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      tbox::Pointer<tbox::Timer> t_log_dump = tbox::TimerManager::getManager()->
	 getTimer("apps::main::VisIt & log dump");

      /*
       * Set number of patch data components supported if necessary.  Note this is 
       * kinda cludgy as is, but it seems OK for now.
       */
      int PX = 0;
      int PY = 0;
      int PZ = 0;
      int Pproduct = 0;
      if (input_db->keyExists("PX")) {
	 PX = input_db->getInteger("PX");
	 if (PX > 0) Pproduct = PX;
      }
      if (input_db->keyExists("PY")) {
	 PY = input_db->getInteger("PY");
	 if (PY > 0) Pproduct *= PY;
      }
#if (NDIM > 2)
      if (input_db->keyExists("PZ")) {
	 PZ = input_db->getInteger("PZ");
	 if (PZ > 0) Pproduct *= PZ;
      }
#endif
      int func_factor = 0;
      string func_dist;
      if (input_db->keyExists("FUNC_DIST")) {
	 func_dist = input_db->getString("FUNC_DIST"); 
	 if (func_dist == "FINE_PATCH_EVEN_ODD") {
	    func_factor = 16;
	 }
	 if (func_dist == "FINE_PATCH_MOD2") {
	    func_factor = 32;
	 }
	 if (func_dist == "FINE_PATCH_ALL") {
	    func_factor = 64;
	 }
      }
      int fcn_total = Pproduct * func_factor;
      if (Pproduct * func_factor > 100) {
	 tbox::SAMRAIManager::setMaxNumberPatchDataEntries(fcn_total);
      }

      tbox::pout << "\nConstructing objects..." << endl;
      tbox::pout << "   PX, PY, PZ = " << PX << " , " << PY << " , " << PZ << endl;
      tbox::pout << "   func_dist : func_factor = " << func_dist << " : " << func_factor << endl;
      tbox::pout << "   fcn_total = " << fcn_total << endl;

      // This is here to ensure that the locally-active variable database is
      // the Singleton object that is created.
      (void) hier::LocallyActiveVariableDatabase<NDIM>::getDatabase();

      char test_nproc_buffer[5];
      sprintf(test_nproc_buffer, "%d%s", tbox::SAMRAI_MPI::getNodes(), "\0");

      string nproc_string(test_nproc_buffer);

      string grid_geometry_input_string = string("CartesianGeometry_" + nproc_string);
      if (!input_db->isDatabase(grid_geometry_input_string)) {
	 grid_geometry_input_string = "CartesianGeometry";
      }
      tbox::Pointer< geom::CartesianGridGeometry<NDIM> > geometry = 
	 new geom::CartesianGridGeometry<NDIM>(grid_geometry_input_string,
					       input_db->getDatabase(grid_geometry_input_string));

      bool register_for_restart = false;
      tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy =
	 new hier::PatchHierarchy<NDIM>("PatchHierarchy",
					geometry,
					register_for_restart);

      tbox::Pointer< appu::VisItDataWriter<NDIM> > visit_writer;
      if (dump_viz_data) {
	 visit_writer = 
	    new appu::VisItDataWriter<NDIM>("SAMRES Visit writer",
					    visit_dump_dirname,
					    visit_number_procs_per_file);
      }

      string tester_string = string("LocallyActiveDataTester_" + nproc_string);
      if (!input_db->isDatabase(tester_string)) {
	 tester_string = "LocallyActiveDataTester";
      }
      LocallyActiveDataTester function_manager(
	 "LocallyActiveDataTester",
	 input_db->getDatabase(tester_string),
	 hierarchy,
	 visit_writer);

      if (dump_databases_to_log) {
	 t_log_dump->start();

	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	 tbox::plog << "\nCheck input data and variables before tests:" << endl;
	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

	 tbox::plog << "Input database..." << endl;
	 input_db->printClassData(tbox::plog);

	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

	 tbox::plog << "\nLocally-active Variable<NDIM> database..." << endl;
	 hier::LocallyActiveVariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	 tbox::plog << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
      
	 t_log_dump->stop();
      }

      tbox::Pointer<tbox::Timer> t_run_test = tbox::TimerManager::getManager()->
	 getTimer("apps::main::RunTest");

      t_run_test->start();
      tbox::pout << "\nBuilding patch hierarchy..." << endl;
      function_manager.buildPatchHierarchy();

      tbox::pout << "\nSetting active data on patch hierarchy..." << endl;
      function_manager.setActivePatchesOnHierarchy();
      if (get_total_work_units) {
	 int total_work_units_on_hierarchy = 0;
	 for (int ln = 0; ln < hierarchy->getNumberLevels(); ln++) { 
	    int level_work_units = function_manager.getLevelWorkUnits(ln);
	    tbox::pout << "   work units on level " << ln << "= " << level_work_units << endl;
	    total_work_units_on_hierarchy += level_work_units;
	 }
	 tbox::pout << "Total work units on hierarchy = " << total_work_units_on_hierarchy << endl;
      }

      tbox::pout << "\nInitializing function data..." << endl;
      function_manager.initializeFunctionData();

      int viz_file_count = 0;

      /*
       * Dump all function data before solution.
       */
      t_log_dump->start();
      if (dump_hierarchy_to_log) {
	 tbox::plog << "Printing hierarchy data after initialization..." << endl;
	 function_manager.printHierarchyData(tbox::plog, dump_patchdata_to_log);
      }

      if (dump_viz_data) {
	 tbox::pout << "Writing Visit data to file BEFORE tests..." << endl;
	 visit_writer->writePlotData(hierarchy, viz_file_count);
      }
      t_log_dump->stop();

      tbox::pout << "\nSetting up communication..." << endl;
      function_manager.setupCommunication();
 
      tbox::pout << "Performing test..." << endl;
      function_manager.performTest();

      /*
       * Dump all function data.
       */
      t_log_dump->start();
      if (dump_hierarchy_to_log) {
	 tbox::plog << "\nPrinting hierarchy data after test ..." << endl;
	 function_manager.printHierarchyData(tbox::plog, dump_patchdata_to_log);
      }

      if (dump_viz_data) {
	 viz_file_count++; 
	 tbox::pout << "Writing Visit data to file AFTER tests: " << endl;
	 visit_writer->writePlotData(hierarchy, viz_file_count);
      }
      t_log_dump->stop();

      tbox::pout << "\nChecking test results..." << endl;
      fail_count = function_manager.checkTestResult(tbox::perr);
      t_run_test->stop();

      tbox::TimerManager::getManager()->print(tbox::plog);

      tbox::pout << "\n\n\nDone." << endl;

      if ( fail_count == 0 ) {
	 tbox::pout << "\nPASSED:  locally_active" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(fail_count); 

}

