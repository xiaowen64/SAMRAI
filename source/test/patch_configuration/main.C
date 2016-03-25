//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/patch_configuration/main.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Tests patch configuration utilities
//

#include "SAMRAI_config.h"


/*
 * Header files for classes particular to this application
 */
#include "PatchConfigurationTester.h"

/*
 * Header files for SAMRAI library classes
 */
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/IOStream.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAIManager.h"

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

      bool dump_hierarchy_to_log = 
	 main_db->getBoolWithDefault("dump_hierarchy_to_log", false);

      tbox::pout << "\nConstructing objects..." << endl;

      string nproc_string( tbox::Utilities::intToString(tbox::SAMRAI_MPI::getNodes(), 5) );

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

      string tester_string = string("PatchConfigurationTester" + nproc_string);
      if (!input_db->isDatabase(tester_string)) {
	 tester_string = "PatchConfigurationTester";
      }
      PatchConfigurationTester patch_tester(
	 "PatchConfigurationTester",
	 input_db->getDatabase(tester_string),
	 hierarchy);

      tbox::pout << "\nBuilding patch hierarchy..." << endl;
      patch_tester.buildPatchHierarchy();

      if (dump_hierarchy_to_log) {
	 tbox::plog << "Printing hierarchy after initialization..." << endl;
	 patch_tester.printHierarchyData(tbox::plog);
      }

      tbox::pout << "\nBuilding patch configuration data..." << endl;
      patch_tester.setupPatchConfiguration();

      tbox::pout << "Performing tests..." << endl;
      fail_count = patch_tester.checkPatchConfiguration(tbox::plog);

      tbox::pout << "\n\n\nDone." << endl;

      if ( fail_count == 0 ) {
	 tbox::pout << "\nPASSED:  patch configuration" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(fail_count);
}

