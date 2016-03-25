//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/mblkcomm/main.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Main program for patch data communication tests.
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"

#include "MultiblockTester.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "StandardTagAndInitialize.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"

// Different component tests available
#include "CellMultiblockTest.h"
#include "EdgeMultiblockTest.h"
#include "FaceMultiblockTest.h"
#include "NodeMultiblockTest.h"
#include "SideMultiblockTest.h"

using namespace SAMRAI;

/*
************************************************************************
*                                                                      *
* This is the driver program to test and time patch data               *
* communication operations on an SAMR patch hierarchy using            *
* SAMRAI.  MultiblockTester is the primary object used in these              *
* processes.  It constructs the patch hierarchy based on               * 
* input file information and invokes the communcation operations       *
* specified in the input file.  The implementation of data type        *
* specific operations (defining variables, initializing data,          *
* defining coarsen/refine operations, and verifying the results)       *
* are provided in a class implemented for the test to be performed.    *
* This test-specific class is derived from the PatchMultiblockTestStrategy   *
* base class which declares the interface between the MultiblockTester       *
* and the test.                                                        *
*                                                                      *
* Input data file sections and keys are defined as follows:            *
*                                                                      *
*    o Main program...                                                 *
*                                                                      *
*      Main {                                                          *
*         log_file_name  = <string> [name of log file]                 *
*                          (optional - "component_test.log" is default)*
*         log_all_nodes  = <bool> [log all nodes or node 0 only?]      *
*                          (optional - FALSE is default)               *
*         ntimes_run     = <int> [how many times to perform test]      *
*                          (optional - 1 is default)                   *
*         test_to_run    = <string> [name of test] (required)          *
*            Available tests are:                                      *
*               "CellMultiblockTest"                                         *
*               "EdgeMultiblockTest"                                         *
*               "FaceMultiblockTest"                                         *
*               "NodeMultiblockTest"                                         *
*               "OuterodeMultiblockTest"                                     *
*               "SideMultiblockTest"                                         *
*               "MultiVariableMultiblockTest"                                *
*         do_refine      = <bool> [test refine operation?]             *
*                          (optional - FALSE is default)               *
*         do_coarsen     = <bool> [test coarsen operation?]            *
*                          (optional - FALSE is default)               *
*         NOTE: Only refine or coarsen test can be run, but not both.  *
*               If both are TRUE, only refine operations will execute. *
*         refine_option  = <string> [how interior of destination       *
*                                    level is filled during refine]    *
*            Options are:                                              *
*               "INTERIOR_FROM_SAME_LEVEL"                             *
*               "INTERIOR_FROM_COARSER_LEVEL"                          *
*               (default is "INTERIOR_FROM_SAME_LEVEL")                *
*      }                                                               *
*                                                                      *
*    o Timers...                                                       *
*                                                                      *
*      TimerManager {                                                  *
*         timer_list = <string array> [names of timers to run]         *
*            Available timers are:                                     *
*               "test::main::createRefineSchedule"                     *
*               "test::main::performRefineOperations"                  *
*               "test::main::createCoarsenSchedule"                    *
*               "test::main::performCoarsenOperations"                 *
*      }                                                               *
*                                                                      *
*    o Patch data tests...                                             *
*                                                                      *
*      Each test defines the input parameters it needs.  Consult the   *
*      documentation in each class for details.  Default operations    *
*      for reading variable data and mesh refinement information       *
*      for data tests are provided in the PatchMultiblockTestStrategy        *
*      base class.  These input are typically read from the input      *
*      file section for each test.  In this case, the input data       *
*      keys are similar to the following example:                      *
*                                                                      *
*      VariableData {  // The variable sub-database                    *
*                      // (key name VariableData not optional)         *
*                                                                      *
*         variable_1 { // sub-database for first variable              *
*                      // (Key name for each variable can be anything. *
*                      //  Key names for variable parameters are       *
*                      //  not optional. However, only name data is    *
*                      //  required)                                   *
*            name = "var1"    // <string> variable name (required)     *
*            depth = 1        // <int> variable depth (opt. - def is 1)*
*            src_ghosts = 0,0,0 // <int array> for ghost width of      *
*                                  source data (opt. - def is 0,0,0)   *
*            dst_ghosts = 1,1,1 // <int array> for ghost width of      *
*                                  dest data (opt. - def is 0,0,0)     *
*            coarsen_operator = "CONSERVATIVE_COARSEN"                 *
*            refine_operator = "LINEAR_REFINE"                         *
*            // Interlevel transfer operator name strings are optional *
*            // Default are "NO_COARSEN", and "NO_REFINE", resp.       * 
*         }                                                            *
*                                                                      *
*         // data for other variables as needed...                     *
*                                                                      *
*      }                                                               *
*                                                                      *
*      RefinementData {  // The variable sub-database                  *
*                        // (key name RefinementData not optional)     *
*                                                                      *
*         // Lists of boxes to refine on each level.  Names of box     *
*         // arrays may be anything.  For example,                     * 
*                                                                      *
*           level0_boxes = [ (1,2,3) , (3,3,5) ]                       *
*           level1_boxes = [ (8,10,6) , (10,10,12) ]                   *
*           // other level box information as needed...                *
*      }                                                               *
*                                                                      *
*  NOTES:                                                              *
*                                                                      *
*     o The MultiblockTester uses the GriddingAlgorithm, and                 *
*       LoadBalancer class to construct the patch hierarchy.           *
*       Appropriate input sections must be provided for these objects  *
*       as needed.                                                     *
*                                                                      *
*     o Each test must register a GridGeometry object with the         *
*       PatchMultiblockTestStrategy base class so the hierarchy can be       *
*       constructed.  Consult the constructor of each test class       *
*       for inforamation about which geomteyr object is constructed,   *
*       and thus which input data is required to initialize the geom.  * 
*                                                                      *
************************************************************************
*/

int main( int argc, char *argv[] )
{

   int return_val = 1;

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
       * Process command line arguments.  For each run, the input 
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       *
       */
      string input_filename;

      if (argc != 2)  {
	 TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n" 
		    << "  options:\n"
		    << "  none at this time" << endl);
      } else {
	 input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */

      tbox::Pointer<tbox::Database> input_db =
	 new tbox::InputDatabase("input_db");
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section from input database.  Set log file 
       * parameters, number of times to run tests (for performance
       * analysis), and read in test information.
       */

      tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

      string log_file_name = "component_test.log";
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

      int ntimes_run = 1;
      if (main_db->keyExists("ntimes_run")) {
	 ntimes_run = main_db->getInteger("ntimes_run");
      }

      string test_to_run;
      if (main_db->keyExists("test_to_run")) {
	 test_to_run = main_db->getString("test_to_run");
      } else {
	 TBOX_ERROR("Error in Main input: no test specified." << endl);
      }
   
      bool do_refine = true;  
      bool do_coarsen = false;

      string refine_option = "INTERIOR_FROM_SAME_LEVEL";

      tbox::plog << "\nPerforming refine data test..." << endl;
      if (main_db->keyExists("refine_option")) {
	 refine_option = main_db->getString("refine_option");
      }

      /*
       * Create communication tester and patch data test object
       */
   
      PatchMultiblockTestStrategy* patch_data_test = NULL;

      if (test_to_run == "CellMultiblockTest") {
	 patch_data_test = new CellMultiblockTest("CellMultiblockTest",
						  input_db, 
						  do_refine,
						  do_coarsen,
						  refine_option);

      } else if (test_to_run == "EdgeMultiblockTest") {
	 patch_data_test = new EdgeMultiblockTest("EdgeMultiblockTest",
						  input_db,
						  do_refine,
						  do_coarsen,
						  refine_option);
      } else if (test_to_run == "FaceMultiblockTest") {
	 patch_data_test = new FaceMultiblockTest("FaceMultiblockTest",
						  input_db,
						  do_refine,
						  do_coarsen,
						  refine_option);
      } else if (test_to_run == "NodeMultiblockTest") {
	 patch_data_test = new NodeMultiblockTest("NodeMultiblockTest",
						  input_db,
						  do_refine,
						  do_coarsen,
						  refine_option); 
      } else if (test_to_run == "SideMultiblockTest") {
	 patch_data_test = new SideMultiblockTest("SideMultiblockTest",
						  input_db,
						  do_refine,
						  do_coarsen,
						  refine_option);
      } else if (test_to_run == "MultiVariableMultiblockTest") {
	 TBOX_ERROR("Error in Main input: no multi-variable test yet." << endl);
      } else {
	 TBOX_ERROR("Error in Main input: illegal test = " <<
		    test_to_run << endl);
      } 

      tbox::Pointer<MultiblockTester> comm_tester = new MultiblockTester("MultiblockTester",
									 input_db,
									 patch_data_test,
									 do_refine,
									 do_coarsen, 
									 refine_option); 

      tbox::Pointer <mesh::StandardTagAndInitialize<NDIM> > cell_tagger =
	 new mesh::StandardTagAndInitialize<NDIM>(
	    "StandardTagggingAndInitializer",
	    comm_tester,
	    input_db->getDatabase("StandardTaggingAndInitializer"));

      comm_tester->setupHierarchy(input_db, cell_tagger);

      tbox::plog << "Specified input file is: " << input_filename << endl;

      tbox::plog << "\nInput file data is ...." << endl;
      input_db->printClassData(tbox::plog); 
   
      tbox::plog << "\nCheck Variable database..." << endl;
      hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

      /*
       * Create timers from input data to check performance of comm. operations.
       */
   
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

      tbox::TimerManager* time_man = tbox::TimerManager::getManager();

      tbox::Pointer<tbox::Timer> refine_create_time = 
	 time_man->getTimer("test::main::createRefineSchedule");
      tbox::Pointer<tbox::Timer> refine_comm_time = 
	 time_man->getTimer("test::main::performRefineOperations");

      tbox::Pointer<tbox::Timer> coarsen_create_time = 
	 time_man->getTimer("test::main::createCoarsenSchedule");
      tbox::Pointer<tbox::Timer> coarsen_comm_time =
	 time_man->getTimer("test::main::performCoarsenOperations"); 

      tbox::TimerManager::getManager()->resetAllTimers();
   
      /*
       * Create communication schedules and perform communication operations.
       */ 

      tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> > patch_hierarchy =
	 comm_tester->getPatchHierarchy();

      int nlevels = 0;

      for (int nb = 0; nb < patch_hierarchy->getNumberOfBlocks(); nb++) {
	 int block_levels = patch_hierarchy->getHierarchy(nb)->getNumberOfLevels();
	 if (block_levels > nlevels) {
	    nlevels = block_levels;
	 }
      }

      for (int n = 0; n < ntimes_run; n++) {

	 /*
	  * Create communication schedules for data refine tests.
	  */
	 refine_create_time->start();
	 for (int i = 0; i < nlevels; i++) {
	    comm_tester->createRefineSchedule(i);
	 }
	 refine_create_time->stop();

	 /*
	  * Perform refine data communication operations.
	  */
	 refine_comm_time->start();
	 for (int j = 0; j < nlevels; j++) {
	    comm_tester->performRefineOperations(j);
	 }
	 refine_comm_time->stop();

      }

      bool test_passed = comm_tester->verifyCommunicationResults();

      /*
       * Deallocate objects when done.
       */

      if (patch_data_test) delete patch_data_test;

      comm_tester.setNull();

      tbox::TimerManager::getManager()->print(tbox::plog);


      if (test_passed) {
	 tbox::pout << "\nPASSED:  mblkcomm" << endl;
	 return_val = 0;
      }

   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(return_val); 
}
