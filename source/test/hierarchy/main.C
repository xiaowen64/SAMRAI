//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/hierarchy/main.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Main program for hierachy coarsen/refine tests.
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"

#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include <string>
using namespace std;
#include "VariableDatabase.h"

#include "HierarchyTester.h"

using namespace SAMRAI;

/*
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*/

int main( int argc, char *argv[] )
{

   int fail_count;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {


      string input_filename;

      if (argc != 2)  {
	 TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n" 
		    << "  options:\n"
		    << "  none at this time" << endl);
      } else {
	 input_filename = string(argv[1]);
      }

      tbox::plog << "\n Starting hierarchy refine/coarsen test..." << endl;
      tbox::plog << "Specified input file is: " << input_filename << endl;

      Pointer<Database> input_db = new InputDatabase("input_db");
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "GlobalInputs" section of the input database and set
       * values accordingly.
       */

      if (input_db->keyExists("GlobalInputs")) {
	 tbox::Pointer<tbox::Database> global_db =
	    input_db->getDatabase("GlobalInputs");
	 if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
	    bool flag = global_db->
	       getBool("call_abort_in_serial_instead_of_exit");
	    tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
	 }
      }

      Pointer<Database> main_db = input_db->getDatabase("Main");

      string log_file_name = "hierarchy_test.log";
      if (main_db->keyExists("log_file_name")) {
	 log_file_name = main_db->getString("log_file_name");
      }
      bool log_all_nodes = false;
      if (main_db->keyExists("log_all_nodes")) {
	 log_all_nodes = main_db->getBool("log_all_nodes");
      }
      if (log_all_nodes) {
	 PIO::logAllNodes(log_file_name);
      } else {
	 PIO::logOnlyNodeZero(log_file_name);
      }

      Pointer<HierarchyTester> hierarchy_tester = 
	 new HierarchyTester("HierarchyTester",
			     input_db->getDatabase("HierarchyTest"));

      hierarchy_tester->setupInitialHierarchy(input_db);

      tbox::plog << "\nInput file data is ...." << endl;
      input_db->printClassData(plog); 

      tbox::plog << "\nVariable database..." << endl;
      hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog, false);

      tbox::plog << "\nPerforming refine/coarsen patch hierarchy test..." << endl;
      fail_count = hierarchy_tester->runHierarchyTestAndVerify();
      tbox::plog << "\n Ending hierarchy refine/coarsen test..." << endl;

      if ( fail_count == 0 ) {
	 tbox::pout << "\nPASSED:  hierarchy tester" << endl;
      } 
   }

   SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(fail_count); 
}
