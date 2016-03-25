//
// File:        main.C
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Main program for hierachy coarsen/refine tests.
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"

#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/MPI.h"
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

   tbox::MPI::init(&argc, &argv);
   SAMRAIManager::startup();

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
         tbox::MPI::setCallAbortInSerialInsteadOfExit(flag);
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

   tbox::plog << "\nCreating initial patch hierarchy..." << endl;
   hierarchy_tester->setupInitialHierarchy(input_db);
   tbox::plog << "\nInitial patch hierarchy created..." << endl;

   tbox::plog << "\nInput file data is ...." << endl;
   input_db->printClassData(plog); 

   tbox::plog << "\nVariable database..." << endl;
   hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog, false);

   tbox::plog << "\nPerforming refine/coarsen patch hierarchy test..." << endl;
   hierarchy_tester->runHierarchyTestAndVerify();
   tbox::plog << "\n Ending hierarchy refine/coarsen test..." << endl;

   tbox::pout << "\nPASSED:  hierarchy tester" << endl;

   SAMRAIManager::shutdown();
   tbox::MPI::finalize();

   return(0); 
}
