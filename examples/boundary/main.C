//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/boundary/main.C $
// Package:     SAMRAI example
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2141 $
// Modified:    $LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description: Example program to demonstrate boundary utilities.
//

#include "SAMRAI_config.h"


#include <string>
using namespace std;


#if (NDIM < 2)
   application must be 2d or 3d
#endif

// Headers for basic SAMRAI objects used in this code.
#include "tbox/SAMRAIManager.h"

#include "BoxArray.h"
#include "BoxList.h"
#include "BoxUtilities.h"
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "PatchHierarchy.h"
#include "ProcessorMapping.h"
#include "tbox/InputManager.h"
#include "IntVector.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

// Headers for classes specific to this example
#include "BoundaryDataTester.h"

using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   int fail_count = -1;

   /*
    * Initialize tbox::MPI and SAMRAI, enable logging, and process command line.
    * Note this example is set up to run in serial only.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv); 

   tbox::SAMRAIManager::startup();

   {

   if (argc != 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input filename> "
           << "<restart dir> <restore number> [options]\n"
           << "  options:\n"
           << "  none at this time"
           << endl);
      return (-1);
   }


   /* 
    * This test only is valid on 1 processor.
    * This should be removed.
    */
   if(tbox::SAMRAI_MPI::getNodes() != 1) {
      tbox::pout << "This test is valid for 1 processor only" << endl;
      tbox::pout << "\nPASSED:  boundary" << endl;
      tbox::SAMRAIManager::shutdown();
      tbox::SAMRAI_MPI::finalize();  
      return 0;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * This should never be true.
    */
   TBOX_ASSERT(tbox::SAMRAI_MPI::getNodes() == 1);
#endif


   string input_filename = argv[1];

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
      if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
         bool flag = global_db->
            getBool("call_abort_in_serial_instead_of_exit");
         tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
      }
   }

   /*
    * Read "Main" input data. 
    */

   tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

   string log_file_name = "boundary.log";
   if (main_db->keyExists("log_file_name")) {
      log_file_name = main_db->getString("log_file_name");
   }
   tbox::PIO::logOnlyNodeZero(log_file_name);

   hier::IntVector<NDIM> num_boxes(1);
   if (main_db->keyExists("num_domain_boxes")) {
      int* tmp_arr = num_boxes;
      main_db->getIntegerArray("num_domain_boxes", tmp_arr, NDIM);
   }

   /*
    * Create objects used in boundary data test.  Then, print out
    * state of BoundaryDataTester to log file for checking.
    */

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
      new geom::CartesianGridGeometry<NDIM>("CartesianGridGeometry",
                                input_db->getDatabase("CartesianGridGeometry"));

   tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
      new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   BoundaryDataTester* btester = 
      new BoundaryDataTester("BoundaryDataTester",
                              input_db->getDatabase("BoundaryDataTester"),
                              grid_geometry);

   tbox::plog << "\nPRINTING BoundaryDataTester object state after initialization..." 
        << endl; 
   btester->printClassData(tbox::plog);

   /*
    * For simplicity, we manually create a hierachy with a single patch level. 
    */

   tbox::pout << "\nBuilding patch hierarchy..." << endl;

   const hier::BoxArray<NDIM>& domain = grid_geometry->getPhysicalDomain();
   hier::BoxList<NDIM> boxes(domain);
   if ( (domain.getNumberOfBoxes() == 1) && (num_boxes != hier::IntVector<NDIM>(1)) ) {
      const hier::Box<NDIM>& dbox = domain[0];
      hier::IntVector<NDIM> max_size = dbox.numberCells();
      hier::IntVector<NDIM> min_size = dbox.numberCells() / num_boxes;
      hier::IntVector<NDIM> cut_factor(1);
      hier::IntVector<NDIM> bad_interval(1);
      hier::BoxUtilities<NDIM>::chopBoxes(boxes,
                              max_size,
                              min_size,
                              cut_factor,
                              bad_interval,
                              domain);
   }
   hier::BoxArray<NDIM> patch_boxes(boxes);

   hier::ProcessorMapping mapping(patch_boxes.getNumberOfBoxes());

   for (int ib = 0; ib < patch_boxes.getNumberOfBoxes(); ib++) {
      mapping.setProcessorAssignment(ib, 0);
   }

   patch_hierarchy->makeNewPatchLevel(0, hier::IntVector<NDIM>(1), patch_boxes, mapping);

   /*
    * Allocate data on hierarchy and set variable data on patch interiors
    * to input values.
    */

   tbox::pout << "\nAllocate and initialize data on patch hierarchy..." << endl;

   btester->initializeDataOnPatchInteriors(patch_hierarchy, 0);
   
   tbox::pout << "Performing tests..." << endl;
 
   fail_count = btester->runBoundaryTest(patch_hierarchy, 0);

   tbox::pout << "\n\n\nDone." << endl;

   /*
    * At conclusion of test, deallocate objects.
    */
   patch_hierarchy.setNull();
   grid_geometry.setNull();

   if (btester) delete btester;

   if ( fail_count == 0 ) {
      tbox::pout << "\nPASSED:  boundary test" << endl;
   }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(fail_count);
}


