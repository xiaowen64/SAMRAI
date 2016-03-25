//
// File:        main.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Release:     $Name$
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Main program for test of hierarchy sum
//

#include "SAMRAI_config.h"

// Headers for basic SAMRAI objects
#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/MPI.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"

// Headers for major algorithm/data structure objects
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "GriddingAlgorithm.h"
#include "StandardTagAndInitialize.h"
#include "PatchHierarchy.h"
#include "LoadBalancer.h"

// Header for application-specific algorithm/data structure object
#include "HierSumTest.h"

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace geom;
using namespace mesh;
using namespace appu;



/**
 * This is the main program for an example case that tests SAMRAI
 * hierarchy sum classes for FE-type operations with node data, 
 * and level edge sum operations with edge data. 
 *
 * The main program constructs the various SAMRAI gridding objects 
 * and performs the time-stepping loop.  The program should be 
 * executed as:
 *
 *    executable <input file name>
 */

int main( int argc, char *argv[] )
{
   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

   tbox::MPI::init(&argc, &argv);
   SAMRAIManager::startup();

   string input_filename;

   if ( argc != 2 ) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
           << "[options]\n" << endl;
      tbox::MPI::abort();
      return (-1);
   } else {
      input_filename = argv[1];
   }

   tbox::plog << "input_filename = " << input_filename << endl;

   /****************************************************************
    *                                                              *
    *  PROBLEM SETUP                                               *
    *                                                              *
    ****************************************************************
    *                                                              *
    *  Read data from input file and initialize SAMRAI classes     *
    *                                                              *
    ****************************************************************/

   /*
    * Create input database and parse all data in input file. 
    */

   Pointer<Database> input_db = new tbox::InputDatabase("input_db");
   InputManager::getManager()->parseInputFile(input_filename,input_db);

   /*
    * Retrieve "Main" section of the input database. 
    */

   Pointer<Database> main_db = input_db->getDatabase("Main");

   /*
    * Determine if we are doing node sum tests, edge sum tests, 
    * or both.
    */
   bool do_node_sum = false;
   if (main_db->keyExists("do_node_sum")) {
      do_node_sum = main_db->getBool("do_node_sum");
   }
   bool do_edge_sum = false;
   if (main_db->keyExists("do_edge_sum")) {
      do_edge_sum = main_db->getBool("do_edge_sum");
   }

   int nsteps = 1;
   if (main_db->keyExists("nsteps")) {
      nsteps = main_db->getInteger("nsteps");
   }
   
   string log_file_name = "hiersumtest.log";
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

   string visit_dump_dirname = "visit_data";
   int visit_number_procs_per_file = 1;
   if (main_db->keyExists("visit_dump_dirname")) {
      visit_dump_dirname = main_db->getString("visit_dump_dirname");
   }
   if (main_db->keyExists("visit_number_procs_per_file")) {
      visit_number_procs_per_file = 
         main_db->getInteger("visit_number_procs_per_file");
   }

   /*
    * The grid geometry defines the grid type (e.g. cartesian, spherical,
    * etc.).  Because SAMRAI operates on block structured indices, it can
    * support any grid geometry that may be represented as an orthogonal
    * grid.
    */
   Pointer<CartesianGridGeometry<NDIM> > grid_geometry =
      new CartesianGridGeometry<NDIM>("CartesianGeometry",
                                input_db->getDatabase("CartesianGeometry"));

   /*
    * The patch hierarchy defines the adaptive grid system.
    */
   Pointer<PatchHierarchy<NDIM> > patch_hierarchy = 
      new PatchHierarchy<NDIM>("PatchHierarchy<NDIM>", grid_geometry);

   /*
    * Set up Visualization writer.  
    */
   Pointer<VisItDataWriter<NDIM> > visit_data_writer =
      new VisItDataWriter<NDIM>("HierSumTest VisIt Writer",
                          visit_dump_dirname,
                          visit_number_procs_per_file);

   /*
    * This is our problem class.  See the class header for comments on it.
    */
   HierSumTest* hier_sum_test = new HierSumTest(
      "HierSumTest",
      input_db->getDatabase("HierSumTest"),
      visit_data_writer);

   /*
    * The StandardTagAndInitialize<NDIM> class performs a variety of operations
    * with user-specified parameters related to adptive gridding.  For example,
    * it manages initialization of a level, cell tagging using a gradient
    * detector, and methods to reset data after the hierarchy has been 
    * regridded.  
    */
   Pointer<StandardTagAndInitialize<NDIM> > tag_and_init_ops =
      new StandardTagAndInitialize<NDIM>(
         "StandardTagAndInitialize",
         hier_sum_test,
         input_db->getDatabase("StandardTagAndInitialize"));

   /*
    * The gridding algorithm manages adaptive gridding.  It expects a 
    * clustering scheme (i.e. how to cluster tagged-cells into patches),
    * and a load balance scheme to distribute work to processors.  In general
    * the baseline classes provided in SAMRAI should suffice for most
    * problems. It also requires a class that defines the particular tag 
    * and initialization ops that correlate with the users problem.  For
    * this, we use the "tag_and_init_ops" above, which references our 
    * "wave_eqn_model" problem class to define the user-specific operations.
    */
   Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

   Pointer<LoadBalancer<NDIM> > load_balancer = 
      new LoadBalancer<NDIM>("LoadBalancer", 
                             input_db->getDatabase("LoadBalancer"));

   Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = 
      new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            input_db->getDatabase("GriddingAlgorithm"),
                            tag_and_init_ops,
                            box_generator,
                            load_balancer);

   /*
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents to
    * the log file.
    */

   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(plog);
   tbox::plog << "\nVariable database..." << endl;
   VariableDatabase<NDIM>::getDatabase()->printClassData(plog);

   

   /****************************************************************
    *                                                              *
    *  INITIALIZE DATA ON PATCHES                                  *
    *                                                              *
    ****************************************************************
    *                                                              *
    *  Build patch hierarchy and initialize the data on the patches*
    *  in the hierarchy.
    *  1) Create a "tag_buffer" for each level in the Hierarchy.   *
    *  2) Create the coarse (i.e. level 0) grid.                   *
    *  3) Cycle through levels 1-max_levels, initializing data     *
    *     on each.  The makeFinerLevel method calls the error      *
    *     estimator (remember, it was registered with the          *
    *     gridding algorithm object) and tags cells for refinement *
    *     as it generates patches on the finer levels.             *
    *                                                              *
    ****************************************************************/


   double loop_time = 0.;
   tbox::Array<int> tag_buffer_array(gridding_algorithm->getMaxLevels());
   for (int il = 0; il < gridding_algorithm->getMaxLevels(); il++) {
      tag_buffer_array[il] = 1;
   }
   gridding_algorithm->makeCoarsestLevel(patch_hierarchy,loop_time);
   
   bool done = false;
   bool initial_time = true;
   for (int ln = 0; 
        gridding_algorithm->levelCanBeRefined(ln) && !done; 
        ln++) {
      gridding_algorithm->makeFinerLevel(patch_hierarchy,
                                         loop_time,
                                         initial_time,
                                         tag_buffer_array[ln]);
      done = !(patch_hierarchy->finerLevelExists(ln));
   }

   int nlevels = patch_hierarchy->getNumberLevels();

   for (int pln = 0; pln <= patch_hierarchy->getFinestLevelNumber(); pln++) {
      Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(pln);

      tbox::plog << "\n PRINTING PATCHES ON LEVEL " << pln << endl;

      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         tbox::plog << "patch # " << ip() << " : " 
              << level->getPatch(ip())->getBox() << endl;
      }
   }
   tbox::plog << endl;
 

   /*******************************************************************
    *                                                                 *
    * Test hier sum operation                                         *
    *                                                                 *
    *******************************************************************
    *                                                                 *
    *  1) Set node values (initial)                                   *
    *  2) Do hierarchy sum operation                                  *
    *  3) Check result                                                *
    *                                                                 *
    ******************************************************************/

   /*
    * Setup the sum operation(s)
    */
   if (do_node_sum) {
      hier_sum_test->setupOuternodeSum(patch_hierarchy);
   }
   if (do_edge_sum) {
      for (int ln = 0; ln < nlevels; ln++) {
         hier_sum_test->setupOuteredgeSum(patch_hierarchy,
                                          ln);
      }
   }

   for (int i = 0; i < nsteps; i++) {

      /*
       * In the process of constructing the hierarchy, we set cell values
       * to their proper weights.  Now go in and set the node/edge values. 
       * (write data to VisIt once it is set).
       */
      if (do_node_sum) {
         hier_sum_test->setInitialNodeValues(patch_hierarchy);
      }
      if (do_edge_sum) {
         for (int ln = 0; ln < nlevels; ln++) {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            hier_sum_test->setInitialEdgeValues(level);
         }
      }

      /*
       * Write the pre-summed cell/node data to VisIt
       */
      visit_data_writer->writePlotData(patch_hierarchy,i,loop_time);

      /*
       * Perform the sum operation(s)
       */
      if (do_node_sum) {
         hier_sum_test->doOuternodeSum();
      }
      if (do_edge_sum) {
         for (int ln = 0; ln < nlevels; ln++) {
            hier_sum_test->doOuteredgeSum(ln);
         }
      }
   }
      
   /*
    * Check result
    */
   if (do_node_sum) {
      hier_sum_test->checkNodeResult(patch_hierarchy);
   }

   tbox::pout << "\n" << endl;
   
   if (do_edge_sum) {
      for (int ln = 0; ln < nlevels; ln++) {
         Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
         hier_sum_test->checkEdgeResult(level);
      }
   }

   /*
    * Write the post-summed cell/node data to VisIt
    */
   visit_data_writer->writePlotData(patch_hierarchy,nsteps+1,loop_time);
   
   /*
    * At conclusion of simulation, deallocate objects.
    */
   visit_data_writer.setNull();
   gridding_algorithm.setNull();
   load_balancer.setNull();
   box_generator.setNull(); 
   tag_and_init_ops.setNull();

   if (hier_sum_test) delete hier_sum_test;

   patch_hierarchy.setNull();
   grid_geometry.setNull();

   input_db.setNull();
   main_db.setNull();
  
   tbox::pout << "\nPASSED:  patchbdrysum" << endl;
 
   SAMRAIManager::shutdown();
   tbox::MPI::finalize();
 
   return(0); 
}



