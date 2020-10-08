/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2020 Lawrence Livermore National Security, LLC
 * Description:   Main program for patch data communication tests.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <string>
#include <memory>

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "CommTester.h"
#include "test/testlib/DerivedVisOwnerData.h"

#include "SAMRAI/hier/BlueprintUtils.h"
#include "SAMRAI/tbox/ConduitDatabase.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"

// Different component tests available
#include "CellDataTest.h"
#include "EdgeDataTest.h"
#include "FaceDataTest.h"
#include "NodeDataTest.h"
#include "OuternodeDataTest.h"
#include "SideDataTest.h"
#include "OutersideDataTest.h"
#include "OuterfaceDataTest.h"
//#include "MultiVariableDataTest.h"


#ifdef SAMRAI_HAVE_CONDUIT
#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mesh.hpp"
#include "conduit_blueprint_mpi_mesh.hpp"
#include "conduit_relay.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace SAMRAI;

void makeQuadElemConnectivity(
   std::vector<int64_t>& connect,
   int64_t element,
   int64_t iwidth);

void transformToPolyhedral(conduit::Node& out, const conduit::Node& in, const std::string& name);
/*
 ************************************************************************
 *
 * This is the driver program to test and time patch data
 * communication operations on an SAMR patch hierarchy using
 * SAMRAI.  CommTester is the primary object used in these
 * processes.  It constructs the patch hierarchy based on
 * input file information and invokes the communcation operations
 * specified in the input file.  The implementation of data type
 * specific operations (defining variables, initializing data,
 * defining coarsen/refine operations, and verifying the results)
 * are provided in a class implemented for the test to be performed.
 * This test-specific class is derived from the PatchDataTestStrategy
 * base class which declares the interface between the CommTester
 * and the test.
 *
 * Input data file sections and keys are defined as follows:
 *
 *    o Main program...
 *
 *      Main {
 *         log_file_name  = <string> [name of log file]
 *                          (optional - "component_test.log" is default)
 *         log_all_nodes  = <bool> [log all nodes or node 0 only?]
 *                          (optional - FALSE is default)
 *         ntimes_run     = <int> [how many times to perform test]
 *                          (optional - 1 is default)
 *         test_to_run    = <string> [name of test] (required)
 *            Available tests are:
 *               "CellDataTest"
 *               "EdgeDataTest"
 *               "FaceDataTest"
 *               "NodeDataTest"
 *               "OuterodeDataTest"
 *               "SideDataTest"
 *               "MultiVariableDataTest"
 *         do_refine      = <bool> [test refine operation?]
 *                          (optional - FALSE is default)
 *         do_coarsen     = <bool> [test coarsen operation?]
 *                          (optional - FALSE is default)
 *         NOTE: Only refine or coarsen test can be run, but not both.
 *               If both are TRUE, only refine operations will execute.
 *         refine_option  = <string> [how interior of destination
 *                                    level is filled during refine]
 *            Options are:
 *               "INTERIOR_FROM_SAME_LEVEL"
 *               "INTERIOR_FROM_COARSER_LEVEL"
 *               (default is "INTERIOR_FROM_SAME_LEVEL")
 *      }
 *
 *    o Timers...
 *
 *      tbox::TimerManager {
 *         timer_list = <string array> [names of timers to run]
 *            Available timers are:
 *               "test::main::createRefineSchedule"
 *               "test::main::performRefineOperations"
 *               "test::main::createCoarsenSchedule"
 *               "test::main::performCoarsenOperations"
 *      }
 *
 *    o hier::Patch data tests...
 *
 *      Each test defines the input parameters it needs.  Consult the
 *      documentation in each class for details.  Default operations
 *      for reading variable data and mesh refinement information
 *      for data tests are provided in the PatchDataTestStrategy
 *      base class.  These input are typically read from the input
 *      file section for each test.  In this case, the input data
 *      keys are similar to the following example:
 *
 *      VariableData {  // The variable sub-database
 *                      // (key name VariableData not optional)
 *
 *         variable_1 { // sub-database for first variable
 *                      // (Key name for each variable can be anything.
 *                      //  Key names for variable parameters are
 *                      //  not optional. However, only name data is
 *                      //  required)
 *            name = "var1"    // <string> variable name (required)
 *            depth = 1        // <int> variable depth (opt. - def is 1)
 *            src_ghosts = 0,0,0 // <int array> for ghost width of
 *                                  source data (opt. - def is 0,0,0)
 *            dst_ghosts = 1,1,1 // <int array> for ghost width of
 *                                  dest data (opt. - def is 0,0,0)
 *            coarsen_operator = "CONSERVATIVE_COARSEN"
 *            refine_operator = "LINEAR_REFINE"
 *            // Interlevel transfer operator name strings are optional
 *            // Default are "NO_COARSEN", and "NO_REFINE", resp.
 *         }
 *
 *         // data for other variables as needed...
 *
 *      }
 *
 *      RefinementData {  // The variable sub-database
 *                        // (key name RefinementData not optional)
 *
 *         // Lists of boxes to refine on each level.  Names of box
 *         // arrays may be anything.  For example,
 *
 *           level0_boxes = [ (1,2,3) , (3,3,5) ]
 *           level1_boxes = [ (8,10,6) , (10,10,12) ]
 *           // other level box information as needed...
 *      }
 *
 *  NOTES:
 *
 *     o The CommTester uses the mesh::GriddingAlgorithm, and
 *       mesh::LoadBalancer class to construct the patch hierarchy
 *       Appropriate input sections must be provided for these objects
 *       as needed.
 *
 *     o Each test must register a hier::BaseGridGeometry object
 *       with the
 *       PatchDataTestStrategy base class so the hierarchy can be
 *       constructed.  Consult the constructor of each test class
 *       for inforamation about which geomteyr object is constructed,
 *       and thus which input data is required to initialize the geom.
 *
 ************************************************************************
 */

int main(
   int argc,
   char* argv[])
{
   int return_val = 1;

   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Make block to force Pointers to be deallocated to prevent memory
    * leaks.
    */
   {

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       *
       */
      std::string input_filename;

      if (argc != 2) {
         TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                               << "  options:\n"
                               << "  none at this time" << std::endl);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Create timers from input data to check performance of comm. operations.
       */

      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

      /*
       * Retrieve "GlobalInputs" section of the input database and set
       * values accordingly.
       */

      if (input_db->keyExists("GlobalInputs")) {
         std::shared_ptr<tbox::Database> global_db(
            input_db->getDatabase("GlobalInputs"));
         if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
            bool flag = global_db->
               getBool("call_abort_in_serial_instead_of_exit");
            tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
         }
      }

      /*
       * Retrieve "Main" section from input database.  Set log file
       * parameters, number of times to run tests (for performance
       * analysis), and read in test information.
       */

      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      const std::string base_name =
         main_db->getStringWithDefault("base_name", "component_test");

      std::string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      if (main_db->keyExists("log_all_nodes")) {
         log_all_nodes = main_db->getBool("log_all_nodes");
      }
      if (log_all_nodes) {
         tbox::PIO::logAllNodes(log_file_name);
      } else {
         tbox::PIO::logOnlyNodeZero(log_file_name);
      }

#ifdef _OPENMP
      tbox::plog << "Compiled with OpenMP version " << _OPENMP
                 << ".  Running with " << omp_get_max_threads() << " threads."
                 << std::endl;
#else
      tbox::plog << "Compiled without OpenMP.\n";
#endif

      int ntimes_run = 1;
      if (main_db->keyExists("ntimes_run")) {
         ntimes_run = main_db->getInteger("ntimes_run");
      }

      std::string test_to_run;
      if (main_db->keyExists("test_to_run")) {
         test_to_run = main_db->getString("test_to_run");
      } else {
         TBOX_ERROR("Error in Main input: no test specified." << std::endl);
      }

      bool do_refine = false;
      bool do_coarsen = false;
      std::string refine_option = "INTERIOR_FROM_SAME_LEVEL";
      if (main_db->keyExists("do_refine")) {
         do_refine = main_db->getBool("do_refine");
         if (do_refine) {
            tbox::plog << "\nPerforming refine data test..." << std::endl;
            if (main_db->keyExists("refine_option")) {
               refine_option = main_db->getString("refine_option");
            }
            tbox::plog << "\nRefine data option = " << refine_option << std::endl;

         }
      }

      if (!do_refine) {
         if (main_db->keyExists("do_coarsen")) {
            do_coarsen = main_db->getBool("do_coarsen");
         }
         if (do_coarsen) {
            tbox::plog << "\nPerforming coarsen data test..." << std::endl;
         }
      }

      /*
       * Create communication tester and patch data test object
       */

      PatchDataTestStrategy* patch_data_test = 0;

      if (test_to_run == "CellDataTest") {
         patch_data_test = new CellDataTest("CellDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);

      } else if (test_to_run == "EdgeDataTest") {
         patch_data_test = new EdgeDataTest("EdgeDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "FaceDataTest") {
         patch_data_test = new FaceDataTest("FaceDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "OuterfaceDataTest") {
         patch_data_test = new OuterfaceDataTest("OuterfaceDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "NodeDataTest") {
         patch_data_test = new NodeDataTest("NodeDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "OuternodeDataTest") {
         patch_data_test = new OuternodeDataTest("OuternodeDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "SideDataTest") {
         patch_data_test = new SideDataTest("SideDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "OutersideDataTest") {
         patch_data_test = new OutersideDataTest("OutersideDataTest",
               dim,
               input_db,
               do_refine,
               do_coarsen,
               refine_option);
      } else if (test_to_run == "MultiVariableDataTest") {
         TBOX_ERROR("Error in Main input: no multi-variable test yet." << std::endl);
      } else {
         TBOX_ERROR(
            "Error in Main input: illegal test = " << test_to_run << std::endl);
      }

      std::shared_ptr<CommTester> comm_tester(
         new CommTester(
            "CommTester",
            dim,
            input_db,
            patch_data_test,
            do_refine,
            do_coarsen,
            refine_option));

      std::shared_ptr<mesh::StandardTagAndInitialize> cell_tagger(
         new mesh::StandardTagAndInitialize(
            "StandardTaggingAndInitializer",
            comm_tester.get(),
            input_db->getDatabase("StandardTaggingAndInitializer")));

      comm_tester->setupHierarchy(input_db, cell_tagger);

      tbox::plog << "Specified input file is: " << input_filename << std::endl;

      tbox::plog << "\nInput file data is ...." << std::endl;
      input_db->printClassData(tbox::plog);

      tbox::plog << "\nCheck hier::Variable database..." << std::endl;
      hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

      tbox::TimerManager* time_man = tbox::TimerManager::getManager();

      std::shared_ptr<tbox::Timer> refine_create_time(
         time_man->getTimer("test::main::createRefineSchedule"));
      std::shared_ptr<tbox::Timer> refine_comm_time(
         time_man->getTimer("test::main::performRefineOperations"));

      std::shared_ptr<tbox::Timer> coarsen_create_time(
         time_man->getTimer("test::main::createCoarsenSchedule"));
      std::shared_ptr<tbox::Timer> coarsen_comm_time(
         time_man->getTimer("test::main::performCoarsenOperations"));

      const bool plot = main_db->getBoolWithDefault("plot", false);
      DerivedVisOwnerData vdd;
      if (plot) {
#ifdef HAVE_HDF5
         const std::string visit_filename = base_name + ".visit";
         /* Create the VisIt data writer. */
         std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
            new appu::VisItDataWriter(
               dim,
               "VisIt Writer",
               visit_filename));
         /*
          * The VisItDataWriter requires some value to be plotted.
          * We are registering the owner value just so we can plot.
          */
         visit_data_writer->registerDerivedPlotQuantity("Owner", "SCALAR", &vdd);
         /* Write the plot file. */
         visit_data_writer->writePlotData(
            comm_tester->getPatchHierarchy(), 0);
#else
         TBOX_WARNING("Cannot write VisIt file--not configured with HDF5.");
#endif
      }

      tbox::TimerManager::getManager()->resetAllTimers();

      /*
       * Create communication schedules and perform communication operations.
       */

      std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         comm_tester->getPatchHierarchy());
      patch_hierarchy->recursivePrint(tbox::plog,
         "H-> ",
         3);
      const int nlevels = patch_hierarchy->getNumberOfLevels();

      if (do_refine) {

         for (int n = 0; n < ntimes_run; ++n) {

            /*
             * Create communication schedules for data refine tests.
             */
            refine_create_time->start();
            for (int i = 0; i < nlevels; ++i) {
               comm_tester->createRefineSchedule(i);
            }
            refine_create_time->stop();

            /*
             * Perform refine data communication operations.
             */
            refine_comm_time->start();
            for (int j = 0; j < nlevels; ++j) {
               comm_tester->performRefineOperations(j);
            }
            refine_comm_time->stop();

         }

      }

      if (do_coarsen) {

         for (int n = 0; n < ntimes_run; ++n) {

            /*
             * Create communication schedules for data coarsen tests.
             */
            coarsen_create_time->start();
            for (int i = nlevels - 1; i > 0; --i) {
               comm_tester->createCoarsenSchedule(i);
            }
            coarsen_create_time->stop();

            /*
             * Perform coarsen data communication operations.
             */
            coarsen_comm_time->start();
            for (int j = nlevels - 1; j > 0; --j) {
               comm_tester->performCoarsenOperations(j);
            }
            coarsen_comm_time->stop();

         }

      }

      bool composite_test_passed = true;
      if (do_refine) {
         for (int i = 0; i < nlevels; ++i) {
            composite_test_passed = comm_tester->performCompositeBoundaryComm(i);
         }
      }

#ifdef SAMRAI_HAVE_CONDUIT
      std::shared_ptr<tbox::MemoryDatabase> memory_db(
         new tbox::MemoryDatabase("mem_hierarchy"));

      patch_hierarchy->putToRestart(memory_db);

      std::shared_ptr<tbox::ConduitDatabase> conduit_db(
         new tbox::ConduitDatabase("conduit_hierarchy"));

      std::shared_ptr<tbox::ConduitDatabase> flat_db(
         new tbox::ConduitDatabase("flat_hierarchy"));

      hier::BlueprintUtils bp_utils(comm_tester.get(), true);
      patch_hierarchy->makeBlueprintDatabase(conduit_db, bp_utils);
      patch_hierarchy->makeFlattenedBlueprintDatabase(flat_db, bp_utils);

      conduit::Node bp_node;
      conduit_db->toConduitNode(bp_node);
      conduit::Node flatn;
      flat_db->toConduitNode(flatn);

      std::vector<int> first_patch_id;
      first_patch_id.push_back(0);

      int patch_count = 0;
      for (int i = 1; i <  patch_hierarchy->getNumberOfLevels(); ++i) {
         patch_count += patch_hierarchy->getPatchLevel(i-1)->getNumberOfPatches();
         first_patch_id.push_back(patch_count);
      }

      int num_hier_patches = 0; 
      for (int i = 0; i < patch_hierarchy->getNumberOfLevels(); ++i) {
         const std::shared_ptr<hier::PatchLevel>& level =  patch_hierarchy->getPatchLevel(i);
         num_hier_patches += patch_hierarchy->getPatchLevel(i)->getNumberOfPatches();

         for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
              ++p) {

            const std::shared_ptr<hier::Patch>& patch = *p;
            const hier::BoxId& box_id = patch->getBox().getBoxId();
            const hier::LocalId& local_id = box_id.getLocalId();
  
            int mesh_id = first_patch_id[i] + local_id.getValue();

            if (test_to_run == "CellDataTest") {
               CellDataTest* cell_test = (CellDataTest*)patch_data_test;

               cell_test->addFields(bp_node, mesh_id, patch); 
            } else if (test_to_run == "NodeDataTest") {
               NodeDataTest* node_test = (NodeDataTest*)patch_data_test;

               node_test->addFields(bp_node, mesh_id, patch);
            }
         }
      }

      if (test_to_run == "CellDataTest") {
         bp_utils.flattenFields(flatn, bp_node);
      }

      bp_utils.writeBlueprintMesh(
         bp_node,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         num_hier_patches,
         "amr_mesh",
         "celldata",
         "bpindex.root",
         "json");
/*
      bp_utils.writeBlueprintMesh(
         flatn,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         num_hier_patches,
         "flat_mesh",
         "flat",
         "flat.root",
         "json");
*/
      conduit::Node info;
      TBOX_ASSERT(conduit::blueprint::verify("mesh", bp_node, info));
      TBOX_ASSERT(conduit::blueprint::verify("mesh", flatn, info));

      conduit::Node poly;

      if (dim.getValue() == 2) {
         conduit::blueprint::mpi::mesh::to_poly(flatn, poly, "mesh");
      } else {
         conduit::blueprint::mpi::mesh::to_polyhedral(flatn, poly, "mesh");
      }
      //transformToPolyhedral(poly, flatn, "mesh");
      TBOX_ASSERT(conduit::blueprint::verify("mesh", poly, info));

      bp_utils.writeBlueprintMesh(
         poly,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         num_hier_patches,
         "mesh",
         "polyhedral",
         "poly.root",
         "json");

/*
      conduit::Node polyuni = poly;

      conduit::NodeConstIterator itr = poly.children();

      conduit::Node polystruct = polyuni;

      itr = polyuni.children();
      while(itr.has_next())
      {
         itr.next();
         std::string domain_name = itr.name();

         conduit::Node dest; 
         conduit::Node cdest; 
         conduit::blueprint::mesh::topology::rectilinear::to_structured(
            polyuni[domain_name]["topologies/mesh"], dest, cdest);

         polystruct[domain_name]["coordsets/coords"].reset();
         polystruct[domain_name]["coordsets/coords"] = cdest;
         polystruct[domain_name]["topologies/mesh"].reset();
         polystruct[domain_name]["topologies/mesh"] = dest;
         polystruct[domain_name]["topologies/mesh/coordset"] = "coords";
      }

      TBOX_ASSERT(conduit::blueprint::verify("mesh", polystruct, info));

      bp_utils.writeBlueprintMesh(
         polystruct,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         num_hier_patches,
         "mesh",
         "polystruct",
         "polystruct.root",
         "json");
*/
#endif

      bool test1_passed = comm_tester->verifyCommunicationResults();
      if (do_refine) {

         for (int n = 0; n < ntimes_run; ++n) {

            /*
             * Create communication schedules for data refine tests.
             */
            refine_create_time->start();
            for (int i = 0; i < nlevels; ++i) {
               comm_tester->resetRefineSchedule(i);
            }
            refine_create_time->stop();

            /*
             * Perform refine data communication operations.
             */
            refine_comm_time->start();
            for (int j = 0; j < nlevels; ++j) {
               comm_tester->performRefineOperations(j);
            }
            refine_comm_time->stop();

         }

      }

      if (do_coarsen) {

         for (int n = 0; n < ntimes_run; ++n) {

            /*
             * Create communication schedules for data coarsen tests.
             */
            coarsen_create_time->start();
            for (int i = nlevels - 1; i > 0; --i) {
               comm_tester->resetCoarsenSchedule(i);
            }
            coarsen_create_time->stop();

            /*
             * Perform coarsen data communication operations.
             */
            coarsen_comm_time->start();
            for (int j = nlevels - 1; j > 0; --j) {
               comm_tester->performCoarsenOperations(j);
            }
            coarsen_comm_time->stop();

         }

      }

      bool test2_passed = comm_tester->verifyCommunicationResults();
      /*
       * Deallocate objects when done.
       */

      if (patch_data_test) delete patch_data_test;

      tbox::TimerManager::getManager()->print(tbox::plog);

      tbox::plog << "\nInput file data at end of run is ...." << std::endl;
      input_db->printClassData(tbox::plog);

      if (test1_passed && test2_passed && composite_test_passed) {
         tbox::pout << "\nPASSED:  communication" << std::endl;
         return_val = 0;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   // 0 if passed, 1 otherwise
   return return_val;
}

void transformToPolyhedral(conduit::Node& out, const conduit::Node& in, const std::string& name)
{
   conduit::NodeConstIterator itr = in.children();

   std::map<int, std::map<int, std::vector<int64_t> > > poly_elems_map;

   while(itr.has_next())
   {
      const conduit::Node& child = itr.next();
      std::string domain_name = itr.name();
      out[domain_name]["state"] = child["state"];
      const conduit::Node& in_coords = child["coordsets/coords"];

      int64_t domain_id = child["state/domain_id"].as_int64();
      std::string win_name =
         "window_" + tbox::Utilities::intToString(domain_id, 6);

      const conduit::Node& in_topo = child["topologies"][name];

      int64_t niwidth = in_topo["elements/dims/i"].as_int64() + 1;

      int64_t i_lo = in_topo["elements/origin/i0"].as_int64();
      int64_t j_lo = in_topo["elements/origin/j0"].as_int64();

      const conduit::Node* in_parent = child.parent();

      if (child.has_path("adjsets/adjset/groups")) {
         const conduit::Node& in_groups = child["adjsets/adjset/groups"];
         conduit::NodeConstIterator grp_itr = in_groups.children();
         while(grp_itr.has_next())
         {
            const conduit::Node& group = grp_itr.next();
            std::string grp_name = grp_itr.name();

            if (group.has_child("neighbors")) {
               conduit::int64_array neighbors = group["neighbors"].as_int64_array();

               int nbr_id = neighbors[1];

               if (group.has_child("windows")) {
                  const conduit::Node& in_windows = group["windows"];
                  std::string nbr_win_name =
                     "window_" + tbox::Utilities::intToString(nbr_id, 6);

                  const conduit::Node& ref_win = in_windows[win_name];
                  const conduit::Node& nbr_win = in_windows[nbr_win_name];
                  if (nbr_win["level_id"].as_int64() < ref_win["level_id"].as_int64()) {

                     int64_t ref_size_i = ref_win["dims/i"].as_int64();
                     int64_t ref_size_j = ref_win["dims/j"].as_int64();
                     int64_t ref_size = ref_size_i * ref_size_j;

                     int64_t nbr_size_i = nbr_win["dims/i"].as_int64();
                     int64_t nbr_size_j = nbr_win["dims/j"].as_int64();
                     int64_t nbr_size = nbr_size_i * nbr_size_j;

                     std::string nbr_name =
                        "domain_" + tbox::Utilities::intToString(nbr_id, 6);

                     if (nbr_size < ref_size && !in_parent->has_child(nbr_name)) {
                        // do sends
                        std::vector<double> xbuffer;
                        std::vector<double> ybuffer;
                        const conduit::Node& fcoords =
                           in_coords["values"];
                        const conduit::double_array& xarray =
                           fcoords["x"].as_double_array();
                        const conduit::double_array& yarray =
                           fcoords["y"].as_double_array();

                        int64_t origin_i = ref_win["origin/i"].as_int64();
                        int64_t origin_j = ref_win["origin/j"].as_int64();

                        if (ref_size_i == 1) {
                           int icnst = origin_i - i_lo;
                           int jstart = origin_j - j_lo;
                           int jend = jstart + ref_size_j;
                           for (int jidx = jstart; jidx < jend; ++jidx) {
                              int offset = jidx * niwidth + icnst;
                              xbuffer.push_back(xarray[offset]);
                              ybuffer.push_back(yarray[offset]);
                           }
                        } else if (ref_size_j == 1) {
                           int jcnst = origin_j - j_lo;
                           int istart = origin_i - i_lo;
                           int iend = istart + ref_size_i;
                           for (int iidx = istart; iidx < iend; ++iidx) {
                              int offset = jcnst * niwidth + iidx;
                              xbuffer.push_back(xarray[offset]);
                              ybuffer.push_back(yarray[offset]);
                           }
                        }
                        int64_t nbr_rank = group["rank"].as_int64();
                        MPI_Send(&xbuffer[0],
                                 xbuffer.size(),
                                 MPI_DOUBLE,
                                 nbr_rank,
                                 domain_id,
                                 MPI_COMM_WORLD);
                        MPI_Send(&ybuffer[0],
                                 ybuffer.size(),
                                 MPI_DOUBLE,
                                 nbr_rank,
                                 domain_id,
                                 MPI_COMM_WORLD);
                     }
                  }
               }
            }
         }
      }
   }

   std::map<int, std::map<int, std::vector<double> > > dom_to_nbr_to_xbuffer;
   std::map<int, std::map<int, std::vector<double> > > dom_to_nbr_to_ybuffer;

   itr = in.children();
   while(itr.has_next())
   {
      const conduit::Node& child = itr.next();
      std::string domain_name = itr.name();

      int64_t domain_id = child["state/domain_id"].as_int64();
      std::string win_name =
         "window_" + tbox::Utilities::intToString(domain_id, 6);

      const conduit::Node* in_parent = child.parent();

      auto& nbr_to_xbuffer = dom_to_nbr_to_xbuffer[domain_id];
      auto& nbr_to_ybuffer = dom_to_nbr_to_ybuffer[domain_id];

      if (child.has_path("adjsets/adjset/groups")) {
         const conduit::Node& in_groups = child["adjsets/adjset/groups"];
         conduit::NodeConstIterator grp_itr = in_groups.children();
         while(grp_itr.has_next())
         {
            const conduit::Node& group = grp_itr.next();
            std::string grp_name = grp_itr.name();

            if (group.has_child("neighbors")) {
               conduit::int64_array neighbors = group["neighbors"].as_int64_array();

               int nbr_id = neighbors[1];
               if (group.has_child("windows")) {
                  const conduit::Node& in_windows = group["windows"];
                  std::string nbr_win_name =
                     "window_" + tbox::Utilities::intToString(nbr_id, 6);

                  const conduit::Node& ref_win = in_windows[win_name];
                  const conduit::Node& nbr_win = in_windows[nbr_win_name];
                  if (nbr_win["level_id"].as_int64() > ref_win["level_id"].as_int64()) {
                     int64_t ref_size_i = ref_win["dims/i"].as_int64();
                     int64_t ref_size_j = ref_win["dims/j"].as_int64();
                     int64_t ref_size = ref_size_i * ref_size_j;

                     int64_t nbr_size_i = nbr_win["dims/i"].as_int64();
                     int64_t nbr_size_j = nbr_win["dims/j"].as_int64();
                     int64_t nbr_size = nbr_size_i * nbr_size_j;

                     if (nbr_size > ref_size) {

                        std::string nbr_name =
                           "domain_" + tbox::Utilities::intToString(nbr_id, 6);

                        auto& xbuffer = nbr_to_xbuffer[nbr_id];
                        auto& ybuffer = nbr_to_ybuffer[nbr_id];

                        if (!in_parent->has_child(nbr_name)) {

                           if (nbr_size_i == 1) {
                              xbuffer.resize(nbr_size_j);
                              ybuffer.resize(nbr_size_j);
                           } else if (nbr_size_j == 1) {
                              xbuffer.resize(nbr_size_i);
                              ybuffer.resize(nbr_size_i);
                           }

                           int64_t nbr_rank = group["rank"].as_int64();
                           MPI_Recv(&xbuffer[0],
                                    xbuffer.size(),
                                    MPI_DOUBLE,
                                    nbr_rank,
                                    nbr_id,
                                    MPI_COMM_WORLD,
                                    MPI_STATUS_IGNORE);
                           MPI_Recv(&ybuffer[0],
                                    ybuffer.size(),
                                    MPI_DOUBLE, nbr_rank,
                                    nbr_id, MPI_COMM_WORLD,
                                    MPI_STATUS_IGNORE);

                        } else {

                           const conduit::Node& nbr_dom =
                              (*in_parent)[nbr_name];
                           const conduit::Node& nbr_coords =
                              nbr_dom["coordsets/coords"];

                           const conduit::Node& ntopo =
                              nbr_dom["topologies"][name];
                           int64_t ni_lo =
                              ntopo["elements/origin/i0"].as_int64();
                           int64_t nj_lo =
                              ntopo["elements/origin/j0"].as_int64();
                           int64_t nbr_iwidth =
                              ntopo["elements/dims/i"].as_int64() + 1;

                           const conduit::Node& fcoords =
                              nbr_coords["values"];
                           const conduit::double_array& xarray =
                              fcoords["x"].as_double_array();
                           const conduit::double_array& yarray =
                              fcoords["y"].as_double_array();

                           int64_t origin_i = nbr_win["origin/i"].as_int64();
                           int64_t origin_j = nbr_win["origin/j"].as_int64();

                           if (nbr_size_i == 1) {
                              int icnst = origin_i - ni_lo;
                              int jstart = origin_j - nj_lo;
                              int jend = jstart + nbr_size_j;
                              for (int jidx = jstart; jidx < jend; ++jidx) {
                                 int offset = jidx * nbr_iwidth + icnst;
                                 xbuffer.push_back(xarray[offset]);
                                 ybuffer.push_back(yarray[offset]);
                              }
                           } else if (nbr_size_j == 1) {
                              int jcnst = origin_j - nj_lo;
                              int istart = origin_i - ni_lo;
                              int iend = istart + nbr_size_i;
                              for (int iidx = istart; iidx < iend; ++iidx) {
                                 int offset = jcnst * nbr_iwidth + iidx;
                                 xbuffer.push_back(xarray[offset]);
                                 ybuffer.push_back(yarray[offset]);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   itr = in.children();
   while(itr.has_next())
   {
      const conduit::Node& child = itr.next();
      std::string domain_name = itr.name();
      conduit::Node& out_coords = out[domain_name]["coordsets/coords"];
      const conduit::Node& in_coords = child["coordsets/coords"];

      int64_t domain_id = child["state/domain_id"].as_int64();
      std::string win_name =
         "window_" + tbox::Utilities::intToString(domain_id, 6);

      conduit::Node& out_values = out_coords["values"];
      if (in_coords["type"].as_string() == "uniform") {
         conduit::blueprint::mesh::coordset::uniform::to_explicit(in_coords, out_coords);
      } else {
         out_coords["type"] = in_coords["type"];
         const conduit::Node& in_values = in_coords["values"];
         out_values = in_values;
      }

      auto& nbr_to_xbuffer = dom_to_nbr_to_xbuffer[domain_id];
      auto& nbr_to_ybuffer = dom_to_nbr_to_ybuffer[domain_id];

      const conduit::Node& in_topo = child["topologies"][name];

      int64_t iwidth = in_topo["elements/dims/i"].as_int64();
      int64_t jwidth = in_topo["elements/dims/j"].as_int64();

      int64_t i_lo = in_topo["elements/origin/i0"].as_int64();
      int64_t j_lo = in_topo["elements/origin/j0"].as_int64();

      std::map<int, std::vector<int64_t> >& poly_elems = poly_elems_map[domain_id];

      if (child.has_path("adjsets/adjset/groups")) {
         const conduit::Node& in_groups = child["adjsets/adjset/groups"];
         conduit::NodeConstIterator grp_itr = in_groups.children();
         while(grp_itr.has_next())
         {
            const conduit::Node& group = grp_itr.next();
            std::string grp_name = grp_itr.name();

            if (group.has_child("neighbors")) {
               conduit::int64_array neighbors = group["neighbors"].as_int64_array();

               int nbr_id = neighbors[1];
               if (group.has_child("windows")) {
                  const conduit::Node& in_windows = group["windows"];
                  std::string nbr_win_name =
                     "window_" + tbox::Utilities::intToString(nbr_id, 6);

                  const conduit::Node& ref_win = in_windows[win_name]; 
                  const conduit::Node& nbr_win = in_windows[nbr_win_name]; 
                  if (nbr_win["level_id"].as_int64() > ref_win["level_id"].as_int64()) {
                     int64_t ratio_i = nbr_win["ratio/i"].as_int64();
                     int64_t ratio_j = nbr_win["ratio/j"].as_int64();

                     int64_t ref_size_i = ref_win["dims/i"].as_int64();
                     int64_t ref_size_j = ref_win["dims/j"].as_int64();
                     int64_t ref_size = ref_size_i * ref_size_j;

                     int64_t nbr_size_i = nbr_win["dims/i"].as_int64();
                     int64_t nbr_size_j = nbr_win["dims/j"].as_int64();
                     int64_t nbr_size = nbr_size_i * nbr_size_j;

                     if (ref_size < nbr_size) {

                        int64_t origin_iref = ref_win["origin/i"].as_int64();
                        int64_t origin_jref = ref_win["origin/j"].as_int64();

                        const conduit::double_array& out_x = out_values["x"].as_double_array();
                        const conduit::double_array& out_y = out_values["y"].as_double_array();
                        int64_t new_vertex = out_x.number_of_elements();

                        if (ref_size_i == 1) {
                           int jstart = origin_jref - j_lo;
                           int jend = origin_jref - j_lo + ref_size_j - 1;
                           if (origin_iref == i_lo) {
                              for (int jidx = jstart; jidx < jend; ++jidx) {
                                 int offset = jidx * iwidth;
                                 std::vector<int64_t>& elem_conn = poly_elems[offset];
                                 if (elem_conn.empty()) {
                                    makeQuadElemConnectivity(
                                       elem_conn,
                                       offset,
                                       iwidth);
                                 }
                              }
                           } else {
                              for (int jidx = jstart; jidx < jend; ++jidx) {
                                 int offset = jidx * iwidth +
                                    (origin_iref - i_lo - 1);
                                 std::vector<int64_t>& elem_conn = poly_elems[offset];
                                 if (elem_conn.empty()) {
                                    makeQuadElemConnectivity(
                                       elem_conn,
                                       offset,
                                       iwidth);
                                 }
                              }
                           }
                        } else if (ref_size_j == 1) {
                           int istart = origin_iref - i_lo;
                           int iend = origin_iref - i_lo + ref_size_i - 1;
                           if (origin_jref == j_lo) {
                              for (int iidx = istart; iidx < iend; ++iidx) {
                                 std::vector<int64_t>& elem_conn = poly_elems[iidx];
                                 if (elem_conn.empty()) {
                                    makeQuadElemConnectivity(
                                       elem_conn,
                                       iidx,
                                       iwidth);
                                 }
                              }
                           } else {
                              for (int iidx = istart; iidx < iend; ++iidx) {
                                 int offset = iidx +
                                    ((origin_jref - j_lo - 1) * iwidth);
                                 std::vector<int64_t>& elem_conn = poly_elems[offset];
                                 if (elem_conn.empty()) {
                                    makeQuadElemConnectivity(
                                       elem_conn,
                                       offset,
                                       iwidth);
                                 }
                              }
                           }
                        }

                        int use_ratio = 0;
                        if (nbr_size_j == 1) {
                           use_ratio = ratio_i;
                        }
                        else if (nbr_size_i == 1) {
                           use_ratio = ratio_j;
                        }

                        auto& xbuffer = nbr_to_xbuffer[nbr_id];
                        auto& ybuffer = nbr_to_ybuffer[nbr_id];

                        size_t added = 0;
                        if (nbr_size_j == 1) {
                           added = xbuffer.size() - ref_size_i;
                        } else if (nbr_size_i == 1) {
                           added = ybuffer.size() - ref_size_j;
                        }

                        size_t out_x_size = out_x.number_of_elements();
                        size_t out_y_size = out_y.number_of_elements();

                        std::vector<double> new_x;
                        std::vector<double> new_y;
                        new_x.reserve(out_x_size + added); 
                        new_y.reserve(out_y_size + added); 
                        const double* out_x_ptr = static_cast<const double*>(out_x.element_ptr(0));
                        const double* out_y_ptr = static_cast<const double*>(out_y.element_ptr(0));

                        new_x.insert(new_x.end(), out_x_ptr, out_x_ptr + out_x_size);
                        new_y.insert(new_y.end(), out_y_ptr, out_y_ptr + out_y_size);

                        if ((xbuffer.size()-1)%use_ratio) {
                           new_x.reserve(out_x_size + added*2); 
                        }
                        for (unsigned int ni = 0; ni < xbuffer.size(); ++ni) {
                           if (ni % use_ratio) {
                              new_x.push_back(xbuffer[ni]);
                              new_y.push_back(ybuffer[ni]);
                           }
                        }

                        //out_values["x"].reset(); 
                        //out_values["y"].reset(); 
                        out_values["x"].set(new_x);
                        out_values["y"].set(new_y);

                        if (ref_size_i == 1) {
                           int jstart = origin_jref - j_lo;
                           int jend = origin_jref - j_lo + ref_size_j - 1;
                           if (origin_iref == i_lo) {
                              for (int jidx = jstart; jidx < jend; ++jidx) {
                                 int offset = jidx * (iwidth);
                                 auto& elem_conn = poly_elems[offset];
                                 if (use_ratio > 1) {
                                    for (int nr = use_ratio-1; nr > 0; --nr) {
                                       elem_conn.push_back(new_vertex+nr-1);
                                    }
                                    elem_conn[0] += use_ratio - 1;
                                    new_vertex += use_ratio - 1;
                                 }
                              }
                           } else {
                              for (int jidx = jstart; jidx < jend; ++jidx) {
                                 int offset = jidx * iwidth +
                                    (origin_iref - i_lo - 1);
                                 auto& elem_conn = poly_elems[offset];
                                 if (use_ratio > 1) {
                                    size_t new_size = elem_conn.size() +
                                       use_ratio - 1;
                                    elem_conn.resize(new_size);
                                    int corner = 2;
                                    if (elem_conn[2] - elem_conn[1] != 1) {
                                       int64_t ioff = offset % iwidth;
                                       int64_t joff = offset / iwidth;
                                       int64_t target =
                                          (iwidth+1)*joff + ioff + 1;
                                       for (int nr = 2; nr < 2+use_ratio; ++nr) {
                                          if (elem_conn[nr] == target) {
                                             corner = nr;
                                             break;
                                          }
                                       }
                                    }
                                    for (int nr = new_size-1; nr > corner+use_ratio-1; --nr) {
                                       elem_conn[nr] = elem_conn[nr-use_ratio+1];
                                    }
                                    for (int nr = corner+1; nr < corner+use_ratio; ++nr) {
                                       elem_conn[nr] = new_vertex;
                                       ++new_vertex;
                                    }
                                    elem_conn[0] += use_ratio - 1;
                                 }
                              }
                           }
                        } else if (ref_size_j == 1) {
                           int istart = origin_iref - i_lo;
                           int iend = origin_iref - i_lo + ref_size_i - 1;
                           if (origin_jref == j_lo) {
                              for (int iidx = istart; iidx < iend; ++iidx) {
                                 auto& elem_conn = poly_elems[iidx];
                                 if (use_ratio > 1) {
                                    size_t new_size = elem_conn.size() +
                                       use_ratio - 1;
                                    elem_conn.resize(new_size);
                                    for (int nr = new_size-1; nr > 2; --nr) {
                                       elem_conn[nr] = elem_conn[nr-use_ratio+1];
                                    }
                                    for (int nr = 2; nr <= use_ratio; ++nr) {
                                       elem_conn[nr] = (new_vertex+nr-2);  
                                    }
                                    elem_conn[0] += use_ratio - 1;
                                    new_vertex += use_ratio - 1;
                                 }
                              }
                           } else {
                              for (int iidx = istart; iidx < iend; ++iidx) {
                                 int offset = iidx +
                                    ((origin_jref - j_lo - 1) * iwidth);
                                 auto& elem_conn = poly_elems[offset];
                                 if (use_ratio > 1) {
                                    size_t new_size = elem_conn.size() +
                                       use_ratio - 1;
                                    elem_conn.resize(new_size);
                                    int corner = 3;
                                    if (elem_conn[0] != 4) {
                                       int64_t ioff = offset % iwidth;
                                       int64_t joff = offset / iwidth;
                                       int64_t target =
                                          (iwidth+1)*(joff+1) + ioff + 1;
                                       for (int nr = 3; nr < 3+use_ratio; ++nr) {
                                          if (elem_conn[nr] == target) {
                                             corner = nr;
                                             break;
                                          }
                                       }
                                    }
                                    for (int nr = new_size-1; nr > corner+use_ratio-1; --nr) {
                                       elem_conn[nr] = elem_conn[nr-use_ratio+1];
                                    }
                                    for (int nr = corner+use_ratio-1; nr > corner; --nr) {
                                       elem_conn[nr] = new_vertex;
                                       ++new_vertex;
                                    }
                                    elem_conn[0] += use_ratio - 1;
                                 }
                              }
                           }
                        }
                     } 
                  }
               }
            }
         }
      }

      std::string coords =
         child["topologies"][name]["coordset"].as_string(); 
      out[domain_name]["topologies"][name]["coordset"] = coords;

      conduit::Node& topo = out[domain_name]["topologies"][name];

      topo["type"] = "unstructured";
      topo["elements/shape"] = "polygonal";

      int64_t elemsize = iwidth*jwidth;
 
      std::vector<int64_t> connect;
      for (int elem = 0; elem < elemsize; ++elem) {
         auto elem_itr = poly_elems.find(elem);
         if (elem_itr == poly_elems.end()) {
            makeQuadElemConnectivity(connect, elem, iwidth);
         } else {
            std::vector<int64_t>& poly_elem = elem_itr->second;
            connect.insert(connect.end(), poly_elem.begin(), poly_elem.end());
         }
      }

      topo["elements/connectivity"].set(connect);

   }
}

void makeQuadElemConnectivity(
   std::vector<int64_t>& connect,
   int64_t element,
   int64_t iwidth)
{
   int64_t i = element % iwidth;
   int64_t j = element / iwidth;
   int64_t LLi = i;
   int64_t LRi = i + 1;
   int64_t ULi = i;
   int64_t URi = i + 1;
   int64_t LLj = j;
   int64_t LRj = j;
   int64_t ULj = j + 1;
   int64_t URj = j + 1;

   int64_t LL = (iwidth+1)*LLj + LLi;
   int64_t LR = (iwidth+1)*LRj + LRi;
   int64_t UL = (iwidth+1)*ULj + ULi;
   int64_t UR = (iwidth+1)*URj + URi;

   connect.push_back(4);
   connect.push_back(LL);
   connect.push_back(LR);
   connect.push_back(UR);
   connect.push_back(UL);
}

