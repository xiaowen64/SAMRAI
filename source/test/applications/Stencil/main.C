/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Main program for SAMRAI Linear Advection example problem.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Header for application-specific algorithm/data structure object

#include "Stencil.h"

// Headers for major algorithm/data structure objects

#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/mesh/TileClustering.h"

// Headers for basic SAMRAI objects

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/shared_ptr.hpp"

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>

#include "SAMRAI/tbox/StartupShutdownManager.h"
#include <cuda_runtime.h>

using namespace std;
using namespace SAMRAI;

class DeviceManager
{
private:
   static void initializeCallback()
   {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      if (mpi.getSize() < 2) return;

      MPI_Comm node_comm;
      int local_size, local_rank, rank;

      MPI_Comm_split(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &node_comm);
      MPI_Comm_size(node_comm, &local_size);
      MPI_Comm_rank(node_comm, &local_rank);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      // Assign device to MPI process
      int num_devices;
      cudaGetDeviceCount(&num_devices);

      const int my_device = local_rank % num_devices;
      cudaSetDevice(my_device);
      std::cout << "My device = " << my_device << std::endl;
   }

   static void finalizeCallback()
   {
      cudaDeviceReset();
   }

   static SAMRAI::tbox::StartupShutdownManager::Handler s_initialize_handler;
};

SAMRAI::tbox::StartupShutdownManager::Handler
DeviceManager::s_initialize_handler(
   DeviceManager::initializeCallback,
   0,
   0,
   DeviceManager::finalizeCallback,
   tbox::StartupShutdownManager::priorityLogger);

int main(
    int argc,
    char* argv[])
{

  /*
   * Initialize tbox::MPI.
   */

  tbox::SAMRAI_MPI::init(&argc, &argv);
  tbox::SAMRAIManager::initialize();

  /*
   * Initialize SAMRAI, enable logging, and process command line.
   */
  tbox::SAMRAIManager::startup();
  const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

  {
    string input_filename;
    string restart_read_dirname;
    int restore_num = 0;

    bool is_from_restart = false;

    if ((argc != 2) && (argc != 4)) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
        << "<restart dir> <restore number> [options]\n"
        << "  options:\n"
        << "  none at this time"
        << endl;
      tbox::SAMRAI_MPI::abort();
      return -1;
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

    boost::shared_ptr<tbox::InputDatabase> input_db(
        new tbox::InputDatabase("input_db"));
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    /*
     * Retrieve "Main" section of the input database.  First, read
     * dump information, which is used for writing plot files.
     * Second, if proper restart information was given on command
     * line, and the restart interval is non-zero, create a restart
     * database.
     */

    boost::shared_ptr<tbox::Database> main_db(
        input_db->getDatabase("Main"));

    const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

    const std::string base_name =
      main_db->getStringWithDefault("base_name", "unnamed");

    const std::string log_filename =
      main_db->getStringWithDefault("log_filename", base_name + ".log");

    bool log_all_nodes = false;
    if (main_db->keyExists("log_all_nodes")) {
      log_all_nodes = main_db->getBool("log_all_nodes");
    }
    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_filename);
    } else {
      tbox::PIO::logOnlyNodeZero(log_filename);
    }

    int viz_dump_interval = 1;
    if (main_db->keyExists("viz_dump_interval")) {
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
    }

    const std::string viz_dump_dirname =
      main_db->getStringWithDefault("viz_dump_dirname", base_name + ".visit");
    int visit_number_procs_per_file = 1;

    const bool viz_dump_data = (viz_dump_interval > 0);

    int restart_interval = 0;
    if (main_db->keyExists("restart_interval")) {
      restart_interval = main_db->getInteger("restart_interval");
    }

    const std::string restart_write_dirname =
      main_db->getStringWithDefault("restart_write_dirname",
          base_name + ".restart");

    bool use_refined_timestepping = false;

    if (main_db->keyExists("timestepping")) {
      string timestepping_method = main_db->getString("timestepping");
      if (timestepping_method == "REFINED") {
        use_refined_timestepping = true;
      }
    }

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
            mpi.getSize());
    }

    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

    /*
     * Create major algorithm and data objects which comprise application.
     * Each object will be initialized either from input data or restart
     * files, or a combination of both.  Refer to each class constructor
     * for details.  For more information on the composition of objects
     * for this application, see comments at top of file.
     */

    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
        new geom::CartesianGridGeometry(
          dim,
          "CartesianGeometry",
          input_db->getDatabase("CartesianGeometry")));

    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
        new hier::PatchHierarchy(
          "PatchHierarchy",
          grid_geometry,
          input_db->getDatabase("PatchHierarchy")));

    Stencil* stencil_model = new Stencil(
        "Stencil",
        dim,
        input_db->getDatabase("Stencil"),
        grid_geometry);

    boost::shared_ptr<algs::HyperbolicLevelIntegrator> hyp_level_integrator(
        new algs::HyperbolicLevelIntegrator(
          "HyperbolicLevelIntegrator",
          input_db->getDatabase("HyperbolicLevelIntegrator"),
          stencil_model,
          use_refined_timestepping));

    boost::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
        new mesh::StandardTagAndInitialize(
          "StandardTagAndInitialize",
          hyp_level_integrator.get(),
          input_db->getDatabase("StandardTagAndInitialize")));

    bool use_tile_clustering = false;
    use_tile_clustering = main_db->getBoolWithDefault(
        "use_tile_clustering",
        use_tile_clustering);

    boost::shared_ptr<mesh::BoxGeneratorStrategy> box_generator;

    if (use_tile_clustering) {
      box_generator = boost::make_shared<mesh::TileClustering>(
          dim,
          input_db->getDatabaseWithDefault(
            "TileClustering",
            boost::shared_ptr<tbox::Database>()));
    } else {
      box_generator = boost::make_shared<mesh::BergerRigoutsos>(
            dim,
            input_db->getDatabaseWithDefault(
              "BergerRigoutsos",
              boost::shared_ptr<tbox::Database>()));
      //box_generator->useDuplicateMPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
    }

    boost::shared_ptr<mesh::CascadePartitioner> load_balancer(
        new mesh::CascadePartitioner(
          dim,
          "LoadBalancer",
          input_db->getDatabase("LoadBalancer")));
    load_balancer->setSAMRAI_MPI(
        tbox::SAMRAI_MPI::getSAMRAIWorld());

    boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
        new mesh::GriddingAlgorithm(
          patch_hierarchy,
          "GriddingAlgorithm",
          input_db->getDatabase("GriddingAlgorithm"),
          error_detector,
          box_generator,
          load_balancer));

    boost::shared_ptr<algs::TimeRefinementIntegrator> time_integrator(
        new algs::TimeRefinementIntegrator(
          "TimeRefinementIntegrator",
          input_db->getDatabase("TimeRefinementIntegrator"),
          patch_hierarchy,
          hyp_level_integrator,
          gridding_algorithm));

    boost::shared_ptr<appu::VisItDataWriter> visit_data_writer(
       new appu::VisItDataWriter(
          dim,
          "Euler VisIt Writer",
          viz_dump_dirname,
          visit_number_procs_per_file));
    stencil_model->registerVisItDataWriter(visit_data_writer);

    /*
     * Initialize hierarchy configuration and data on all patches.
     * Then, close restart file and write initial state for visualization.
     */

    double dt_now = time_integrator->initializeHierarchy();

    tbox::RestartManager::getManager()->closeRestartFile();

    /*
     * After creating all objects and initializing their state, we
     * print the input database and variable database contents
     * to the log file.
     */

    tbox::plog << "\nCheck input data and variables before simulation:"
      << endl;
    tbox::plog << "Input database..." << endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

    if ((viz_dump_interval > 0))
      visit_data_writer->writePlotData(
          patch_hierarchy,
          time_integrator->getIntegratorStep(),
          time_integrator->getIntegratorTime());

    // tbox::plog << "\nCheck Linear Advection data... " << endl;
    // linear_advection_model->printClassData(tbox::plog);

    /*
     * Time step loop.  Note that the step count and integration
     * time are maintained by algs::TimeRefinementIntegrator.
     */

    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();

    int iteration_num = time_integrator->getIntegratorStep();

    while ((loop_time < loop_time_end) &&
        time_integrator->stepsRemaining()) {

      iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1
        << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;

      double dt_new = time_integrator->advanceHierarchy(dt_now);

      loop_time += dt_now;
      dt_now = dt_new;

      tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;

       if ((viz_dump_interval > 0)
           && (iteration_num % viz_dump_interval) == 0) {
          visit_data_writer->writePlotData(patch_hierarchy,
             iteration_num,
             loop_time);
       }
    }


    /*
     * Output timer results.
     */
    tbox::TimerManager::getManager()->print(tbox::pout);

    /*
     * At conclusion of simulation, deallocate objects.
     */

    time_integrator.reset();
    gridding_algorithm.reset();
    load_balancer.reset();
    box_generator.reset();
    error_detector.reset();
    hyp_level_integrator.reset();

    if (stencil_model) delete stencil_model;

    patch_hierarchy.reset();
    grid_geometry.reset();

    input_db.reset();
    main_db.reset();

  }

  tbox::SAMRAIManager::shutdown();

  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();

}
