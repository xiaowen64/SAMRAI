/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Main program to test index data operations 
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Application classes
#include "SAMRAI/appu/EmbeddedBoundaryGeometry.h"
#include "SampleApp.h"

// SAMRAI classes
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisItDataWriter.h"

//#undef RECORD_STATS
#define RECORD_STATS

#ifdef RECORD_STATS
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/Statistic.h"
#endif

using namespace SAMRAI;
using namespace tbox;
using namespace geom;
using namespace appu;
using namespace mesh;

int main(
   int argc,
   char* argv[]) {

   int num_fails = 0;

   /*
    * Initialize MPI and SAMRAI, enable logging, and process command line.
    */

   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      const tbox::Dimension dim(NDIM);

      PIO::logAllNodes("ebgeom.log");

      string input_filename;
      string restart_read_dirname;
      int restore_num = 0;

      bool is_from_restart = false;

      if ((argc != 2) && (argc != 4)) {
         pout << "USAGE:  " << argv[0] << " <input filename> "
                                       <<
         "<restart dir> <restore number> [options]\n"
                                       << "  options:\n"
                                       << "  none at this time"
                                       << endl;
         SAMRAI_MPI::abort();
         return -1;
      } else {
         input_filename = argv[1];
         if (argc == 4) {
            restart_read_dirname = argv[2];
            restore_num = atoi(argv[3]);

            is_from_restart = true;
         }
      }

      plog << "input_filename = " << input_filename << endl;
      plog << "restart_read_dirname = " << restart_read_dirname << endl;
      plog << "restore_num = " << restore_num << endl;

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

      Pointer<Database> input_db = new InputDatabase("input_db");
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section of the input database.
       */

      Pointer<Database> main_db = input_db->getDatabase("Main");

      double time = 0.;
      if (main_db->keyExists("time")) {
         time = main_db->getDouble("time");
      }
      double dt = 0.;
      if (main_db->keyExists("dt")) {
         dt = main_db->getDouble("dt");
      }
      int nsteps = 0;
      if (main_db->keyExists("nsteps")) {
         nsteps = main_db->getInteger("nsteps");
      }

      int tag_buffer = 0;
      if (main_db->keyExists("tag_buffer")) {
         tag_buffer = main_db->getInteger("tag_buffer");
      }

      bool print_boundarynode_data = false;
      if (main_db->keyExists("print_boundarynode_data")) {
         print_boundarynode_data = main_db->getBool("print_boundarynode_data");
      }

      int visit_dump_interval = 0;
      if (main_db->keyExists("visit_dump_interval")) {
         visit_dump_interval = main_db->getInteger("visit_dump_interval");
      }
      string visit_dump_dirname;
      int visit_number_procs_per_file = 1;
      if (visit_dump_interval > 0) {
         if (main_db->keyExists("visit_dump_dirname")) {
            visit_dump_dirname = main_db->getString("visit_dump_dirname");
         }
         if (main_db->keyExists("visit_number_procs_per_file")) {
            visit_number_procs_per_file =
               main_db->getInteger("visit_number_procs_per_file");
         }
      }

#if (TESTING == 1)
      bool check_result = false;
      double correct_volume = 0.;
      if (main_db->keyExists("AutoTesting")) {
         check_result = true;
         Pointer<Database> test_db = main_db->getDatabase("AutoTesting");
         correct_volume = test_db->getDouble("correct_volume");
      }
#endif

      TimerManager::createManager(input_db->getDatabase("TimerManager"));

#ifdef RECORD_STATS
      tbox::Statistician* statistician = tbox::Statistician::createStatistician();
#endif

      /*
       * The grid geometry defines the grid type (e.g. cartesian, spherical,
       * etc.).  Because SAMRAI operates on block structured indices, it can
       * support any grid geometry that may be represented as an orthogonal
       * grid.
       */
      Pointer<CartesianGridGeometry> grid_geometry =
         new CartesianGridGeometry("CartesianGeometry",
            input_db->getDatabase("CartesianGeometry"));

      /*
       * The patch hierarchy defines the adaptive grid system.
       */
      Pointer<PatchHierarchy> patch_hierarchy =
         new PatchHierarchy("PatchHierarchy", grid_geometry);

      /*
       * Set up Visualization writer.
       */
      Pointer<VisItDataWriter> visit_data_writer =
         new VisItDataWriter("EBGeom VisIt Writer",
            visit_dump_dirname,
            visit_number_procs_per_file);

      IntVector<NDIM> ghosts(1);
      Pointer<EmbeddedBoundaryGeometry> eb_geom =
         new EmbeddedBoundaryGeometry("EmbeddedBoundaryGeometry",
            input_db->getDatabase("EmbeddedBoundaryGeometry"),
            grid_geometry,
            ghosts);
      /*
       * Dummy problem class, intended to mimic a real application
       * that would use an EmbeddedBoundaryGeometry.
       */
      SampleApp* app_model = new SampleApp("SampleApp",
            input_db->getDatabase("SampleApp"),
            grid_geometry,
            eb_geom,
            visit_data_writer);

      /*
       * The StandardTagAndInitialize class performs a variety of operations
       * with user-specified parameters related to adptive gridding.  For example,
       * it manages initialization of a level, cell tagging using a gradient
       * detector, and methods to reset data after the hierarchy has been
       * regridded.
       */
      Pointer<StandardTagAndInitialize> tag_and_init_ops =
         new StandardTagAndInitialize(
            "StandardTagAndInitialize",
            app_model,
            input_db->getDatabase("StandardTagAndInitialize"));

      /*
       * The gridding algorithm manages adaptive gridding.  It expects a
       * clustering scheme (i.e. how to cluster tagged-cells into patches),
       * and a load balance scheme to distribute work to processors.  In general
       * the baseline classes provided in SAMRAI should suffice for most
       * problems. It also requires a class that defines the particular tag
       * and initialization ops that correlate with the users problem.  For
       * this, we use the "tag_and_init_ops" above, which references our
       * "dummy_app" problem class to define the user-specific operations.
       */
      Pointer<BergerRigoutsos> box_generator = new BergerRigoutsos();

      Pointer<TreeLoadBalancer> load_balancer =
         new TreeLoadBalancer(dim, "LoadBalancer",
            input_db->getDatabase("LoadBalancer"));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getCommunicator());

#ifdef USE_NONUNIFORM_LB
      // set the workload factor in the load balancer
      int weight_id = app_model->getWorkloadIndex();
      load_balancer->setWorkloadPatchDataIndex(weight_id);
#endif

      Pointer<GriddingAlgorithm> gridding_algorithm =
         new GriddingAlgorithm(
            dim,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            tag_and_init_ops,
            box_generator,
            load_balancer);

/*
 ************************************************************************
 *
 * Build hierarchy, and build embedded boundary on each level.
 *
 ************************************************************************
 */

      /*
       * Create timers.
       */
      Pointer<Timer> t_create_eb = TimerManager::getManager()->
         getTimer("apps::main::create_eb");
      Pointer<Timer> t_write_viz = TimerManager::getManager()->
         getTimer("apps::main::write_viz");
      Pointer<Timer> t_init_lev = TimerManager::getManager()->
         getTimer("apps::main::init_lev");
      Pointer<Timer> t_regrid = TimerManager::getManager()->
         getTimer("apps::main::regrid");

      /*
       * Set parameters in the embedded boundary geometry object.
       */
      eb_geom->registerVisItDataWriter(visit_data_writer);

      /*
       * Build embedded boundary on coarsest level
       */
      gridding_algorithm->makeCoarsestLevel(patch_hierarchy, time);
      Pointer<PatchLevel> coarse_level = patch_hierarchy->getPatchLevel(0);

      t_create_eb->start();
      eb_geom->buildEmbeddedBoundaryOnLevel(coarse_level);
      //app_model->printBoundaryNodeData(coarse_level, eb_geom);
      t_create_eb->stop();

#ifdef USE_NONUNIFORM_LB
      // set weights of cut cells
      for (PatchLevel::Iterator ip(coarse_level); ip; ip++) {
         Pointer<Patch> patch = coarse_level->getPatch(ip());
         app_model->setWeightOnPatch(patch,
            time,
            true);
      }
#endif

      /*
       * Build embedded boundary on finer levels
       */
      bool done = false;
      bool initial_time = true;
      int init_tag_buffer = tag_buffer;
      int ln;
      for (ln = 0; gridding_algorithm->levelCanBeRefined(ln) && !done;
           ln++) {
         gridding_algorithm->makeFinerLevel(patch_hierarchy,
            time,
            initial_time,
            init_tag_buffer);

         if (patch_hierarchy->finerLevelExists(ln)) {
            Pointer<PatchLevel> fine_level = patch_hierarchy->getPatchLevel(
                  ln + 1);
            t_create_eb->start();
            eb_geom->buildEmbeddedBoundaryOnLevel(fine_level, patch_hierarchy);
            t_create_eb->stop();

#ifdef USE_NONUNIFORM_LB
            // set weights of cut cells
            for (PatchLevel::Iterator ip(fine_level); ip; ip++) {
               Pointer<Patch> patch = fine_level->getPatch(ip());
               app_model->setWeightOnPatch(patch,
                  time,
                  initial_time);
            }
#endif

            if (print_boundarynode_data) {
               app_model->printBoundaryNodeData(fine_level, eb_geom);
            }

         }

         done = !(patch_hierarchy->finerLevelExists(ln));
      }

      pout.precision(12);
      for (ln = 0; ln < patch_hierarchy->getNumberOfLevels(); ln++) {
         Pointer<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
         double total_vol = eb_geom->computeTotalVolumeOnLevel(level);

#if (TESTING == 1)
         if (check_result) {
            if (!(tbox::MathUtilities<double>::equalEps(total_vol,
                     correct_volume))) {
               pout << "emb_bdry: Volume Test FAILED"
                    << "\n\tcomputed: " << total_vol
                    << "\n\tcorrect:  " << correct_volume << endl;
            }
            ++num_fails;
         }
#else
         pout << "Level: " << ln << "\tTotal Volume: " << total_vol << endl;
#endif

      }

      if (visit_dump_interval > 0) {
         t_write_viz->start();
         visit_data_writer->writePlotData(patch_hierarchy, 0, time);
         t_write_viz->stop();
      }

      /*****************************************************************
       *
       * Re-generate grids (if desired)
       *
       ****************************************************************/
      int nlevels = patch_hierarchy->getNumberOfLevels();
      Array<int> regrid_tag_buffer(nlevels);
      for (ln = 0; ln < nlevels; ln++) {
         regrid_tag_buffer[ln] = tag_buffer;
      }

      for (int nstep = 1; nstep < nsteps; nstep++) {
         time = time + dt;

         /*
          * Tag data on coarsest level
          */
         Pointer<PatchLevel> coarsest_level =
            patch_hierarchy->getPatchLevel(0);
         t_init_lev->start();
         app_model->initializeLevelData(patch_hierarchy,
            0,
            time,
            true,
            false,
            coarsest_level,
            false);
         t_init_lev->stop();

         /*
          * Regrid
          */
         pout << "\n\n############################################" << endl;
         pout << "              REGRIDDING STEP: " << nstep << endl;
         t_regrid->start();
         gridding_algorithm->regridAllFinerLevels(patch_hierarchy,
            0,
            time,
            regrid_tag_buffer);
         t_regrid->stop();
         pout << "############################################\n\n" << endl;

         /*
          * Reset EB on all levels
          */
         for (ln = 0; ln < nlevels; ln++) {
            Pointer<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            t_create_eb->start();
            eb_geom->buildEmbeddedBoundaryOnLevel(level);
            t_create_eb->stop();

#ifdef USE_NONUNIFORM_LB
            // set weights of cut cells
            for (PatchLevel::Iterator ip(level); ip; ip++) {
               Pointer<Patch> patch = level->getPatch(ip());
               app_model->setWeightOnPatch(patch,
                  time,
                  false);
            }
#endif

         }

         pout.precision(12);
         for (ln = 0; ln < patch_hierarchy->getNumberOfLevels(); ln++) {
            Pointer<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            double total_vol = eb_geom->computeTotalVolumeOnLevel(level);
            pout << "Level: " << ln << "\tTotal Volume: " << total_vol << endl;
         }

         t_write_viz->start();
         visit_data_writer->writePlotData(patch_hierarchy, nstep, time);
         t_write_viz->stop();

      } // loop over nsteps

      /*
       * Output timer results.
       */
#if (TESTING != 1)
      TimerManager::getManager()->print(pout);
#endif

#ifdef RECORD_STATS
      statistician->finalize();
      statistician->printAllSummedGlobalStatData("stats.txt");
//   statistician->finalize();
      statistician->printSpreadSheetOutput("stats");
#endif

      grid_geometry.setNull();
      patch_hierarchy.setNull();
      eb_geom.setNull();

      if (app_model) delete app_model;

#if (TESTING == 1)
      if (num_fails == 0) {
         tbox::pout << "\nPASSED:  Embedded Boundary" << endl;
      }
#endif

   }

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   return num_fails;
}
