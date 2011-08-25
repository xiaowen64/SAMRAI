/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Test program for performance of tree search algorithm.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/TransferOperatorRegistry.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <algorithm>
#include <vector>
#include <iomanip>

using namespace SAMRAI;
using namespace tbox;


/*
 * Break up boxes in the given MappedBoxLevel.  This method is meant
 * to create a bunch of small boxes from the user-input boxes in order
 * to set up a non-trivial mesh configuration.
 *
 * The database may have the following items: min_box_size,
 * max_box_size, refinement_ratio.
 *
 * The mapped_box_level will be refined by refinement_ratio and
 * partitioned to generate a non-trivial configuration for testing.
 */
void breakUpBoxes(
   hier::MappedBoxLevel &mapped_box_level,
   const hier::MappedBoxLevel &domain_mapped_box_level,
   const tbox::Pointer<tbox::Database> &database);

void alterAndGenerateMapping(
   hier::MappedBoxLevel &mapped_box_level_c,
   hier::Connector &b_to_c,
   hier::Connector &c_to_b,
   const hier::MappedBoxLevel &mapped_box_level_b,
   const tbox::Pointer<tbox::Database> &database );

/*
************************************************************************
*
* This is an accuracy test for the MappingConnectorAlgorithm class:
*
* 1. Read in user-specified GridGeometry.
*
* 2. Build a domain MappedBoxLevel from GridGeometry.
*
* 3. Build MappedBoxLevel A by refining and partitioning domain.
*
* 4. Build MappedBoxLevel B by refining and partitioning domain.
*    Compute overlap Connectors A<==>B.
*
* 5. Build MappedBoxLevel C by changing B based on some simple formula.
*    Generate mapping Connectors B<==>C.
*
* 6. Apply mapping B<==>C to update A<==>B.
*
* 7. Check correctness of updated A<==>B.
*
*************************************************************************
*/


int main(
   int argc,
   char* argv[])
{
   /*
    * Initialize MPI, SAMRAI.
    */

   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   const int rank = mpi.getRank();
   size_t fail_count = 0;

   {

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       * executable <input file name>
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

      Pointer<Database> input_db(new InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      Pointer<Database> main_db = input_db->getDatabase("Main");

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Start logging.
       */
      const std::string log_filename = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_filename);
      } else {
         PIO::logOnlyNodeZero(log_filename);
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);


      /*
       * Generate the GridGeometry.
       */
      tbox::ConstPointer<hier::GridGeometry> grid_geometry;
      if (main_db->keyExists("GridGeometry")) {
         grid_geometry = new hier::GridGeometry(
            dim,
            "GridGeometry",
            tbox::Pointer<hier::TransferOperatorRegistry>(),
            main_db->getDatabase("GridGeometry"));
      } else {
         TBOX_ERROR("Multiblock tree search test: could not find entry GridGeometry"
                    << "\nin input.");
      }



      /*
       * Print input database again to fully show usage.
       */
      plog << "Input database after running..." << std::endl;
      input_db->printClassData(plog);

      const hier::IntVector &one_vector(hier::IntVector::getOne(dim));
      const hier::IntVector &zero_vector(hier::IntVector::getZero(dim));

      hier::BoxSet multiblock_boxes;
      grid_geometry->computePhysicalDomain(
         multiblock_boxes,
         hier::IntVector::getOne(dim) );

      hier::MappedBoxLevel domain_mapped_box_level(
         multiblock_boxes,
         one_vector,
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         hier::MappedBoxLevel::GLOBALIZED);


      /*
       * Generate MappedBoxLevel A from the multiblock domain description
       * using input database MappedBoxLevelA.
       */
      hier::MappedBoxLevel mapped_box_level_a(domain_mapped_box_level);
      Pointer<Database> a_db = main_db->getDatabase("MappedBoxLevelA");
      breakUpBoxes(mapped_box_level_a, domain_mapped_box_level, a_db);
      mapped_box_level_a.cacheGlobalReducedData();
      // tbox::pout << "mapped box level a:\n" << mapped_box_level_a.format("A: ",2) << std::endl;


      /*
       * Generate MappedBoxLevel B from the multiblock domain description
       * using input database MappedBoxLevelB.
       */
      hier::MappedBoxLevel mapped_box_level_b(domain_mapped_box_level);
      Pointer<Database> b_db = main_db->getDatabase("MappedBoxLevelB");
      breakUpBoxes(mapped_box_level_b, domain_mapped_box_level, b_db);
      mapped_box_level_b.cacheGlobalReducedData();
      // tbox::pout << "mapped box level b:\n" << mapped_box_level_b.format("B: ",2) << std::endl;


      /*
       * Generate Connector A<==>B, to be modified by the mapping
       * operation.
       */

      hier::IntVector base_width_a(zero_vector);
      hier::IntVector base_width_b(zero_vector);
      if ( main_db->isInteger("base_width_a") ) {
         main_db->getIntegerArray("base_width_a", &base_width_a[0], dim.getValue());
      }
      base_width_b = hier::Connector::convertHeadWidthToBase(
         mapped_box_level_b.getRefinementRatio(),
         mapped_box_level_a.getRefinementRatio(),
         base_width_a );

      hier::Connector a_to_b( mapped_box_level_a,
                              mapped_box_level_b,
                              base_width_a );
      hier::Connector b_to_a( mapped_box_level_b,
                              mapped_box_level_a,
                              base_width_b );

      hier::OverlapConnectorAlgorithm oca;
      oca.findOverlaps(a_to_b);
      oca.findOverlaps(b_to_a);
      // tbox::pout << "a_to_b:\n" << a_to_b.format("AB: ",2) << std::endl;
      // tbox::pout << "b_to_a:\n" << b_to_a.format("BA: ",2) << std::endl;

      a_to_b.checkConsistencyWithBase();
      a_to_b.checkConsistencyWithHead();
      b_to_a.checkConsistencyWithBase();
      b_to_a.checkConsistencyWithHead();

      oca.checkOverlapCorrectness(b_to_a);
      oca.checkOverlapCorrectness(b_to_a);


      /*
       * Generate MappedBoxLevel C by altering B based on a simple formula.
       * Generate the mapping Connectors B<==>C.
       */

      hier::MappedBoxLevel mapped_box_level_c(dim);
      hier::Connector b_to_c, c_to_b;
      Pointer<Database> alteration_db = main_db->getDatabase("Alteration");

      alterAndGenerateMapping(
         mapped_box_level_c,
         b_to_c,
         c_to_b,
         mapped_box_level_b,
         alteration_db );
      mapped_box_level_c.cacheGlobalReducedData();
      // tbox::pout << "mapped box level c:\n" << mapped_box_level_c.format("C: ",2) << std::endl;
      // tbox::pout << "b_to_c:\n" << b_to_c.format("BC: ",2) << std::endl;
      // tbox::pout << "c_to_b:\n" << c_to_b.format("CB: ",2) << std::endl;

      hier::MappingConnectorAlgorithm mca;
      mca.modify( a_to_b,
                  b_to_a,
                  b_to_c,
                  c_to_b,
                  &mapped_box_level_b,
                  &mapped_box_level_c );
      // tbox::pout << "mapped box level b after modify:\n" << mapped_box_level_b.format("B: ",2) << std::endl;

      // tbox::pout << "checking a--->b consistency with base:" << std::endl;
      a_to_b.checkConsistencyWithBase();
      // tbox::pout << "checking a--->b consistency with head:" << std::endl;
      a_to_b.checkConsistencyWithHead();
      // tbox::pout << "checking b--->a consistency with base:" << std::endl;
      b_to_a.checkConsistencyWithBase();
      // tbox::pout << "checking b--->a consistency with head:" << std::endl;
      b_to_a.checkConsistencyWithHead();

      tbox::pout << "Checking for a--->b errors:" << std::endl;
      const size_t a_to_b_errors = oca.checkOverlapCorrectness(a_to_b);
      if ( a_to_b_errors ) {
         tbox::pout << "... " << a_to_b_errors << " errors." << std::endl;
      }
      else {
         tbox::pout << "... none." << std::endl;
      }

      tbox::pout << "Checking for b--->a errors:" << std::endl;
      const size_t b_to_a_errors = oca.checkOverlapCorrectness(b_to_a);
      if ( b_to_a_errors ) {
         tbox::pout << "... " << b_to_a_errors << " errors." << std::endl;
      }
      else {
         tbox::pout << "... none." << std::endl;
      }

      fail_count += a_to_b_errors + b_to_a_errors;



      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  Connector modify" << std::endl;
      }

      input_db.setNull();
      main_db.setNull();

      /*
       * Exit properly by shutting down services in correct order.
       */
      tbox::plog << "\nShutting down..." << std::endl;

   }

   /*
    * Shut down.
    */
   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();

   if (fail_count == 0) {
      SAMRAI_MPI::finalize();
   } else {
      tbox::pout << "Process " << std::setw(5) << rank << " aborting."
                 << std::endl;
      SAMRAI::tbox::Utilities::abort("Aborting due to nonzero fail count",
         __FILE__, __LINE__);
   }

   tbox::plog << "Process " << std::setw(5) << rank << " exiting." << std::endl;
   return int(fail_count);
}





/*
 * Break up boxes in the given MappedBoxLevel.  This method is meant
 * to create a bunch of small boxes from the user-input boxes in order
 * to set up a non-trivial mesh configuration.
 *
 * 1. Refine the boxes according to refinement_ratio in the database.
 * 2. Partition according to min and max box sizes in the database.
 */
void breakUpBoxes(
   hier::MappedBoxLevel &mapped_box_level,
   const hier::MappedBoxLevel &domain_mapped_box_level,
   const tbox::Pointer<tbox::Database> &database) {

   const tbox::Dimension &dim(mapped_box_level.getDim());

   hier::IntVector refinement_ratio(hier::IntVector::getOne(dim));
   if ( database->isInteger("refinement_ratio") ) {
      database->getIntegerArray("refinement_ratio", &refinement_ratio[0], dim.getValue());
   }

   if ( refinement_ratio != hier::IntVector::getOne(dim) ) {
      const hier::BoxSet &mapped_boxes(mapped_box_level.getMappedBoxes());
      hier::BoxSet refined_mapped_boxes;
      for ( hier::BoxSet::const_iterator bi=mapped_boxes.begin();
            bi!=mapped_boxes.end(); ++bi ) {
         hier::Box refined_mapped_box(*bi);
         refined_mapped_box.refine(refinement_ratio);
         refined_mapped_boxes.insert( refined_mapped_boxes.end(), refined_mapped_box );
      }
      mapped_box_level.swapInitialize(
         refined_mapped_boxes,
         mapped_box_level.getRefinementRatio()*refinement_ratio,
         mapped_box_level.getGridGeometry(),
         mapped_box_level.getMPI() );
   }

   hier::IntVector max_box_size(dim,tbox::MathUtilities<int>::getMax());
   if ( database->isInteger("max_box_size") ) {
      database->getIntegerArray("max_box_size", &max_box_size[0], dim.getValue());
   }

   hier::IntVector min_box_size(hier::IntVector::getOne(dim));
   if ( database->isInteger("min_box_size") ) {
      database->getIntegerArray("min_box_size", &min_box_size[0], dim.getValue());
   }

   mesh::TreeLoadBalancer load_balancer( mapped_box_level.getDim() );

   const tbox::Pointer<hier::PatchHierarchy> hierarchy;
   const int level_number(0);

   hier::Connector dummy_connector;

   const hier::IntVector bad_interval(dim,1);
   const hier::IntVector cut_factor(dim,1);

   load_balancer.loadBalanceMappedBoxLevel(
      mapped_box_level,
      dummy_connector,
      dummy_connector,
      tbox::Pointer<hier::PatchHierarchy>(),
      level_number,
      dummy_connector,
      dummy_connector,
      min_box_size,
      max_box_size,
      domain_mapped_box_level,
      bad_interval,
      cut_factor);

   return;
}


/*
 * Generate MappedBoxLevel C by altering B based on a simple formula.
 * Generate the mapping Connectors B<==>C.
 */
void alterAndGenerateMapping(
   hier::MappedBoxLevel &mapped_box_level_c,
   hier::Connector &b_to_c,
   hier::Connector &c_to_b,
   const hier::MappedBoxLevel &mapped_box_level_b,
   const tbox::Pointer<tbox::Database> &database )
{
   const tbox::Dimension dim(mapped_box_level_b.getDim());

   /*
    * Increment for changing the LocalIds.
    * Set to zero to disable.
    */
   const int local_id_increment =
      database->getIntegerWithDefault("local_id_increment", 0);

   const hier::BoxSet mapped_boxes_b(mapped_box_level_b.getMappedBoxes());
   hier::BoxSet mapped_boxes_c;
   hier::NeighborhoodSet b_eto_c;
   hier::NeighborhoodSet c_eto_b;
   for ( hier::BoxSet::const_iterator bi(mapped_boxes_b.begin());
         bi!=mapped_boxes_b.end(); ++bi ) {
      const hier::Box &mapped_box_b(*bi);
      hier::Box mapped_box_c(mapped_box_b,
                                   mapped_box_b.getLocalId() + local_id_increment,
                                   mapped_box_b.getOwnerRank(),
                                   mapped_box_b.getBlockId(),
                                   mapped_box_b.getPeriodicId());
      mapped_boxes_c.insert(mapped_box_c);
      b_eto_c[mapped_box_b.getId()].insert(mapped_box_c);
      c_eto_b[mapped_box_c.getId()].insert(mapped_box_b);
   }

   mapped_box_level_c.initialize( mapped_boxes_c,
                                  mapped_box_level_b.getRefinementRatio(),
                                  mapped_box_level_b.getGridGeometry(),
                                  mapped_box_level_b.getMPI() );

   b_to_c.initialize(mapped_box_level_b,
                     mapped_box_level_c,
                     hier::IntVector::getZero(dim),
                     b_eto_c);
   c_to_b.initialize(mapped_box_level_c,
                     mapped_box_level_b,
                     hier::IntVector::getZero(dim),
                     c_eto_b);

   b_to_c.checkConsistencyWithBase();
   b_to_c.checkConsistencyWithHead();
   c_to_b.checkConsistencyWithBase();
   c_to_b.checkConsistencyWithHead();

   hier::MappingConnectorAlgorithm mca;
   mca.assertMappingValidity(b_to_c);
   mca.assertMappingValidity(c_to_b);

   return;
}
