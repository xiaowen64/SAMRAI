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
#include "SAMRAI/hier/MappedBoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappedBoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/MultiblockMappedBoxTree.h"
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
 * Partition boxes in the given MappedBoxLevel.  This method is meant
 * to partition and also to create a bunch of small boxes from the
 * user-input boxes in order to set up a non-trivial mesh
 * configuration.
 */
void partitionBoxes(
   hier::MappedBoxLevel &mapped_box_level,
   const hier::IntVector &max_box_size,
   const hier::IntVector &min_box_size);

/*
 * Shrink a MappedBoxLevel by the given amount.
 */
void shrinkMappedBoxLevel( hier::MappedBoxLevel &small_mapped_box_level,
                           const hier::MappedBoxLevel &big_mapped_box_level,
                           const hier::IntVector &shrinkage,
                           const tbox::Array<int> &unshrunken_blocks );

/*
 * Refine a MappedBoxLevel by the given ratio.
 */
void refineMappedBoxLevel( hier::MappedBoxLevel &mapped_box_level,
                           const hier::IntVector &ratio );

/*
************************************************************************
*
* This is an correctness test for the
* computeExternalPartsForMultiblock and
* computeExternalPartsForMultiblock method in
* MappedBoxLevelConnectorUtils.
*
* 1. Set up GridGeometry.
*
* 2. Build a big MappedBoxLevel.
*
* 3. Build a small MappedBoxLevel by shrinking the big one slightly.
*
* 4. Partition the two MappedBoxLevels.
*
* 5. Check the internal and external parts the big and small
*    MappedBoxLevels with respect to each other.
*
* Inputs:
*
*   GridGeometry { ... }
*
*   Main {
*      // Domain dimension
*      dim = 2
*
*      // Base name of output files generated.
*      base_name = "..."
*
*      // Logging option
*      log_all_nodes  = TRUE
*
*      // For breaking up user-specified domain boxes.
*      max_box_size = 7, 7
*      min_box_size = 2, 2
*
*      // Index space to exclude from the big MappedBoxLevel.
*      exclude2 = [(0,0,4), (19,19,7)]
*
*      // Blocks not to be shrunken when generating small MappedBoxLevel
*      // This allows some small blocks to touch the domain boundary.
*      unshrunken_boxes = 1, 3
*
*      // Refinement ratio of big MappedBoxLevel.
*      big_refinement_ratio = 2, 2
*
*      // Refinement ratio of small MappedBoxLevel.
*      small_refinement_ratio = 6, 6
*   }
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
   int fail_count = 0;

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


      /*
       * Generate the GridGeometry.
       */
      tbox::ConstPointer<hier::GridGeometry> grid_geometry;
      if (input_db->keyExists("GridGeometry")) {
         grid_geometry = new hier::GridGeometry(
            dim,
            "GridGeometry",
            tbox::Pointer<hier::TransferOperatorRegistry>(),
            input_db->getDatabase("GridGeometry"));
      } else {
         TBOX_ERROR("MappedBoxLevelConnectorUtils test: could not find entry GridGeometry"
                    << "\nin input.");
      }
      grid_geometry->printClassData(tbox::plog);

      /*
       * Empty blocks are blocks not to have any MappedBoxes.
       */
      tbox::Array<int> empty_blocks;
      if ( main_db->isInteger("empty_blocks") ) {
         empty_blocks = main_db->getIntegerArray("empty_blocks");
      }


      /*
       * Unshrunken blocks are blocks in which the small MappedBoxLevel
       * has the same index space as the big MappedBoxLevel.
       */
      tbox::Array<int> unshrunken_blocks;
      if ( main_db->isInteger("unshrunken_blocks") ) {
         unshrunken_blocks = main_db->getIntegerArray("unshrunken_blocks");
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      const hier::IntVector &one_vector(hier::IntVector::getOne(dim));
      const hier::IntVector &zero_vector(hier::IntVector::getZero(dim));

      /*
       * How much to shrink the big MappedBoxLevel to get the small one.
       */
      const hier::IntVector shrinkage(dim, 1);


      hier::MappedBoxLevelConnectorUtils mblcu;


      /*
       * Set up the domain.
       */

      hier::MappedBoxSet domain_mapped_boxes;
      grid_geometry->computePhysicalDomain(
         domain_mapped_boxes, one_vector );
      tbox::plog << "domain_mapped_boxes:\n"
                 << domain_mapped_boxes.format()
                 << std::endl;

      const hier::MultiblockMappedBoxTree domain_mapped_box_tree(
         grid_geometry, domain_mapped_boxes );


      /*
       * Construct the big_mapped_box_level.  It is a refinement of
       * the domain but without exclude* boxes.
       */

      hier::MappedBoxSet big_mapped_boxes;
      const std::string exclude("exclude");
      for ( int bn=0; bn<grid_geometry->getNumberBlocks(); ++bn ) {

         const hier::BlockId block_id(bn);

         const std::string exclude_boxes_name = exclude + tbox::Utilities::intToString(bn);
         if (main_db->keyExists(exclude_boxes_name)) {

            /*
             * Block bn boxes in big_mapped_box_level are the
             * block_domain \ exclude_boxes.
             */

            hier::BoxList block_domain(dim);
            grid_geometry->computePhysicalDomain(block_domain,
                                                 one_vector,
                                                 block_id);

            hier::BoxList exclude_boxes(dim);
            exclude_boxes = main_db->getDatabaseBoxArray(exclude_boxes_name);
            block_domain.removeIntersections(exclude_boxes);

            hier::LocalId last_local_id(-1);
            for ( hier::BoxList::Iterator bi(block_domain); bi; bi++ ) {
               big_mapped_boxes.insert(big_mapped_boxes.end(),
                                       hier::Box(*bi,
                                                       ++last_local_id,
                                                       0,
                                                       block_id));
            }

         }
         else {

            /*
             * Block bn boxes in big_mapped_box_level are the same as
             * block bn domain.
             */
            for ( hier::MappedBoxSetSingleBlockIterator bi(domain_mapped_boxes,block_id);
                  bi.isValid(); ++bi ) {
               big_mapped_boxes.insert(big_mapped_boxes.end(), *bi);
            }

         }

      }
      hier::MappedBoxLevel big_mapped_box_level(
         big_mapped_boxes,
         one_vector,
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld());

      const hier::MappedBoxSet &big_mapped_box_set(
         big_mapped_box_level.getMappedBoxes() );


      /*
       * Generate the "small" MappedBoxLevel by shrinking the big one
       * back at its boundary.
       */
      hier::MappedBoxLevel small_mapped_box_level(dim);
      shrinkMappedBoxLevel( small_mapped_box_level,
                            big_mapped_box_level,
                            shrinkage,
                            unshrunken_blocks );


      /*
       * Refine MappedBoxlevels as user specified.
       */
      if ( main_db->isInteger("big_refinement_ratio") ) {
         hier::IntVector big_refinement_ratio(dim);
         main_db->getIntegerArray("big_refinement_ratio", &big_refinement_ratio[0], dim.getValue());
         refineMappedBoxLevel( big_mapped_box_level,
                               big_refinement_ratio );
      }

      if ( main_db->isInteger("small_refinement_ratio") ) {
         hier::IntVector small_refinement_ratio(dim);
         main_db->getIntegerArray("small_refinement_ratio", &small_refinement_ratio[0], dim.getValue());
         refineMappedBoxLevel( small_mapped_box_level,
                               small_refinement_ratio );
      }


      /*
       * Partition the big and small MappedBoxLevels.
       *
       * Limit box sizes to make the configuration more complex than
       * the domain description.  Default is not to limit box size.
       */
      hier::IntVector max_box_size(dim,tbox::MathUtilities<int>::getMax());
      if ( main_db->isInteger("max_box_size") ) {
         main_db->getIntegerArray("max_box_size", &max_box_size[0], dim.getValue());
      }
      hier::IntVector min_box_size(dim,2);
      if ( main_db->isInteger("min_box_size") ) {
         main_db->getIntegerArray("min_box_size", &min_box_size[0], dim.getValue());
      }
      partitionBoxes(small_mapped_box_level, max_box_size, min_box_size);
      partitionBoxes(big_mapped_box_level, max_box_size, min_box_size);

      big_mapped_box_level.cacheGlobalReducedData();
      small_mapped_box_level.cacheGlobalReducedData();

      tbox::plog << "\nbig_mapped_box_level:\n"
                 << big_mapped_box_level.format("",2)
                 << '\n'
                 << "small_mapped_box_level:\n"
                 << small_mapped_box_level.format("",2)
                 << '\n'
         ;


      const hier::MappedBoxSet &small_mapped_box_set( small_mapped_box_level.getMappedBoxes() );

      const hier::MultiblockMappedBoxTree small_box_tree( grid_geometry, small_mapped_box_level.getGlobalizedVersion().getGlobalMappedBoxes() );



      /*
       * Connectors between big and small MappedBoxLevels.
       */

      const hier::Connector &small_to_big(
         small_mapped_box_level.getPersistentOverlapConnectors().createConnector(
            big_mapped_box_level,
            shrinkage) );
      small_to_big.cacheGlobalReducedData();

      const hier::Connector &big_to_small(
         big_mapped_box_level.getPersistentOverlapConnectors().createConnector(
            small_mapped_box_level,
            shrinkage) );
      big_to_small.cacheGlobalReducedData();

      tbox::plog << "\nsmall_to_big:\n"
                 << small_to_big.format("",2)
                 << '\n'
                 << "big_to_small:\n"
                 << big_to_small.format("",2)
                 << '\n';



      /*
       * Setup is complete.  Begin testing.
       */

      {
         /*
          * small_mapped_box_level nests inside big_mapped_box_level
          * by the shrinkage amount.  Verify that
          * computeInternalPartsForMultiblock finds all
          * small_mapped_box_level to be internal to
          * big_mapped_box_level and that
          * computeExternalPartsForMultiblock finds none of
          * small_mapped_box_level to be external to
          * big_mapped_box_level.
          *
          * small_mapped_box_level's internal parts should include
          * everything.  Thus small_to_everything should only map
          * small_mapped_box_level to its own index space (no more and
          * no less).
          *
          * small_mapped_box_level's external parts should include
          * nothing.  Thus small_to_nothing should map
          * small_mapped_box_level to nothing.
          */
         hier::MappedBoxLevel everything(dim), nothing(dim);
         hier::Connector small_to_everything, small_to_nothing;
         mblcu.computeExternalPartsForMultiblock(
            nothing,
            small_to_nothing,
            small_to_big,
            -shrinkage,
            domain_mapped_box_tree );
         mblcu.computeInternalPartsForMultiblock(
            everything,
            small_to_everything,
            small_to_big,
            -shrinkage,
            domain_mapped_box_tree );
         tbox::plog << "\nsmall_to_nothing:\n"
                    << small_to_nothing.format("",2) << '\n'
                    << "\nnothing:\n"
                    << nothing.format("",2) << '\n'
                    << "small_to_everything:\n"
                    << small_to_everything.format("",2) << '\n'
                    << "\neverything:\n"
                    << everything.format("",2) << '\n'
            ;

         for ( hier::MappedBoxSet::const_iterator bi=small_mapped_box_set.begin();
               bi!=small_mapped_box_set.end(); ++bi ) {
            const hier::Box &small_mapped_box = *bi;

            if ( small_to_everything.hasNeighborSet(small_mapped_box.getId()) ) {
               const hier::MappedBoxSet &neighbors = small_to_everything.getNeighborSet(small_mapped_box.getId());

               hier::BoxList neighbor_box_list;
               for ( hier::MappedBoxSet::const_iterator na=neighbors.begin();
                     na!=neighbors.end(); ++na ) {
                  neighbor_box_list.unionBoxes(*na);
                  if ( ! small_mapped_box.contains(*na) ) {
                     tbox::perr << "Mapping small_to_everyting erroneously mapped "
                                << small_mapped_box << " to:\n" << *na
                                << " which is outside itself.\n";
                     ++fail_count;
                  }
               }

               hier::BoxList tmp_box_list(small_mapped_box);
               tmp_box_list.removeIntersections(neighbor_box_list);
               if ( tmp_box_list.size() != 0 ) {
                  tbox::perr << "Mapping small_to_everything erroneously mapped "
                             << small_mapped_box << " to something less than itself:\n"
                             << neighbors.format();
               }

            }

            if ( small_to_nothing.hasNeighborSet(small_mapped_box.getId()) ) {
               const hier::MappedBoxSet &neighbors = small_to_nothing.getNeighborSet(small_mapped_box.getId());
               if ( ! neighbors.empty() ) {
                  tbox::perr << "Mapping small_to_nothing erroneously mapped " << small_mapped_box
                             << " to:\n" << neighbors.format()
                             << "\nIt should be mapped to nothing\n";
                  ++fail_count;
               }
            }
            else {
               tbox::perr << "Mapping small_to_nothing is missing a map from "
                          << small_mapped_box << " to nothing.\n";
               ++fail_count;
            }

         }
      }




      {
         /*
          * Compute the parts of big_mapped_box_level that are
          * internal to small_mapped_box_level and check for
          * correctness.
          *
          * To verify that the internal parts are correctly computed:
          *
          * - check that the small_mapped_box_level and the internal
          * parts of big_mapped_box_level have the same index space.
          */

         hier::MappedBoxLevel internal_mapped_box_level(dim);
         hier::Connector big_to_internal;
         mblcu.computeInternalPartsForMultiblock(
            internal_mapped_box_level,
            big_to_internal,
            big_to_small,
            zero_vector );
         const hier::MappedBoxSet &internal_mapped_box_set(internal_mapped_box_level.getMappedBoxes());
         tbox::plog << "internal_mapped_box_level:\n"
                    << internal_mapped_box_level.format("", 2)
                    << '\n'
                    << "big_to_internal:\n"
                    << big_to_internal.format("", 2);

         hier::MultiblockMappedBoxTree internal_box_tree(
            grid_geometry,
            internal_mapped_box_level.getGlobalizedVersion().getGlobalMappedBoxes() );

         for ( hier::MappedBoxSet::const_iterator ni=small_mapped_box_set.begin();
               ni!=small_mapped_box_set.end(); ++ni ) {
            hier::BoxList tmp_box_list(*ni);
            small_to_big.getHeadCoarserFlag() ?
               tmp_box_list.coarsen(small_to_big.getRatio()) :
               tmp_box_list.refine(small_to_big.getRatio());
            tmp_box_list.removeIntersections( ni->getBlockId(),
                                              big_mapped_box_level.getRefinementRatio(),
                                              internal_box_tree );
            if ( tmp_box_list.size() > 0 ) {
               tbox::perr <<"Small box " << *ni << " should fall within "
                          <<"the internal index space, but it doesn't." << std::endl;
               ++fail_count;
            }
         }

         for ( hier::MappedBoxSet::const_iterator ni=internal_mapped_box_set.begin();
               ni!=internal_mapped_box_set.end(); ++ni ) {
            hier::BoxList tmp_box_list(*ni);
            big_to_small.getHeadCoarserFlag() ?
               tmp_box_list.coarsen(big_to_small.getRatio()) :
               tmp_box_list.refine(big_to_small.getRatio());
            tmp_box_list.removeIntersections( ni->getBlockId(),
                                              small_mapped_box_level.getRefinementRatio(),
                                              small_box_tree );
            if ( tmp_box_list.size() > 0 ) {
               tbox::perr <<"Internal box " << *ni << " should fall within "
                          <<"the small index space, but it doesn't." << std::endl;
               ++fail_count;
            }
         }

      }


      {

         /*
          * Compute parts of big_mapped_box_level that are external to
          * small_mapped_box_level.
          *
          * Verify that the external parts are correctly computed:
          *
          * - check that external parts do not overlap small_mapped_box_level.
          *
          * - check that small_mapped_box_level does not overlap external parts.
          *
          * - check that big_mapped_box_level \ { small_mapped_box_level, external parts }
          *   is empty.
          */
         hier::MappedBoxLevel external_mapped_box_level(dim);
         hier::Connector big_to_external;
         mblcu.computeExternalPartsForMultiblock(
            external_mapped_box_level,
            big_to_external,
            big_to_small,
            zero_vector,
            hier::MultiblockMappedBoxTree() );
         const hier::MappedBoxSet &external_mapped_box_set(external_mapped_box_level.getMappedBoxes());
         tbox::plog << "\nexternal_mapped_box_level:\n"
                    << external_mapped_box_level.format("", 2)
                    << '\n'
                    << "big_to_external:\n"
                    << big_to_external.format("", 2);

         hier::MultiblockMappedBoxTree external_box_tree(
            grid_geometry,
            external_mapped_box_level.getGlobalizedVersion().getGlobalMappedBoxes() );

         for ( hier::MappedBoxSet::const_iterator ni=external_mapped_box_set.begin();
               ni!=external_mapped_box_set.end(); ++ni ) {
            hier::BoxList tmp_box_list(*ni);
            big_to_small.getHeadCoarserFlag() ?
               tmp_box_list.coarsen(big_to_small.getRatio()) :
               tmp_box_list.refine(big_to_small.getRatio());
            tmp_box_list.intersectBoxes( ni->getBlockId(),
                                         small_mapped_box_level.getRefinementRatio(),
                                         small_box_tree );
            if ( tmp_box_list.size() != 0 ) {
               tbox::perr <<"External box " << *ni << " should not\n"
                          <<"intersect small_mapped_box_level but does.\n"
                          <<"Intersections:\n";
               tmp_box_list.print(tbox::perr);
               tbox::perr << std::endl;
               ++fail_count;
            }
         }

         for ( hier::MappedBoxSet::const_iterator ni=small_mapped_box_set.begin();
               ni!=small_mapped_box_set.end(); ++ni ) {
            hier::BoxList tmp_box_list(*ni);
            small_to_big.getHeadCoarserFlag() ?
               tmp_box_list.coarsen(small_to_big.getRatio()) :
               tmp_box_list.refine(small_to_big.getRatio());
            tmp_box_list.intersectBoxes( ni->getBlockId(),
                                         big_mapped_box_level.getRefinementRatio(),
                                         external_box_tree );
            if ( tmp_box_list.size() != 0 ) {
               tbox::perr <<"Small box " << *ni << " should not intersect "
                          <<"the external parts but is does.\n"
                          <<"Intersections:\n";
               tmp_box_list.print(tbox::perr);
               tbox::perr << std::endl;
               ++fail_count;
            }
         }

         for ( hier::MappedBoxSet::const_iterator ni=big_mapped_box_set.begin();
               ni!=big_mapped_box_set.end(); ++ni ) {
            hier::BoxList tmp_box_list(*ni);
            big_to_small.getHeadCoarserFlag() ?
               tmp_box_list.coarsen(big_to_small.getRatio()) :
               tmp_box_list.refine(big_to_small.getRatio());
            tmp_box_list.removeIntersections( ni->getBlockId(),
                                              small_mapped_box_level.getRefinementRatio(),
                                              small_box_tree );
            small_to_big.getHeadCoarserFlag() ?
               tmp_box_list.coarsen(small_to_big.getRatio()) :
               tmp_box_list.refine(small_to_big.getRatio());
            tmp_box_list.removeIntersections( ni->getBlockId(),
                                              big_mapped_box_level.getRefinementRatio(),
                                              external_box_tree );
            if ( tmp_box_list.size() > 0 ) {
               tbox::perr <<"Big box " << *ni << " should not be inside "
                          <<"the small MappedBoxLevel and the external parts but is not.\n"
                          <<"Outside parts:\n";
               tmp_box_list.print(tbox::perr);
               tbox::perr << std::endl;
               ++fail_count;
            }
         }


      }




      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  MappedBoxLevelConnector test" << std::endl;
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
   return fail_count;
}





/*
 * Partition boxes in the given MappedBoxLevel.  This method is meant
 * to partition and also to create a bunch of small boxes from the
 * user-input boxes in order to set up a non-trivial mesh
 * configuration.
 */
void partitionBoxes(
   hier::MappedBoxLevel &mapped_box_level,
   const hier::IntVector &max_box_size,
   const hier::IntVector &min_box_size) {

   const tbox::Dimension &dim(mapped_box_level.getDim());

   hier::MappedBoxLevel domain_mapped_box_level(mapped_box_level);
   domain_mapped_box_level.setParallelState(hier::MappedBoxLevel::GLOBALIZED);

   mesh::TreeLoadBalancer load_balancer( mapped_box_level.getDim() );

   const tbox::Pointer<hier::PatchHierarchy> hierarchy;

   hier::Connector dummy_connector;

   const hier::IntVector bad_interval(dim,1);
   const hier::IntVector cut_factor(dim,1);

   load_balancer.loadBalanceMappedBoxLevel(
      mapped_box_level,
      dummy_connector,
      dummy_connector,
      tbox::Pointer<hier::PatchHierarchy>(),
      0,
      dummy_connector,
      dummy_connector,
      min_box_size,
      max_box_size,
      domain_mapped_box_level,
      bad_interval,
      cut_factor);

   return;
}




void shrinkMappedBoxLevel(
   hier::MappedBoxLevel &small_mapped_box_level,
   const hier::MappedBoxLevel &big_mapped_box_level,
   const hier::IntVector &shrinkage,
   const tbox::Array<int> &unshrunken_blocks )
{

   const tbox::ConstPointer<hier::GridGeometry> &grid_geometry(big_mapped_box_level.getGridGeometry());

   const int local_rank = big_mapped_box_level.getMPI().getRank();

   const hier::MappedBoxSet &big_mapped_box_set(big_mapped_box_level.getMappedBoxes());

   const hier::Connector &big_to_big(
      big_mapped_box_level.getPersistentOverlapConnectors().createConnector(
         big_mapped_box_level,
         shrinkage) );

   const hier::NeighborhoodSet &big_eto_big = big_to_big.getNeighborhoodSets();

   hier::MappedBoxSet visible_mapped_boxes(big_mapped_box_set);
   for ( hier::NeighborhoodSet::const_iterator mi=big_eto_big.begin();
         mi!=big_eto_big.end(); ++mi ) {
      visible_mapped_boxes.insert( mi->second.begin(), mi->second.end() );
   }

   std::map<hier::BlockId, hier::BoxList> boundary_boxes;
   for (hier::MappedBoxSet::const_iterator si = visible_mapped_boxes.begin();
        si != visible_mapped_boxes.end(); ++si) {
      boundary_boxes[si->getBlockId()].appendItem(*si);
   }

   hier::MultiblockMappedBoxTree visible_box_tree(
      grid_geometry,
      visible_mapped_boxes );

   hier::MappedBoxLevelConnectorUtils mblcu;

   mblcu.computeBoxesAroundBoundary(
      boundary_boxes,
      big_mapped_box_level.getRefinementRatio(),
      big_mapped_box_level.getGridGeometry() );

   tbox::plog << "shrinkMappedBoxLevel: Boundary plain boxes:\n";
   for ( std::map<hier::BlockId,hier::BoxList>::iterator mi=boundary_boxes.begin();
         mi!=boundary_boxes.end(); ++mi ) {
      tbox::plog << "Block " << mi->first << '\n';
      for ( hier::BoxList::Iterator bi(mi->second); bi; bi++ ) {
         tbox::plog << "  " << *bi << '\t' << (*bi).numberCells() << '\n';
      }
   }


   /*
    * Construct the complement of the small_mapped_box_level by
    * growing the boundary boxes.
    */

   std::vector<hier::Box> complement_mapped_boxes;

   hier::LocalId last_local_id(-1);
   for ( std::map<hier::BlockId,hier::BoxList>::iterator mi=boundary_boxes.begin();
         mi!=boundary_boxes.end(); ++mi ) {

      hier::BlockId block_id(mi->first);
      hier::BoxList &boundary_for_block(mi->second);

      for ( hier::BoxList::Iterator bi(boundary_for_block); bi; bi++ ) {
         hier::Box box(*bi);
         box.grow(shrinkage);
         hier::Box complement_mapped_box(
            box, ++last_local_id, local_rank, block_id );
         complement_mapped_boxes.push_back(complement_mapped_box);
      }

   }

   const hier::MultiblockMappedBoxTree complement_mapped_box_tree(
      grid_geometry,
      complement_mapped_boxes);


   /*
    * Construct the small_mapped_box_level.
    */

   last_local_id = -1;
   hier::MappedBoxSet small_mapped_boxes;
   for ( hier::MappedBoxSet::const_iterator bi=big_mapped_box_set.begin();
         bi!=big_mapped_box_set.end(); ++bi ) {

      const hier::Box &mapped_box = *bi;

      int ix;
      for ( ix=0; ix<unshrunken_blocks.size(); ++ix ) {
         if ( mapped_box.getBlockId() == unshrunken_blocks[ix] ) {
            break;
         }
      }

      if ( ix < unshrunken_blocks.size() ) {
         /*
          * This block should be excluded from shrinking.
          */
         small_mapped_boxes.insert( small_mapped_boxes.end(), mapped_box );
      }

      else {

         hier::BoxList shrunken_boxes(mapped_box);

         shrunken_boxes.removeIntersections( mapped_box.getBlockId(),
                                             big_mapped_box_level.getRefinementRatio(),
                                             complement_mapped_box_tree );
         shrunken_boxes.simplifyBoxes();

         for ( hier::BoxList::Iterator li(shrunken_boxes); li; li++ ) {
            const hier::Box shrunken_mapped_box(
               *li,
               ++last_local_id,
               mapped_box.getOwnerRank(),
               mapped_box.getBlockId() );

            small_mapped_boxes.insert( small_mapped_boxes.end(), shrunken_mapped_box );
         }
      }

   }

   small_mapped_box_level.swapInitialize(
      small_mapped_boxes,
      big_mapped_box_level.getRefinementRatio(),
      grid_geometry,
      big_mapped_box_level.getMPI() );

   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void refineMappedBoxLevel( hier::MappedBoxLevel &mapped_box_level,
                           const hier::IntVector &ratio )
{
   hier::MappedBoxSet refined_mapped_boxes;
   mapped_box_level.getMappedBoxes().refine(refined_mapped_boxes, ratio);
   mapped_box_level.swapInitialize(
      refined_mapped_boxes,
      ratio * mapped_box_level.getRefinementRatio(),
      mapped_box_level.getGridGeometry(),
      mapped_box_level.getMPI());
}
