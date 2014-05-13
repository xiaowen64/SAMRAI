/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Test program to test Connector class
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/AssumedPartition.h"
#include "SAMRAI/geom/GridGeometry.h"

#include "boost/make_shared.hpp"

using namespace std;

using namespace SAMRAI;


/*!
 * @brief Primitive BoxGenerator (independent of mesh package)
 * creating boxes using an AssumedPartition followed by an index
 * filter to keep a subset of boxes.
 */
struct PrimitiveBoxGen {
   boost::shared_ptr<hier::BaseGridGeometry> d_geom;
   hier::AssumedPartition d_ap;
   // Assumed partition parameters.
   int d_rank_begin;
   int d_rank_end;
   int d_index_begin;
   enum IndexFilter { ALL	/* Keep all boxes */,
                      INTERVAL	/* Keep d_num_keep, discard d_num_discard */,
                      LOWER	/* Keep indices below d_frac */,
                      UPPER	/* Keep indices above d_frac */
   };
   // Index filtering parameters.
   int d_index_filter;
   int d_num_keep;
   int d_num_discard;
   double d_frac;
   PrimitiveBoxGen( const boost::shared_ptr<hier::BaseGridGeometry> &geom,
                    int rank_begin,
                    int rank_end,
                    int index_begin ) :
      d_geom(geom),
      d_ap(geom->getPhysicalDomain(), rank_begin, rank_end, index_begin),
      d_index_filter(ALL),
      d_num_keep(1),
      d_num_discard(0),
      d_frac(0.5) {}
   PrimitiveBoxGen( const PrimitiveBoxGen &other ) :
      d_geom(other.d_geom),
      d_ap(other.d_ap),
      d_index_filter(other.d_index_filter),
      d_num_keep(other.d_num_keep),
      d_num_discard(other.d_num_discard),
      d_frac(other.d_frac) {}
   void getBoxLevel( hier::BoxLevel &box_level );
};


struct CommonTestParams {
   PrimitiveBoxGen d_boxes1;
   PrimitiveBoxGen d_boxes2;
   CommonTestParams( const tbox::Dimension &dim );
   CommonTestParams( const CommonTestParams &other );
};


PrimitiveBoxGen getPrimitiveBoxGenFromInput(
   tbox::Database &test_db,
   const boost::shared_ptr<hier::BaseGridGeometry> &grid_geom);


CommonTestParams getTestParametersFromDatabase(
   tbox::Database &test_db );


int main(
   int argc,
   char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   int fail_count = 0;

   /*
    * Process command line arguments.  For each run, the input
    * filename must be specified.  Usage is:
    *
    * executable <input file name>
    */
   std::string input_filename;

   if (argc < 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input file> [case name]\n"
                            << "  options:\n"
                            << "  none at this time" << std::endl);
   } else {
      input_filename = argv[1];
   }

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      /*
       * Create input database and parse all data in input file.
       */
      boost::shared_ptr<tbox::MemoryDatabase> input_db(new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);


      boost::shared_ptr<tbox::Database> main_db = input_db->getDatabase("Main");

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Start logging.
       */
      const std::string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         tbox::PIO::logAllNodes(log_file_name);
      } else {
         tbox::PIO::logOnlyNodeZero(log_file_name);
      }

      const int rank = tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();

      {

         const tbox::Dimension dim(main_db->getInteger("dim"));

         if ( !input_db->isDatabase("BlockGeometry") ) {
            TBOX_ERROR("getTestParametersFromDatabase: You must specify \"BlockGeometry\" in input database.");
         }
         // Note: Using GridGeometry only because BaseGridGeometry can't be instanstiated.
         boost::shared_ptr<hier::BaseGridGeometry> grid_geom =
            boost::make_shared<geom::GridGeometry>(
               dim,
               "BlockGeometry",
               input_db->getDatabase("BlockGeometry"));

         PrimitiveBoxGen pb1 = getPrimitiveBoxGenFromInput(
            *input_db->getDatabase("PrimitiveBoxGen1"),
            grid_geom );
         hier::BoxLevel l1(hier::IntVector::getOne(pb1.d_geom->getDim()), pb1.d_geom);
         pb1.getBoxLevel(l1);

         PrimitiveBoxGen pb2 = getPrimitiveBoxGenFromInput(
            *input_db->getDatabase("PrimitiveBoxGen2"),
            grid_geom );
         hier::BoxLevel l2(hier::IntVector::getOne(pb2.d_geom->getDim()), pb2.d_geom);
         pb2.getBoxLevel(l2);

         /*
          * Rig up edges in l1_to_l2 by various contrivances and
          * compute transpose l2_to_l1.  Then check transpose
          * correctness.
          */

         {
            hier::Connector l1_to_l2(l1, l2, hier::IntVector::getZero(dim));
            for ( int i=pb1.d_ap.beginOfRank(rank); i<pb1.d_ap.endOfRank(rank); ++i ) {
               hier::Box l2box = pb2.d_ap.getBox(i);
               hier::Box l1box = pb1.d_ap.getBox(i);
               l1_to_l2.insertLocalNeighbor(l2box, l1box.getBoxId());
            }

            hier::Connector l2_to_l1(l2, l1, hier::IntVector::getZero(dim));
            l2_to_l1.setToTransposeOf(l1_to_l2);

            size_t test_fail_count = l1_to_l2.checkTransposeCorrectness(l2_to_l1);
            fail_count += test_fail_count;
            if ( test_fail_count ) {
               tbox::pout << "FAILED: " << test_fail_count << " problems in overlap correctness\n";
            }
         }


      }


   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  Connector" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}




/*
 *************************************************************************
 *************************************************************************
 */
PrimitiveBoxGen getPrimitiveBoxGenFromInput(
   tbox::Database &test_db,
   const boost::shared_ptr<hier::BaseGridGeometry> &grid_geom)
{
   int rank_begin = test_db.getIntegerWithDefault("rank_begin", 0);
   int rank_end = test_db.getIntegerWithDefault("rank_end", tbox::SAMRAI_MPI::getSAMRAIWorld().getSize());
   int index_begin = test_db.getIntegerWithDefault("index_begin", 0);
   int index_filter = test_db.getIntegerWithDefault("index_filter", PrimitiveBoxGen::ALL);

   return PrimitiveBoxGen( grid_geom, rank_begin, rank_end, index_begin );
}




/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::getBoxLevel( hier::BoxLevel &box_level )
{
   hier::BoxContainer boxes;

   const int rank = tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();

   if (d_index_filter == ALL) {
      int idbegin = d_ap.beginOfRank( rank );
      int idend = d_ap.endOfRank( rank );
      for ( int id=idbegin; id<idend; ++id ) {
         boxes.push_back( d_ap.getBox(id) );
      }
   }
   else if (d_index_filter == INTERVAL) {
      int idbegin = d_ap.begin();
      int idend = d_ap.end();
      int interval = d_num_keep + d_num_discard;
      for ( int id=idbegin; id<idend; ++id ) {
         int interval_id = id%interval;
         if ( interval_id < d_num_keep && d_ap.getOwner(id) == rank ) {
            boxes.push_back( d_ap.getBox(id) );
         }
      }
   }
   else if (d_index_filter == LOWER) {
      int threshold = d_ap.begin() + static_cast<int>( d_frac*(d_ap.end() - d_ap.begin()) );
      int idbegin = d_ap.beginOfRank( rank );
      int idend = tbox::MathUtilities<int>::Max( threshold, d_ap.endOfRank( rank ) );
      for ( int id=idbegin; id<idend; ++id ) {
         boxes.push_back( d_ap.getBox(id) );
      }
   }
   else if (d_index_filter == UPPER) {
      int threshold = d_ap.begin() + static_cast<int>( d_frac*(d_ap.end() - d_ap.begin()) );
      int idbegin = tbox::MathUtilities<int>::Min( threshold, d_ap.beginOfRank( rank ) );
      int idend = d_ap.endOfRank( rank );
      for ( int id=idbegin; id<idend; ++id ) {
         boxes.push_back( d_ap.getBox(id) );
      }
   }

   box_level.swapInitialize( boxes, box_level.getRefinementRatio(), d_geom );

   return;
}

