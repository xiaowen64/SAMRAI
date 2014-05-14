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
using namespace hier;


/*!
 * @brief Primitive BoxGenerator (independent of mesh package)
 * creating boxes using an AssumedPartition followed by an index
 * filter to keep a subset of boxes.
 */
struct PrimitiveBoxGen {
   boost::shared_ptr<hier::BaseGridGeometry> d_geom;
   hier::AssumedPartition d_ap;
   // Index filtering parameters.
   enum IndexFilter { ALL = 0	/* Keep all boxes */,
                      INTERVAL = 1	/* Keep d_num_keep, discard d_num_discard */,
                      LOWER = 2	/* Keep indices below d_frac */,
                      UPPER = 3	/* Keep indices above d_frac */
   };
   int d_index_filter;
   int d_num_keep;
   int d_num_discard;
   double d_frac;
   PrimitiveBoxGen(
      tbox::Database &database,
      const boost::shared_ptr<hier::BaseGridGeometry> &geom)
      {
         d_geom = geom;
         getFromInput(database);
      }
   PrimitiveBoxGen( const PrimitiveBoxGen &other ) :
      d_geom(other.d_geom),
      d_ap(other.d_ap),
      d_index_filter(other.d_index_filter),
      d_num_keep(other.d_num_keep),
      d_num_discard(other.d_num_discard),
      d_frac(other.d_frac) {}
   void getFromInput( tbox::Database &input_db );
   void getBoxes( hier::BoxContainer &boxes );
};


struct CommonTestParams {
   PrimitiveBoxGen d_boxes1;
   PrimitiveBoxGen d_boxes2;
   CommonTestParams( const tbox::Dimension &dim );
   CommonTestParams( const CommonTestParams &other );
};


CommonTestParams getTestParametersFromDatabase(
   tbox::Database &test_db );


void contriveConnector( Connector &conn,
                        const PrimitiveBoxGen &pb1,
                        const PrimitiveBoxGen &pb2,
                        const std::string &method,
                        tbox::Database &contrivance_db );


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

         PrimitiveBoxGen pb1( *input_db->getDatabase("PrimitiveBoxGen1"), grid_geom );
         BoxContainer boxes1;
         pb1.getBoxes(boxes1);
         hier::BoxLevel l1(boxes1, hier::IntVector::getOne(pb1.d_geom->getDim()), pb1.d_geom);

         PrimitiveBoxGen pb2( *input_db->getDatabase("PrimitiveBoxGen2"), grid_geom );
         BoxContainer boxes2;
         pb1.getBoxes(boxes2);
         hier::BoxLevel l2(boxes2, hier::IntVector::getOne(pb2.d_geom->getDim()), pb2.d_geom);

         /*
          * Rig up edges in l1_to_l2 by various contrivances and
          * compute transpose l2_to_l1.  Then check transpose
          * correctness.
          */
         char *contrivance_methods[] = { "mod", "near" };
         const size_t num_methods = 2;
         for ( size_t mi=0; mi<num_methods; ++mi ) {
            hier::Connector l1_to_l2(l1, l2, hier::IntVector::getZero(dim));
            contriveConnector( l1_to_l2, pb1, pb2, contrivance_methods[mi],
                               *input_db->getDatabase("ConnectorContrivance") );

            tbox::plog << "Testing with connector contrived by " << contrivance_methods[mi]
                       << "\nl1:\n" << l1.format("\t")
                       << "\nl2:\n" << l2.format("\t")
                       << "\nl1_to_l2:\n" << l1_to_l2.format("\t")
                       << std::endl;

            hier::Connector l2_to_l1(l2, l1, hier::IntVector::getZero(dim));
            l2_to_l1.setToTransposeOf(l1_to_l2);

            size_t test_fail_count = l1_to_l2.checkTransposeCorrectness(l2_to_l1);
            fail_count += test_fail_count;
            if ( test_fail_count ) {
               tbox::pout << "FAILED: " << test_fail_count << " with "
                          << contrivance_methods[mi] << "-contrived Connector.\n";
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
void contriveConnector( Connector &conn,
                        const PrimitiveBoxGen &pb1,
                        const PrimitiveBoxGen &pb2,
                        const std::string &method,
                        tbox::Database &contrivance_db )
{
   const int rank = conn.getBase().getMPI().getRank();

   if ( method == "mod" ) {
      const int denom = contrivance_db.getInteger("denom");
      for ( int i=pb1.d_ap.beginOfRank(rank); i<pb1.d_ap.endOfRank(rank); ++i ) {
         hier::Box l1box = pb1.d_ap.getBox(i);
         for ( int j=pb2.d_ap.begin(); j<pb2.d_ap.end(); ++j ) {
            hier::Box l2box = pb2.d_ap.getBox(i);
            if ( (i+j)%denom == 0 ) {
               conn.insertLocalNeighbor(l2box, l1box.getBoxId());
            }
         }
      }
   }

   else if ( method == "near" ) {
      int begin_shift = contrivance_db.getInteger("begin_shift");
      int end_shift = contrivance_db.getInteger("end_shift");
      int inc = contrivance_db.getInteger("inc");
      for ( int i=pb1.d_ap.beginOfRank(rank); i<pb1.d_ap.endOfRank(rank); ++i ) {
         hier::Box l1box = pb1.d_ap.getBox(i);

         int begin = l1box.getLocalId().getValue() + begin_shift;
         int end = l1box.getLocalId().getValue() + end_shift;
         begin = tbox::MathUtilities<int>::Max(begin, pb2.d_ap.begin());
         end = tbox::MathUtilities<int>::Min(end, pb2.d_ap.end());

         for ( int j=begin; j<end; j+=inc ) {
            hier::Box l2box = pb2.d_ap.getBox(i);
            conn.insertLocalNeighbor(l2box, l1box.getBoxId());
         }
      }
   }

   else {
      TBOX_ERROR("Contrivance method must be one of these: mod, forward, backward.");
   }

   return;
}




/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::getFromInput( tbox::Database &test_db )
{
   int rank_begin = 0;
   int rank_end = tbox::SAMRAI_MPI::getSAMRAIWorld().getSize();
   int index_begin = test_db.getIntegerWithDefault("index_begin", 0);
   d_ap.partition( d_geom->getPhysicalDomain(), rank_begin, rank_end, index_begin );

   std::string index_filter = test_db.getStringWithDefault("index_filter", "ALL");
   if ( index_filter == "ALL" ) {
      d_index_filter = PrimitiveBoxGen::ALL;
   }
   else if ( index_filter == "INTERVAL" ) {
      d_index_filter = PrimitiveBoxGen::INTERVAL;
   }
   else if ( index_filter == "LOWER" ) {
      d_index_filter = PrimitiveBoxGen::LOWER;
   }
   else if ( index_filter == "UPPER" ) {
      d_index_filter = PrimitiveBoxGen::UPPER;
   }
   return;
}




/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::getBoxes( hier::BoxContainer &boxes )
{
   const int rank = tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();

   if (d_index_filter == ALL) {
      int idbegin = d_ap.begin();
      int idend = d_ap.end();
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
         if ( interval_id < d_num_keep ) {
            boxes.push_back( d_ap.getBox(id) );
         }
      }
   }
   else if (d_index_filter == LOWER) {
      int threshold = d_ap.begin() + static_cast<int>( d_frac*(d_ap.end() - d_ap.begin()) );
      int idbegin = d_ap.begin();
      int idend = tbox::MathUtilities<int>::Max( threshold, d_ap.end() );
      for ( int id=idbegin; id<idend; ++id ) {
         boxes.push_back( d_ap.getBox(id) );
      }
   }
   else if (d_index_filter == UPPER) {
      int threshold = d_ap.begin() + static_cast<int>( d_frac*(d_ap.end() - d_ap.begin()) );
      int idbegin = tbox::MathUtilities<int>::Min( threshold, d_ap.begin() );
      int idend = d_ap.end();
      for ( int id=idbegin; id<idend; ++id ) {
         boxes.push_back( d_ap.getBox(id) );
      }
   }
   else {
      TBOX_ERROR("Invalid value of index_filter: " << d_index_filter);
   }

   return;
}

