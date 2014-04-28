/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Test program to test the AssumedPartition classes
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/AssumedPartitionBox.h"

using namespace std;

using namespace SAMRAI;


struct CommonTestParams {
   hier::Box box;
   int rank_begin;
   int rank_end;
   int index_begin;
   CommonTestParams( const tbox::Dimension &dim ) :
      box(dim),
      rank_begin(0),
      rank_end(1),
      index_begin(0) {}
   CommonTestParams( const CommonTestParams &other ) :
      box(other.box),
      rank_begin(other.rank_begin),
      rank_end(other.rank_end),
      index_begin(other.index_begin) {}
};

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

      {
         /*
          * Test single-box assumed partitions.
          */

         int test_number = 0;

         while (true) {

            std::string test_name("SingleBoxTest");
            test_name += tbox::Utilities::intToString(test_number, 2);

            boost::shared_ptr<tbox::Database> test_db =
               input_db->getDatabaseWithDefault(test_name, boost::shared_ptr<tbox::Database>());

            if ( !test_db ) {
               break;
            }

            const std::string nickname =
               test_db->getStringWithDefault("nickname", test_name);

            tbox::plog << "\n\n\nStarting test " << test_name << " (" << nickname << ")\n";

            CommonTestParams ctp = getTestParametersFromDatabase( *test_db );

            hier::AssumedPartitionBox apb( ctp.box, ctp.rank_begin, ctp.rank_end, ctp.index_begin );
            fail_count += apb.selfCheck();

            ++test_number;

         }
      }


   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  transformation" << std::endl;
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
CommonTestParams getTestParametersFromDatabase( tbox::Database &test_db )
{
   const tbox::Dimension dim(test_db.getInteger("dim"));
   CommonTestParams ctp( dim );
   ctp.box = test_db.getDatabaseBox("box");
   ctp.rank_begin = test_db.getIntegerWithDefault("rank_begin", ctp.rank_begin);
   ctp.rank_end = test_db.getIntegerWithDefault("rank_end", ctp.rank_end);
   ctp.index_begin = test_db.getIntegerWithDefault("index_begin", ctp.index_begin);
   return ctp;
}
