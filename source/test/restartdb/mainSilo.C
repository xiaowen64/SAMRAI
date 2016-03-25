//
// File:        $URL$
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2117 $
// Modified:    $LastChangedDate: 2008-04-04 14:25:55 -0700 (Fri, 04 Apr 2008) $
// Description: Tests Silo database in SAMRAI
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/DatabaseBox.h"
#include "tbox/Complex.h"
#include "tbox/SiloDatabase.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include <string>

using namespace std;
using namespace SAMRAI;

#include "database_tests.h"

class RestartTester : public tbox::Serializable 
{
public:   

   RestartTester() 
   {
      tbox::RestartManager::getManager()->registerRestartItem("RestartTester",
                                                             this);
   }

   virtual ~RestartTester() {}

   void putToDatabase(tbox::Pointer<tbox::Database> db) 
   {
      writeTestData(db);
   }

   void getFromDatabase() 
   {
      tbox::Pointer<tbox::Database> root_db =
         tbox::RestartManager::getManager()->getRootDatabase();

      tbox::Pointer<tbox::Database> db;
      if (root_db->isDatabase("RestartTester")) {
         db = root_db->getDatabase("RestartTester");
      } 

      readTestData(db);
   }

};


int main(int argc, char *argv[]) 
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("Silotest.log");

#ifdef HAVE_SILO

      tbox::pout << "\n--- Silo database tests BEGIN ---" << endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
      
      RestartTester silo_tester;

      tbox::pout << "\n--- Silo write database tests BEGIN ---" << endl;

      setupTestData();

      tbox::Pointer<tbox::SiloDatabase> database = new tbox::SiloDatabase("SAMRAI Restart");

      database -> create("./restart." + tbox::Utilities::processorToString(tbox::SAMRAI_MPI::getRank()) + ".silo");

      restart_manager -> setRootDatabase(database);
      
      restart_manager->writeRestartToDatabase();

      database -> close();

      tbox::pout << "\n--- Silo write database tests END ---" << endl;

      tbox::pout << "\n--- Silo read database tests BEGIN ---" <<  endl;

      database = new tbox::SiloDatabase("SAMRAI Restart");

      database -> open("./restart." + tbox::Utilities::processorToString(tbox::SAMRAI_MPI::getRank()) + ".silo");

      restart_manager -> setRootDatabase(database);

      silo_tester.getFromDatabase();

      database -> close();

      tbox::pout << "\n--- Silo read database tests END ---" << endl;

      tbox::pout << "\n--- Silo database tests END ---" << endl;

#endif

      if (number_of_failures == 0) {
	 tbox::pout << "\nPASSED:  Silo" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(number_of_failures);

}

