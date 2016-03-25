//
// File:        $URL$
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2123 $
// Modified:    $LastChangedDate: 2008-04-08 16:33:49 -0700 (Tue, 08 Apr 2008) $
// Description: Tests Memory database in SAMRAI
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/DatabaseBox.h"
#include "tbox/Complex.h"
#include "tbox/MemoryDatabase.h"
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

      tbox::PIO::logAllNodes("Memorytest.log");

      tbox::pout << "\n--- Memory database tests BEGIN ---" << endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
      
      RestartTester memory_tester;

      tbox::pout << "\n--- Memory write database tests BEGIN ---" << endl;

      setupTestData();

      tbox::Pointer<tbox::MemoryDatabase> database = new tbox::MemoryDatabase("SAMRAI Restart");

      restart_manager -> setRootDatabase(database);
      
      restart_manager->writeRestartToDatabase();

      tbox::pout << "\n--- Memory write database tests END ---" << endl;

      tbox::pout << "\n--- Memory read database tests BEGIN ---" <<  endl;

      // In this test just read the database stored in memory that 
      // was just created.
      memory_tester.getFromDatabase();

      database -> close();

      tbox::pout << "\n--- Memory read database tests END ---" << endl;

      tbox::pout << "\n--- Memory database tests END ---" << endl;

      if (number_of_failures == 0) {
	 tbox::pout << "\nPASSED:  Memory" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(number_of_failures);

}

