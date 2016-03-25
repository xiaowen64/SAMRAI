//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/restartdb/mainHDF5.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2104 $
// Modified:    $LastChangedDate: 2008-04-01 10:02:41 -0700 (Tue, 01 Apr 2008) $
// Description: Tests HDF database in SAMRAI
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/DatabaseBox.h"
#include "tbox/Complex.h"
#include "tbox/HDFDatabase.h"
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

      tbox::PIO::logAllNodes("HDF5test.log");

#ifdef HAVE_HDF5

      tbox::pout << "\n--- HDF5 database tests BEGIN ---" << endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

      RestartTester hdf_tester;

      tbox::pout << "\n--- HDF5 write database tests BEGIN ---" << endl;

      setupTestData();

      restart_manager->writeRestartFile("test_dir", 0);

      tbox::pout << "\n--- HDF5 write database tests END ---" << endl;

      tbox::pout << "\n--- HDF5 read database tests BEGIN ---" <<  endl;

      restart_manager->closeRestartFile();

      restart_manager->openRestartFile("test_dir", 0, tbox::SAMRAI_MPI::getNodes());

      hdf_tester.getFromDatabase();

      restart_manager->closeRestartFile();

      tbox::pout << "\n--- HDF5 read database tests END ---" << endl;

      tbox::pout << "\n--- HDF5 database tests END ---" << endl;

#endif

      if (number_of_failures == 0) {
	 tbox::pout << "\nPASSED:  HDF5" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(number_of_failures);

}

