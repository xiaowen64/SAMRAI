//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/branches/smith84/source/test/restartdb/mainHDF5.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2169 $
// Modified:    $LastChangedDate: 2008-04-30 16:24:39 -0700 (Wed, 30 Apr 2008) $
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

      tbox::Pointer<tbox::HDFDatabase> database = new tbox::HDFDatabase("SAMRAI Restart");      
      std::string name = "./restart." + tbox::Utilities::processorToString(tbox::SAMRAI_MPI::getRank()) + ".hdf5";
      hid_t file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, 
				H5P_DEFAULT, H5P_DEFAULT);
      hid_t hdf_group = H5Gcreate(file_id, "SAMRAIGroup", 0);
      database -> attachToFile(hdf_group);

      restart_manager -> setRootDatabase(database);

      restart_manager->writeRestartToDatabase();

      tbox::pout << "\n--- HDF5 write database tests END ---" << endl;

      tbox::pout << "\n--- HDF5 read database tests BEGIN ---" <<  endl;

      database -> close();

      restart_manager -> setRootDatabase(NULL);

      H5Fclose(file_id);

      database = new tbox::HDFDatabase("SAMRAI Restart");      
      file_id = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      hdf_group = H5Gopen(file_id, "SAMRAIGroup");
      database -> attachToFile(hdf_group);

      restart_manager -> setRootDatabase(database);

      hdf_tester.getFromDatabase();

      database -> close();

      restart_manager -> setRootDatabase(NULL);

      H5Fclose(file_id);

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

