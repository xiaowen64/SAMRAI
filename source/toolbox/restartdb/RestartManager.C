//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/restartdb/RestartManager.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2129 $
// Modified:	$LastChangedDate: 2008-04-11 16:59:14 -0700 (Fri, 11 Apr 2008) $
// Description:	An restart manager singleton class 
//

#include <string>

#include "tbox/RestartManager.h"
#include "tbox/HDFDatabaseFactory.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/NullDatabase.h"
#include "tbox/Parser.h"
#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"

#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/RestartManager.I"
#endif

namespace SAMRAI {
   namespace tbox {

RestartManager* RestartManager::s_manager_instance = 
                     (RestartManager*)NULL;
bool RestartManager::s_registered_callback = false;

/*
*************************************************************************
*									*
* Basic singleton classes to create, set, and destroy the manager	*
* instance.								*
*									*
*************************************************************************
*/

RestartManager *RestartManager::getManager()
{
   if (!s_manager_instance) {
      s_manager_instance = new RestartManager;
   }
   if (!s_registered_callback) {
      ShutdownRegistry::registerShutdownRoutine(freeManager,
			ShutdownRegistry::priorityRestartManager);
      s_registered_callback = true;
   }
   return(s_manager_instance) ;
}

void RestartManager::freeManager()
{
   if (s_manager_instance) {
      s_manager_instance->clearRestartItems();
      delete s_manager_instance;
   }
   s_manager_instance = ((RestartManager *) NULL);
}

void RestartManager::registerSingletonSubclassInstance(
   RestartManager* subclass_instance)
{
   if (!s_manager_instance) {
      s_manager_instance = subclass_instance;
      if (!s_registered_callback) {
         ShutdownRegistry::registerShutdownRoutine(freeManager,
			ShutdownRegistry::priorityRestartManager);
         s_registered_callback = true;
      }
   } else {
      TBOX_ERROR("RestartManager internal error...\n"
                 << "Attemptng to set Singleton instance to subclass instance,"
                 << "\n but Singleton instance already set." << std::endl);
   }
}

/*
*************************************************************************
*									*
* The constructor and destructor are protected and call only be called	*
* by the singleton class or its subclasses.				*
*									*
*************************************************************************
*/

RestartManager::RestartManager()
{
   d_database_root = new NullDatabase();
#ifdef HAVE_HDF5
   d_database_factory = new HDFDatabaseFactory();
#else
   d_database_factory = NULL;
#endif
   d_is_from_restart = false;
   clearRestartItems();
}


/*
*************************************************************************
*									*
* Mount restart_file to the empty database created in the		*
* constructor and sets d_is_from_restart to true.  			*
* Return d_database_root.                                               *
*									*
*************************************************************************
*/

bool RestartManager::openRestartFile(
   const std::string& root_dirname,
   const int restore_num,
   const int num_nodes)
{
   int proc_num = SAMRAI_MPI::getRank();

   /* create the intermediate parts of the full path name of restart file */
   std::string restore_buf  = "/restore." + tbox::Utilities::intToString(restore_num, 6);
   std::string nodes_buf = "/nodes." + tbox::Utilities::nodeToString(num_nodes);
   std::string proc_buf = "/proc." + tbox::Utilities::processorToString(proc_num);

   /* create full path name of restart file */
   std::string restart_filename = root_dirname + restore_buf +
      nodes_buf + proc_buf; 

   bool open_successful = true;
   /* try to mount restart file */

   if(d_database_factory) {
   
      Pointer<Database> database = d_database_factory -> allocate(restart_filename);

      if (!database->open(restart_filename)){
	 TBOX_ERROR("Error attempting to open restart file " << restart_filename  
		    << "\n   No restart file for processor: " << proc_num
		    << "\n   restart directory name = " << root_dirname
		    << "\n   number of processors   = " << num_nodes
		    << "\n   restore number         = " << restore_num << std::endl);
	 open_successful = false;
      } else {
	 /* set d_database root and d_is_from_restart */
	 d_database_root = database;
	 d_is_from_restart = true;
      }
   } else {
      TBOX_ERROR("No DatabaseFactory supplied to RestartManager for opening " 
		 << restart_filename  << std::endl);
   }
      
   return (open_successful);
}

/*
*************************************************************************
*                                                                       *
* Closes the restart file by unmounting d_database_root and setting it  *
* to be a NullDatabase.                                            *
*                                                                       *
*************************************************************************
*/

void RestartManager::closeRestartFile()
{
   Pointer<Database> temp_database = d_database_root;

   if (!temp_database.isNull()) {
      temp_database->close();
   }

   d_database_root = new NullDatabase();
}

/*
*************************************************************************
*									*
* Registers the object for restart by adding it to			*
* d_restart_items_list.							*
*									*
*************************************************************************
*/
void RestartManager::registerRestartItem(
   const std::string& name,
   Serializable* obj)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(obj != ((Serializable*)NULL));
#endif
   /*
    * Run through list to see if there is another object registered
    * with the specified name.
    */
   List< RestartManager::RestartItem >::Iterator
      iter(d_restart_items_list);

   bool found_item = false;
   for ( ; !found_item && iter; iter++) {
      found_item = (iter().name == name);
   }

   /*
    * If there are no other items registered with the specified name,
    * add the object to the restart list.  Otherwise, throw an
    * error.
    */
   if (!found_item) {
      RestartItem r_obj;
      r_obj.name = name;
      r_obj.obj= obj;

      d_restart_items_list.appendItem(r_obj);

   } else {
      TBOX_ERROR("Register restart item error..."
		 << "\n   Multiple objects with name `" << name << "' registered "
		 << "with restart manager." << std::endl);
   }
}

/*
*************************************************************************
*									*
* Removes the object with the specified name from d_restart_items_list. *
*									*
*************************************************************************
*/
void RestartManager::unregisterRestartItem(const std::string& name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
#endif

   List< RestartManager::RestartItem >::Iterator 
      iter(d_restart_items_list);

   bool found_item = false;
   for ( ; !found_item && iter; iter++) {
      if (iter().name == name) {
         d_restart_items_list.removeItem(iter);
         found_item = true;
      }
   }
}

/*
*************************************************************************
*									*
* Remove all items from the restart item list.                          *
*									*
*************************************************************************
*/
void RestartManager::clearRestartItems()
{
   d_restart_items_list.clearItems();
}

/*
*************************************************************************
*									*
* Creates a new file with the given name and writes out the current 	*
* simulation state to the file by invoking the writeRestartFile() 	*
* method for all objects contained in	d_restart_objects_list.		*
*									*
*************************************************************************
*/
void RestartManager::writeRestartFile(
   const std::string& root_dirname, 
   int restore_num)
{
   /* Create necessary directories and cd proper directory for writing */
   std::string restart_dirname = createDirs(root_dirname,restore_num);

   /* Create full path name of restart file */

   int proc_rank = SAMRAI_MPI::getRank();

   std::string restart_filename_buf = 
      "/proc." + tbox::Utilities::processorToString(proc_rank);

   std::string restart_filename = restart_dirname + restart_filename_buf;

   if(d_database_factory) {

      Pointer<Database> new_restartDB = d_database_factory -> allocate(restart_filename);

      new_restartDB->create(restart_filename);

      writeRestartFile(new_restartDB);
 
      new_restartDB->close();

      new_restartDB.setNull();
   } else {
      TBOX_ERROR("No DatabaseFactory supplied to RestartManager for writeRestartFile " 
		 << restart_filename  << std::endl);
   }
}

/*
*************************************************************************
*									*
* Write simulation state to supplied database.                          *
*									*
*************************************************************************
*/
void RestartManager::writeRestartFile(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(database);
#endif

   List<RestartManager::RestartItem>::Iterator i(d_restart_items_list);
   for ( ; i; i++) {
      Pointer<Database> obj_db = 
         database->putDatabase(i().name);
      (i().obj)->putToDatabase(obj_db);
   }
}

/*
*************************************************************************
*									*
* Write simulation state to root database                               *
*									*
*************************************************************************
*/
void RestartManager::writeRestartToDatabase() 
{
   if(d_database_root) {
      writeRestartFile(d_database_root);
   } else {
      TBOX_ERROR("writeRestartToDatabase has no database to write to" 
		 << std::endl);
   }
}

/*
*************************************************************************
*									*
* Creates the directory structure for the data files if they have not	*
* already been created.  						*
*									*
*************************************************************************
*/

std::string RestartManager::createDirs(
   const std::string& root_dirname,
   int restore_num)
{
   int num_procs = SAMRAI_MPI::getNodes();

   std::string restore_buf = "/restore." + tbox::Utilities::intToString(restore_num, 6);
   std::string nodes_buf ="/nodes." + tbox::Utilities::processorToString(num_procs);

   std::string full_dirname = root_dirname + restore_buf + nodes_buf;

   tbox::Utilities::recursiveMkdir(full_dirname);
  
   return full_dirname;
}

}
}
