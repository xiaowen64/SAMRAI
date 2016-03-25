//
// File:	RestartManager.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An restart manager singleton class 
//

#include "tbox/RestartManager.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "tbox/HDFDatabase.h"
#include "tbox/MPI.h"
#include "tbox/NullDatabase.h"
#include "tbox/Parser.h"
#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"
#include <string>
using namespace std;
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

#ifndef NAME_BUFSIZE
#define NAME_BUFSIZE (32)
#endif

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
                 << "\n but Singleton instance already set." << endl);
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
   const string& root_dirname,
   const int restore_num,
   const int num_nodes)
{
   /* create the intermediate parts of the full path name of restart file */
   char restore_buf[NAME_BUFSIZE];
   char nodes_buf[NAME_BUFSIZE];
   char proc_buf[NAME_BUFSIZE];

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( NAME_BUFSIZE > (1 + 8 + 1 + 5 + 1) );
#endif
   sprintf(restore_buf,"/restore.%05d",restore_num);
   
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( NAME_BUFSIZE > (1 + 5 + 1 + 5 + 1) );
#endif
   sprintf(nodes_buf,"/nodes.%05d",num_nodes);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( NAME_BUFSIZE > (1 + 4 + 1 + 5 + 1) );
#endif
   int proc_num = MPI::getRank();
   sprintf(proc_buf,"/proc.%05d",proc_num);

   /* create full path name of restart file */
   string restart_filename = root_dirname + restore_buf + nodes_buf + proc_buf; 

   bool open_successful = true;
#ifdef HAVE_HDF5
   /* try to mount restart file */
   Pointer<HDFDatabase> temp_HDFDatabase = 
      new HDFDatabase(restart_filename);

   if (temp_HDFDatabase->mount(restart_filename, "R") < 0){
      TBOX_ERROR("Error attempting to open restart file " << restart_filename  
              << "\n   No restart file for processor: " << proc_num
	      << "\n   restart directory name = " << root_dirname
	      << "\n   number of processors   = " << num_nodes
	      << "\n   restore number         = " << restore_num << endl);
      open_successful = false;
   } else {
      /* set d_database root and d_is_from_restart */
      d_database_root = temp_HDFDatabase;
      d_is_from_restart = true;
   }

#else
   open_successful = false;
   TBOX_ERROR("Cannot open restart file...library not compiled with HDF" 
              << endl);
#endif

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
#ifdef HAVE_HDF5
   Pointer<HDFDatabase> temp_database = d_database_root;

   if (!temp_database.isNull()) {
      temp_database->unmount();
   }
#endif

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
   const string& name,
   Serializable* obj)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
   assert(obj != ((Serializable*)NULL));
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
         << "with restart manager." << endl);
   }
}

/*
*************************************************************************
*									*
* Removes the object with the specified name from d_restart_items_list. *
*									*
*************************************************************************
*/
void RestartManager::unregisterRestartItem(const string& name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
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
   const string& root_dirname, 
   int restore_num)
{
#ifdef HAVE_HDF5

   /* Create necessary directories and cd proper directory for writing */
   string restart_dirname = createDirs(root_dirname,restore_num);

   /* Create full path name of restart file */

   int proc_rank = MPI::getRank();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( NAME_BUFSIZE > (1 + 4 + 1 + 5 + 1) );
#endif
   char restart_filename_buf[NAME_BUFSIZE];
   sprintf(restart_filename_buf, "/proc.%05d", proc_rank);

   string restart_filename = restart_dirname + restart_filename_buf;

   HDFDatabase* new_restartDB = new HDFDatabase(restart_filename);
   new_restartDB->mount(restart_filename, "W");

   List<RestartManager::RestartItem>::Iterator 
                                                     i(d_restart_items_list);
   for ( ; i; i++) {
      Pointer<Database> obj_db = 
         new_restartDB->putDatabase(i().name);
      (i().obj)->putToDatabase(obj_db);
   }
 
   new_restartDB->unmount();

   delete new_restartDB;

#else
   TBOX_ERROR("Cannot write restart file...library not compiled with HDF" 
              << endl);
#endif
}

/*
*************************************************************************
*									*
* Creates the directory structure for the data files if they have not	*
* already been created.  						*
*									*
*************************************************************************
*/

string RestartManager::createDirs(
   const string& root_dirname,
   int restore_num)
{
   char restore_buf[NAME_BUFSIZE];
   char nodes_buf[NAME_BUFSIZE];

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( NAME_BUFSIZE> (1 + 8 + 1 + 5 + 1) );
#endif
   sprintf(restore_buf, "/restore.%05d", restore_num);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert( NAME_BUFSIZE > (1 + 5 + 1 + 5 + 1) );
#endif
   int num_procs = MPI::getNodes();
   sprintf(nodes_buf,"/nodes.%05d",num_procs);

   string full_dirname = root_dirname + restore_buf + nodes_buf;

   Utilities::recursiveMkdir(full_dirname);
  
   return full_dirname;
}


}
}
