//
// File:	RestartManager.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An restart manager singleton class 
//

#ifndef included_tbox_RestartManager
#define included_tbox_RestartManager

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class RestartManager coordinates SAMRAI restart files (currently 
 * implemented using the HDF database class) and the objects comprising
 * SAMRAI-based application code.  The manager class orchestrates opening
 * and closing the database, stores data to be written out for restart, 
 * and writes out the restart data to the database.  Note that the restart
 * manager is a Singleton class that acts as a single point of control 
 * for restart capabilities.  As such its constructor and destructor are
 * protected members; they are not to be called outside of this class.
 *
 * The general procedure for starting a simulation from a restart file
 * is as follows.
 *
 * \begin{enumerate}
 * - Open the restart file using openRestartFile("filename").
 * - Get root of restart database using getRootDatabase().
 * - Initialize simulation objects using the restart constructor
 *       for the objects.
 * - Close the restart file using closeRestartFile().
 * \end{enumerate}
 *
 * Technically, there is no need to close the restart file because this will 
 * automatically be taken care of by the destructor for the database object.
 *
 * It is important to note in the initialization process, some objects
 * will need to be constructed in the "empty" state and filled in later
 * using some sort of getFromDatabase() method.
 *
 * The process for writing out state to a restart file is somewhat more
 * complicated.  The following things need to be taken care of.
 * 
 * \begin{enumerate}
 * - Each object that has state that needs to be saved for restart
 *       must be derived from the Serializable class (which
 *       responds to the putToDatabase() method).
 * - Any object that needs to save its state to the restart file 
 *       must be registered with the restart manager using the 
 *       registerRestartItem() method.   NOTE THAT NO TWO RESTARTABLE
 *       OBJECTS ARE ALLOWED TO HAVE THE SAME NAME STRING IDENTIFIER.
 * - The patchdata to be written to restart need to be specified 
 *       using the VariableDatabase::setPatchdataForRestart() method.  
 *       This is usually taken care of by the numerical algorithm object.
 * \end{enumerate}
 * 
 * When all these items are accounted for, writing to the restart file
 * is accomplished using a writeRestartFile() method.  There are
 * two writeRestartFile() methods available.  One takes only
 * a restart directory name as an argument whereas the other takes
 * both a restart directory name and a restore number for its arguments.
 * See comments for member functions for more details.
 *
 * @see tbox::Database
 */

class RestartManager
{
public:
   /**
    * Return a pointer to the single instance of the restart manager.
    * All access to the restart manager object is through getManager().
    *
    * Note that when the manager is accessed for the first time, the
    * Singleton instance is registered with the ShutdownRegistry
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    */
   static RestartManager *getManager();

   /**
    * Deallocate the restart manager instance.  It is not necessary to call
    * this routine at program termination, since it is automatically called
    * by the ShutdownRegistry class.
    */
   static void freeManager();

   /**
    * Returns true if the run is from a restart file (i.e. a restart file
    * has been opened from main()).  Returns false otherwise.
    */
   virtual bool isFromRestart();

   /**
    * Attempts to mount, for reading, the restart file for the processor.  
    * If there is no error opening the file, then the restart manager
    * mounts the restart file.  
    * Returns true if open is successful; false otherwise.
    */
   virtual bool openRestartFile(const string& root_dirname,
                                const int restore_num,
                                const int num_nodes);

   /** 
    * Closes the restart file.
    */
   virtual void closeRestartFile();

   /** 
    * Returns a Pointer to the root of the database.
    */
   virtual Pointer<Database> getRootDatabase();

   /**
    * Registers an object for restart with the given name.
    *
    * When assertion checking is active, an unrecoverable exception
    * will result if either the string is empty or the serializable
    * object pointer is null.
    */
   virtual void registerRestartItem(const string& name, 
                                    Serializable* obj);

   /**
    * Removes the object with the specified name from the list of
    * restartable items.
    *
    * When assertion checking is active, an unrecoverable exception
    * will result if the string is empty.
    */
   virtual void unregisterRestartItem(const string& name);

   /**
    * Clear all restart items managed by the restart manager.
    */
   virtual void clearRestartItems();

   /**
    * Write all objects registered to as restart objects to the 
    * restart database.  The string argument is the name of the
    * root of restart directory.
    * 
    * Note:  This method creates/uses a restart directory structure
    *    with 00000 as the restore number.
    */
   virtual void writeRestartFile(const string& root_dirname);

   /**
    * Write all objects registered to as restart objects to the 
    * restart database.  The string argument is the name of the
    * root of restart directory.  The integer argument is the 
    * identification number associated with the restart files generated.
    */
   virtual void writeRestartFile(const string& root_dirname, 
                                 const int restore_num);

protected:
   /**
    * The constructor for RestartManager is protected.
    * Consistent with the definition of a Singleton class, only the
    * manager object has access to the constructor for the class.
    * 
    * The constructor for RestartManager initializes the root 
    * data base to a NullDatabase and sets the restart flag to false.
    */
   RestartManager();

   /**
    * The destructor for the restart manager is protected, since only the
    * singleton class and subclasses may destroy the manager objects.
    */
   virtual ~RestartManager();

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    */
   void registerSingletonSubclassInstance(
      RestartManager* subclass_instance);

private:
   /* 
    * Create the directory structure for the data files.  
    * The directory structure created is       
    *
    *   restart_dirname/
    *     restore.[restore number]/
    *       nodes.[number of processors]/
    *         proc.[processor number]
    */
   string createDirs(const string& root_dirname, int restore_num);
 
   struct RestartItem {
      string name;
      Serializable* obj;
   };

   static RestartManager *s_manager_instance;
   static bool s_registered_callback;

   /* 
    * list of objects registered to be written to the restart database
    */
#ifdef LACKS_NAMESPACE_IN_DECLARE
   List< RestartItem > d_restart_items_list;
#else
   List< RestartManager::RestartItem > d_restart_items_list;
#endif

   Pointer<Database> d_database_root;
   bool d_is_from_restart;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/RestartManager.I"
#endif
#endif
