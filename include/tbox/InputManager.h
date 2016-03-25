//
// File:	InputManager.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An input manager singleton class that parses input files
//

#ifndef included_tbox_InputManager
#define included_tbox_InputManager

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_InputDatabase
#include "tbox/InputDatabase.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif


namespace SAMRAI {
   namespace tbox {


/**
 * Class InputManager parses an input file and returns the associated
 * database.  This manager class hides the complexity of opening the
 * input file, creating the parser, and populating the database with
 * values.  The input manager is simple enough that it did not need
 * to be implemented as a singleton class; however, it was implemented
 * as a singleton to be consistent with the restart manager class.
 *
 * All processors must call the parsing routines.  Any errors are reported
 * to pout and will result in termination of the program.
 *
 * The input file will generally have the following format.  For each
 * object that requires user specified input, an entry of the following 
 * format should be included in the input file.
 *
 * \verbatim
 * Object_entry {
 *    keyword1 = <value1>     // maybe some end line comments
 *    keyword2 = <value2>
 *    nested_input_keyword {
 *       nested_data_keyword1 = <value1> 
 *    }
 *     \ldots
 * }
 * \endverbatim
 * 
 * For convenience, the input parser also supports C/C++ style comments,
 * "include" files, and some expression evaluation.
 *
 */

class InputManager
{
public:
   /**
    * Return a pointer to the single instance of the input manager.
    * All access to the input manager object is through getManager().
    */
   static InputManager *getManager();

   /**
    * Set a new input manager.  This routine can only be used by subclasses
    * since the constructor and destructor are protected.  The manager object
    * will be deleted at program exit.
    */
   static void setManager(InputManager *manager);

   /**
    * Deallocate the input manager instance.  It is not necessary to call
    * this routine at program termination, since it is automatically called
    * by the SAMRAI shutdown routines.
    */
   static void freeManager();

   /**
    * Return whether or not the manager has read an input database.  If 
    * so, it returns true.  If not, false.
    */
   static bool inputDatabaseExists();

   /**
    * Accessor method for the root input database held by InputManager.  
    * Inputs are read from the input file and held in this database. 
    * This method returns a pointer to the database, allowing any class 
    * in SAMRAI to access the information inside it using standard 
    * database calls.  For example, the following could appear in a 
    * SAMRAI class:
    *
    *       // get root database
    *       Pointer<Database> root_db =
    *          InputManager::getManager()->getInputDatabase();
    *       // get class's sub-database
    *       Pointer<Database> class_db = root_db->getDatabase("MyClass");
    *       // get parameter(s) from sub-database
    *       int dummy = class_db->getInteger("dummy");
    *      
    *
    * where "dummy" was supplied in "MyClass" entry of the input file:
    *
    *       MyClass {
    *          dummy = ...
    *       }
    *
    * This function is intended for SAMRAI classes in which there is no
    * easy or efficient way to supply input parameters.
    */
   static Pointer<Database> getInputDatabase();

   /**
    * Create a new database named "main" from the specified input file.
    */
   virtual Pointer<InputDatabase> parseInputFile(
      const string& filename);

   /**
    * Parse data from the specified file into the existing database.
    */
   virtual void parseInputFile(
      const string& filename, Pointer<InputDatabase> db);

protected:
   /**
    * The constructor is protected, since only subclasses of the singleton
    * may access the constructor.  Others may not explicitly create a
    * singleton class.
    */
   InputManager();

   /**
    * The destructor for the input manager is protected, since only the
    * singleton class and subclasses may destroy the manager objects.
    */
   virtual ~InputManager();

private:

   InputManager(const InputManager&); // not implemented
   void operator=(const InputManager&);    // not implemented

   static InputManager *s_manager_instance;
   static bool s_registered_callback;

   static Pointer<Database> s_input_db;
};


}
}

#endif
