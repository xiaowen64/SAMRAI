/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An input manager singleton class that parses input files
 *
 ************************************************************************/

#include "SAMRAI/tbox/InputManager.h"
#include <stdlib.h>
#include <stdio.h>
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Parser.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

InputManager * InputManager::s_manager_instance = NULL;

Pointer<Database> InputManager::s_input_db =
   Pointer<Database>(NULL);

StartupShutdownManager::Handler InputManager::s_finalize_handler(
   0,
   0,
   0,
   InputManager::finalizeCallback,
   StartupShutdownManager::priorityInputManager);

/*
 *************************************************************************
 *									*
 * Basic singleton classes to create, set, and destroy the manager	*
 * instance.								*
 *									*
 *************************************************************************
 */

InputManager *InputManager::getManager()
{
   if (!s_manager_instance) {
      s_manager_instance = new InputManager;
   }
   return s_manager_instance;
}

void InputManager::setManager(
   InputManager* manager)
{
   if (s_manager_instance) {
      delete s_manager_instance;
   }
   s_manager_instance = manager;
}

void InputManager::finalizeCallback()
{
   if (s_manager_instance) {
      delete s_manager_instance;
      s_manager_instance = ((InputManager *)NULL);
   }

   s_input_db = ((InputDatabase *)NULL);
}

/*
 *************************************************************************
 *									*
 * The constructor and destructor are protected and call only be called	*
 * by the singleton class or its subclasses.				*
 *									*
 *************************************************************************
 */

InputManager::InputManager()
{
}

InputManager::~InputManager()
{
}

/*
 *************************************************************************
 *									*
 * Return whether or not the manager contains an valid input database.   *
 *									*
 *************************************************************************
 */

bool InputManager::inputDatabaseExists()
{
   return !(s_input_db.isNull());
}

/*
 *************************************************************************
 *									*
 * Parse the specified input file and return the new database.		*
 *									*
 *************************************************************************
 */

Pointer<InputDatabase>
InputManager::parseInputFile(
   const std::string& filename)
{
   Pointer<InputDatabase> db(new InputDatabase("main"));
   this->parseInputFile(filename, db);
   return db;
}

/*
 *************************************************************************
 *									*
 * Accessor method for InputManger's root input database.                *
 *									*
 *************************************************************************
 */
Pointer<Database> InputManager::getInputDatabase()
{
   return s_input_db;
}

/*
 *************************************************************************
 *									*
 * Parse the specified input file into the given database.		*
 *									*
 *************************************************************************
 */

void InputManager::parseInputFile(
   const std::string& filename,
   Pointer<InputDatabase> db)
{
   FILE* fstream = NULL;
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() == 0) {
      fstream = fopen(filename.c_str(), "r");
   }
   int worked = (fstream ? 1 : 0);
   mpi.Bcast(&worked, 1, MPI_INT, 0);
   if (!worked) {
      TBOX_ERROR("tbox::InputManager: Could not open input file``"
         << filename.c_str() << "''\n");
   }

   /*
    * Parse input file.
    */
   Parser* parser = new Parser();
   const int errors = parser->parse(filename, fstream, db);
   const int warnings = parser->getNumberWarnings();

   if (errors > 0) {
      TBOX_WARNING(
         "InputManager: Errors = " << errors
                                   << ", Warnings = " << warnings
                                   << "\n when parsing input file = "
                                   << filename << std::endl);
      db->printClassData(plog);
      TBOX_ERROR("InputManager exiting..." << std::endl);
   }
   if (warnings > 0) {
      TBOX_WARNING(
         "InputManager: Warnings  = " << warnings
                                      <<
         "\n when parsing input file = " << filename << std::endl);
   }

   /*
    * Store the root database in the static s_input_db variable.
    */
   s_input_db = db;

   delete parser;
   if (fstream) fclose(fstream);
}

}
}
