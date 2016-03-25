//
// File:	SAMRAIManager.C
// Package:	SAMRAI initialization and shutdown
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 406 $
// Modified:	$Date: 2005-06-01 09:48:43 -0700 (Wed, 01 Jun 2005) $
// Description:	SAMRAI class to manage package startup and shutdown
//

#include "tbox/SAMRAIManager.h"
#include "tbox/IEEE.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#include <new>

namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*									*
* Set default number of patch data components to 256.                   *
* This value can be changed by calling the static function              *
* SAMRAIManager::setMaxNumberPatchDataEntries().  To avoid              *
* potentially errant behavior, the value should be set early            *
* on in the program execution (i.e., before constructing a patch        *
* hierarchy, or performing any operations involving variables).         *
* One set by the user it cannot be rest during program execution.       *
*									*
*************************************************************************
*/

int SAMRAIManager::s_max_patch_data_entries = 256;
bool SAMRAIManager::s_max_patch_data_entries_accessed = false;
bool SAMRAIManager::s_max_patch_data_entries_set_by_user = false;

/*
*************************************************************************
*									*
* Initialize the SAMRAI package.  This routine performs the following	*
* tasks:								*
*									*
* (1) Initialize the SAMRAI MPI package					*
* (2) Initialize the parallel I/O routines				*
* (3) Set up IEEE exception handlers					*
* (4) Set new handler so that an error message is printed if new fails. *
*									*
*************************************************************************
*/

static void badnew()
{
   TBOX_ERROR("operator `new' failed -- program abort!" << endl);
}

void SAMRAIManager::startup()
{
   MPI::initialize();
   PIO::initialize();
   IEEE::setupExceptionHandlers();

#ifndef LACKS_PROPER_MEMORY_HANDLER
   set_new_handler(badnew);
#endif
}

/*
*************************************************************************
*									*
* Shutdown the SAMRAI package.  This routine currently only deallocates	*
* statically allocated memory and finalizes the output streams.		*
*									*
*************************************************************************
*/

void SAMRAIManager::shutdown()
{
   ShutdownRegistry::callRegisteredShutdowns();
   PIO::finalize();
}

/*
*************************************************************************
*									*
* Functions to get ans set the maximum number of patch data components  *
* that will be supported.  Note that the set routine can only be called *
* once during program execution.                                        *
*									*
*************************************************************************
*/

int SAMRAIManager::getMaxNumberPatchDataEntries()
{
   s_max_patch_data_entries_accessed = true;
   return (s_max_patch_data_entries);
}

void SAMRAIManager::setMaxNumberPatchDataEntries(int maxnum)
{
   if (s_max_patch_data_entries_set_by_user) {
      TBOX_ERROR("SAMRAIManager::setMaxNumberPatchDataEntries() has already been \n"
                 << "called.  It cannot be called more than once! -- program abort!"
                 << endl);
   }
   if (s_max_patch_data_entries_accessed) {
      TBOX_ERROR("SAMRAIManager::getMaxNumberPatchDataEntries() has already been \n"
                 << "called and the value cannot be reset after that by calling \n"
                 << "SAMRAIManager::setMaxNumberPatchDataEntries() -- program abort!"
                 << endl);
   }
   s_max_patch_data_entries = ( (maxnum < 0) ? s_max_patch_data_entries : maxnum );
   s_max_patch_data_entries_set_by_user = true; 
}


}
}
