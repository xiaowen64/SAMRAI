//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/toolbox/base/SAMRAIManager.C $
// Package:	SAMRAI initialization and shutdown
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2223 $
// Modified:	$LastChangedDate: 2008-06-19 13:12:20 -0700 (Thu, 19 Jun 2008) $
// Description:	SAMRAI class to manage package startup and shutdown
//

#include "tbox/SAMRAIManager.h"
#include "tbox/IEEE.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"

#include <new>

namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*									*
* Set static members to set maximum number of patch data entries,       *
* statistics, and timers supported by the code.                         *
* These numbers are used to set the sizes of certain array containers   *
* used in the SAMRAI library.  They are set here so that they may be    *
* resized if necessary from a single access point (i.e., via the        *
* SAMRAIManager) if the default sizes (set here) are insufficient.      *
*                                                                       *
* These values can be changed by calling the static functions:          *
* SAMRAIManager::setMaxNumberPatchDataEntries(),                        *
*									*
* To avoid potentially erroneous or unexpected behavior, these          *
* value should be set early on during program execution before they     *
* are accessed in the library. Once accessed within the library,        *
* they cannot be reset during program execution.                        *
*									*
*************************************************************************
*/

int SAMRAIManager::s_max_patch_data_entries = 256;
bool SAMRAIManager::s_max_patch_data_entries_accessed = false;

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
   TBOX_ERROR("operator `new' failed -- program abort!" << std::endl);
}

void SAMRAIManager::startup()
{
   SAMRAI_MPI::initialize();
   PIO::initialize();
   IEEE::setupFloatingPointExceptionHandlers();

#ifndef LACKS_PROPER_MEMORY_HANDLER
   std::set_new_handler(badnew);
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
* Functions to get and set the maximum number of patch data entries     *
* that will be supported.  Note that the set routine cannot be called   *
* after the max number has been accessed.                               *
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
   if (s_max_patch_data_entries_accessed) {
      TBOX_ERROR("SAMRAIManager::setMaxNumberPatchDataEntries() error..."
                 << "\nThe max patch data entries value has already been accessed and cannot"
                 << "\nbe reset after that point by calling this method -- program abort!"
                 << std::endl);
   } else {
      s_max_patch_data_entries = MathUtilities<int>::Max(maxnum,
                                                         s_max_patch_data_entries);
   }
}

}
}
