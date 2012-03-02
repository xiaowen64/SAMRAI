/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   SAMRAI class to manage package startup and shutdown
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
void
SAMRAIManager::startup()
{
   TBOX_ASSERT(s_initialized);
   TBOX_ASSERT(!s_started);

   StartupShutdownManager::startup();

   s_started = true;
}

/*
 *************************************************************************
 *
 * Shutdown the SAMRAI package.  This routine currently only deallocates
 * statically allocated memory and finalizes the output streams.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
SAMRAIManager::shutdown()
{
   TBOX_ASSERT(s_initialized);
   TBOX_ASSERT(s_started);

   StartupShutdownManager::shutdown();

   s_started = false;

}

SAMRAI_INLINE_KEYWORD
void
SAMRAIManager::finalize()
{
   TBOX_ASSERT(s_initialized);

   StartupShutdownManager::finalize();

   PIO::finalize();

   s_initialized = false;
}

SAMRAI_INLINE_KEYWORD
bool
SAMRAIManager::isInitialized() {
   return s_initialized;
}

SAMRAI_INLINE_KEYWORD
bool
SAMRAIManager::isStarted() {
   return s_started;
}

/*
 *************************************************************************
 *
 * Functions to get and set the maximum number of patch data entries
 * that will be supported.  Note that the set routine cannot be called
 * after the max number has been accessed.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
SAMRAIManager::getMaxNumberPatchDataEntries()
{
   s_max_patch_data_entries_accessed = true;
   return s_max_patch_data_entries;
}

}
}
