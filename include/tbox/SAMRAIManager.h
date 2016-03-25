//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/SAMRAIManager.h $
// Package:     SAMRAI initialization and shutdown
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: SAMRAI class to manage package startup and shutdown
//

#ifndef included_tbox_SAMRAIManager
#define included_tbox_SAMRAIManager

#include "SAMRAI_config.h"
#include "tbox/Database.h"


namespace SAMRAI {
   namespace tbox {


/*!
 * @brief Class SAMRAIManager is a utility for managing startup and shutdown 
 * for SAMRAI applications and for changing the maximum number of patch data 
 * components supported by SAMRAI patches.  All applications should call 
 * SAMRAIManager::startup() (or SAMRAIManager::startup()) at the 
 * beginning of the program.  Startup should be called after initializing 
 * MPI but before any SAMRAI objects are used.  SAMRAIManager::shutdown()
 * (or SAMRAIManager:shutdown()) should be called near the end of the program, 
 * but before shutting down MPI and calling exit(0).  Note that the shutdown
 * function does not exit the program; it merely shuts down certain packages 
 * and deallocates memory (mostly objects with static members).
 */

struct SAMRAIManager
{
   /*!
    * Initialize the SAMRAI package.  Depending on the architecture and
    * compile flags, this routine sets up MPI, initializes IEEE exception
    * handlers, and other architecture-specific details.
    */
   static void startup();

   /*!
    * Shutdown the SAMRAI package.  Depending on the compile flags set at
    * compile-time, this routine shuts down MPI and calls registered shutdown
    * handlers.
    */
   static void shutdown();

   /*!
    * Return maximum number of patch data entries supported by SAMRAI.
    * The value is either the default value (256) or the value set by calling
    * the setMaxNumberPatchDataEntries() function.
    */
   static int getMaxNumberPatchDataEntries();

   /*!
    * Set maximum number of patch data entries supported by SAMRAI to the 
    * maximum of the current value and the argument value.
    * 
    * Note that this routine cannot be called anytime after the max patch 
    * data entries value has been accessed via the getMaxNumberPatchDataEntries()
    * function, either by the user or internally within SAMRAI.  Typically, the 
    * first internal access of this value occurs whenever any objects related 
    * to the patch hierarchy or variables are created. 
    */
   static void setMaxNumberPatchDataEntries(int maxnum);

   /*!
    * Return maximum number of timers supported by SAMRAI.
    * The value is either the default value (128) or the value set by
    * calling the setMaxNumberTimers() function.
    */
   static int getMaxNumberTimers();

   /*!
    * Set maximum number of timers supported by SAMRAI to the
    * maximum of the current value and the argument value.
    *
    * Note that this routine cannot be called anytime after the max
    * statistics value has been accessed via the getMaxTimers()
    * function, either by the user or internally within SAMRAI.  Typically,
    * the first internal access of this value occurs whenever the
    * timer manager is accessed.
    */
   static void setMaxNumberTimers(int maxnum);

   /*!
    * Return maximum number of statistics supported by SAMRAI.
    * The value is either the default value (128) or the value set by 
    * calling the setMaxNumberStatistics() function.
    */
   static int getMaxNumberStatistics();

   /*!
    * Set maximum number of statistics supported by SAMRAI to the
    * maximum of the current value and the argument value.
    *
    * Note that this routine cannot be called anytime after the max 
    * statistics value has been accessed via the getMaxStatistics()
    * function, either by the user or internally within SAMRAI.  Typically, 
    * the first internal access of this value occurs whenever the 
    * statistician is accessed. 
    */
   static void setMaxNumberStatistics(int maxnum);

private:
   static int s_max_patch_data_entries;
   static bool s_max_patch_data_entries_accessed;

   static int s_max_timers;
   static bool s_max_timers_accessed;

   static int s_max_statistics;
   static bool s_max_statistics_accessed;
};

}
}

#endif
