//
// File:	ShutdownRegistry.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Registry of shutdown routines to be called at program exit
//

#ifndef included_tbox_ShutdownRegistry
#define included_tbox_ShutdownRegistry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

namespace SAMRAI {
   namespace tbox {

struct ShutdownRegistryItem;

/**
 * @brief Class ShutdownRegistry is a utility for managing callbacks 
 * at program completion.  
 *
 * Classes that wish to have a static function called
 * at program shutdown must register a class via registerShutdownRoutine().
 * The registered shutdown handler is then called at shutdown when the
 * callRegisteredShutdowns() is invoked.  Note that the application program
 * must explicitly call callRegisteredShutdowns().
 *
 * This class is useful in releasing allocated class-specific static memory
 * on program completion to simplify the search for memory leaks.
 */

struct ShutdownRegistry
{
   /**
    * Register a shutdown routine with the ShutdownRegistry.  The
    * shutdown handler must take no arguments and returns no value.
    * 
    * The priority is used to specify the order for shutdowns, higher
    * numbers are shutdown first, lower last (254 first, 0 last).
    * 
    * The priority is required since handlers can depend on one another.
    * Basic handlers like Lists should be last, higher level
    * handlers (which may free Lists) should be shutdown first.
    * 
    * Users are reserved priorities 254 to 127.  Unless there 
    * is a known dependency on a SAMRAI shutdown handler, users 
    * should use these priorities.
    *
    * There is a potential for an infinite loop if routines
    * register themselves with the registry as part of
    * the shutdown.  
    * 
    * @param callback  The function to call on shutdown
    * @param priority  The priority (254 called first, 0 last)
    */
   static void registerShutdownRoutine(void (*callback)(), 
				       unsigned char priority);

   /**
    * Invoke the registered shutdown handlers.  Note that this routine
    * must be explicitly called at the end of the application.  This
    * routine may only be called once, since the registered information
    * is discarded at the end of this call.
    */
   static void callRegisteredShutdowns();

   static const unsigned char priorityArenaManager            = 10;
   static const unsigned char priorityReferenceCounter        = 20;
   static const unsigned char priorityList                    = 30;
   static const unsigned char priorityVariableDatabase        = 40;
   static const unsigned char priorityInputManager            = 50;
   static const unsigned char priorityRestartManager          = 60;
   static const unsigned char priorityStatistician            = 60;
   static const unsigned char priorityBoundaryLookupTable     = 60;
   static const unsigned char priorityHierarchyDataOpsManager = 60;
   static const unsigned char priorityTimerManger             = 70;

private:

   static void initialize();

   static const short s_number_of_priorities = 255;
   static ShutdownRegistryItem* s_shutdown_list[s_number_of_priorities];
   static ShutdownRegistryItem* 
                                s_shutdown_list_last[s_number_of_priorities];
   static int s_num_shutdown_items[s_number_of_priorities];

   static bool s_initialized;

};


}
}

#endif
