/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2017 Lawrence Livermore National Security, LLC
 * Description:   Singleton database class for managing variables and contexts.
 *
 ************************************************************************/

#ifndef included_tbox_AllocatorDatabase
#define included_tbox_AllocatorDatabase

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "umpire/Allocator.hpp"
#include "umpire/TypedAllocator.hpp"

namespace SAMRAI {
namespace tbox {

class AllocatorDatabase
{
public:
   static AllocatorDatabase* getDatabase();

   void initialize();

   umpire::Allocator
   getDevicePool();

   umpire::TypedAllocator<char>
   getStreamAllocator();

   umpire::Allocator
     getTagAllocator();

protected:
   AllocatorDatabase() = default;

   virtual ~AllocatorDatabase();

private:
   static void startupCallback();
   static void shutdownCallback();

   static AllocatorDatabase* s_allocator_database_instance;

   static StartupShutdownManager::Handler
   s_startup_handler;
};

}
}

#endif
