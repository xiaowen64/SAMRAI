/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#include "SAMRAI/tbox/AllocatorDatabase.h"

#include "umpire/strategy/DynamicPool.hpp"
#include "umpire/ResourceManager.hpp"

#include "umpire/resource/MemoryResourceTypes.hpp"
#include "umpire/strategy/DynamicPool.hpp"
#include "umpire/strategy/AllocationAdvisor.hpp"

#include <iostream>

namespace SAMRAI {
namespace tbox {

#if defined(HAVE_UMPIRE)

AllocatorDatabase * AllocatorDatabase::s_allocator_database_instance(0);

StartupShutdownManager::Handler
AllocatorDatabase::s_startup_handler(
    0,
    AllocatorDatabase::startupCallback,
    0,
    0,
    tbox::StartupShutdownManager::priorityArenaManager);

void
AllocatorDatabase::startupCallback()
{
  AllocatorDatabase::getDatabase()->initialize();
}

void
AllocatorDatabase::shutdownCallback()
{
   if (s_allocator_database_instance) {
      delete s_allocator_database_instance;
   }
   s_allocator_database_instance = 0;
}

AllocatorDatabase *
AllocatorDatabase::getDatabase()
{
   if (!s_allocator_database_instance) {
      s_allocator_database_instance = new AllocatorDatabase();
   }
   return s_allocator_database_instance;
}

AllocatorDatabase::~AllocatorDatabase()
{
}

void
AllocatorDatabase::initialize()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  if (!rm.isAllocator("samrai::data_allocator")) {
#if defined(HAVE_CUDA)
    // Internal pool for allocations
#if 1
    auto allocator = rm.makeAllocator<umpire::strategy::AllocationAdvisor>(
        "internal::samrai::um_allocation_advisor",
        rm.getAllocator(umpire::resource::Unified),
        // Set preferred location to GPU
        "PREFERRED_LOCATION");
#endif
    //auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif

    rm.makeAllocator<umpire::strategy::DynamicPool>("samrai::data_allocator", allocator);
  }

  if (!rm.isAllocator("samrai::tag_allocator")) {
    rm.makeAllocator<umpire::strategy::DynamicPool>("samrai::tag_allocator",
        rm.getAllocator(umpire::resource::Host));
  }

  if (!rm.isAllocator("samrai::stream_allocator")) {
#if defined(HAVE_CUDA)
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif

    rm.makeAllocator<umpire::strategy::DynamicPool>("samrai::stream_allocator", allocator);
  }

  if (!rm.isAllocator("samrai::temporary_data_allocator")) {
#if defined(HAVE_CUDA)
    auto allocator = rm.getAllocator(umpire::resource::Device);
    //auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif
    rm.makeAllocator<umpire::strategy::DynamicPool>("samrai::temporary_data_allocator", allocator);
  }
}

umpire::Allocator
AllocatorDatabase::getDevicePool()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("samrai::temporary_data_allocator");
}

umpire::TypedAllocator<char>
AllocatorDatabase::getStreamAllocator()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return umpire::TypedAllocator<char>(rm.getAllocator("samrai::stream_allocator"));
}

umpire::Allocator
AllocatorDatabase::getTagAllocator()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("samrai::tag_allocator");
}

#endif

}
}
