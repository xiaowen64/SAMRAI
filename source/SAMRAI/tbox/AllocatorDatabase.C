/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AllocatorDatabase.h"

#include "umpire/strategy/DynamicPool.hpp"
#include "umpire/ResourceManager.hpp"

#include "umpire/resource/MemoryResourceTypes.hpp"
#include "umpire/strategy/DynamicPool.hpp"
#include "umpire/strategy/AllocationAdvisor.hpp"

#include <iostream>

namespace SAMRAI {
namespace tbox {

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

AllocatorDatabase::AllocatorDatabase()
{
}

AllocatorDatabase::~AllocatorDatabase()
{
}

void
AllocatorDatabase::initialize()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  /*
   * Register internal SAMRAI pool
   */
  auto um_alloc = rm.makeAllocator<umpire::strategy::AllocationAdvisor>(
      "SAMRAI_UM",
      rm.getAllocator(umpire::resource::Unified),
      // Set preferred location to GPU
      "PREFERRED_LOCATION");

  rm.makeAllocator<umpire::strategy::DynamicPool>("SAMRAI_pool", um_alloc);

#if defined(HAVE_CUDA)
  auto device_allocator = rm.getAllocator(umpire::resource::Device);
  auto stream_allocator = rm.getAllocator(umpire::resource::Pinned);
#else
  auto device_allocator = rm.getAllocator(umpire::resource::Host);
  auto stream_allocator = rm.getAllocator(umpire::resource::Host);
#endif

  rm.makeAllocator<umpire::strategy::DynamicPool>("SAMRAI_device_pool", device_allocator);
  rm.makeAllocator<umpire::strategy::DynamicPool>("SAMRAI_stream_pool", stream_allocator);
}

umpire::Allocator
AllocatorDatabase::getDevicePool()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("SAMRAI_device_pool");
}

umpire::TypedAllocator<char>
AllocatorDatabase::getStreamAllocator()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return umpire::TypedAllocator<char>(rm.getAllocator("SAMRAI_stream_pool"));
}

}
}
