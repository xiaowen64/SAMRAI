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
  std::cout << "In startup Callback for AllocatorDatabase" << std::endl;
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
  std::cout << "Initializing database" << std::endl;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  /*
   * Register internal SAMRAI pool
   */
  auto um_alloc = rm.makeAllocator<umpire::strategy::AllocationAdvisor>(
      "SAMRAI_UM",
      rm.getAllocator(umpire::resource::UnifiedMemory),
      // Set preferred location to GPU
      "PREFERRED_LOCATION");

  rm.makeAllocator<umpire::strategy::DynamicPool>(
      "SAMRAI_pool",
      um_alloc);

  rm.makeAllocator<umpire::strategy::DynamicPool>(
      "SAMRAI_device_pool",
      rm.getAllocator(umpire::resource::Device));

  rm.makeAllocator<umpire::strategy::DynamicPool>(
      "SAMRAI_stream_pool",
      rm.getAllocator(umpire::resource::PinnedMemory));
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
