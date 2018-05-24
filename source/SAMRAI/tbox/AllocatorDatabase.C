#include "SAMRAI/tbox/AllocatorDatabase.h"

#include "umpire/ResourceManager.hpp"

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
  rm.makeAllocator("SAMRAI_pool",
      "POOL",
      {0,0,0},
      {rm.getAllocator(umpire::resource::UnifiedMemory)});

  rm.makeAllocator("SAMRAI_device_pool",
      "POOL",
      {0,0,0},
      {rm.getAllocator(umpire::resource::Device)});

  rm.makeAllocator("SAMRAI_stream_pool",
      "POOL",
      {0,0,0},
      {rm.getAllocator(umpire::resource::PinnedMemory)});
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
