//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/ArenaManager.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Arena memory management singleton class
//

#include "tbox/ArenaManager.h"
#include "tbox/Arena.h"
#include "tbox/FixedArena.h"
#include "tbox/Pointer.h"
#include "tbox/ScratchArena.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/StandardArena.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
   namespace tbox {

ArenaManager* ArenaManager::s_arena_instance = 
                   (ArenaManager*)NULL;
bool ArenaManager::s_registered_callback = false;

/*
*************************************************************************
*                                                                       *
* Static manager member functions.                                      *
*                                                                       *
*************************************************************************
*/

ArenaManager* ArenaManager::getManager()
{
   if (!s_arena_instance) {
      s_arena_instance = new ArenaManager;
   }
   if (!s_registered_callback) {
      ShutdownRegistry::registerShutdownRoutine(freeManager,
			     ShutdownRegistry::priorityArenaManager);
      s_registered_callback = true;
   }
   return(s_arena_instance);
}

void ArenaManager::freeManager()
{
   if (s_arena_instance) delete s_arena_instance;
   s_arena_instance = ((ArenaManager *) NULL);
}

void ArenaManager::registerSingletonSubclassInstance(
   ArenaManager* subclass_instance)
{
   if (!s_arena_instance) {
      s_arena_instance = subclass_instance;
      if (!s_registered_callback) {
         ShutdownRegistry::registerShutdownRoutine(freeManager,
			ShutdownRegistry::priorityArenaManager);
         s_registered_callback = true;
      }
   } else {
      TBOX_ERROR("ArenaManager internal error...\n"
                 << "Attemptng to set Singleton instance to subclass instance,"
                 << "\n but Singleton instance already set." << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Manager constructor and destructor.                                   *
*                                                                       *
*************************************************************************
*/

ArenaManager::ArenaManager()
{
}
   
ArenaManager::~ArenaManager()
{
}

Pointer<Arena> ArenaManager::getScratchAllocator()
{
   if (d_scratch_allocator.isNull()) {
      d_scratch_allocator = new ScratchArena;
   }
   return(d_scratch_allocator);
}

Pointer<Arena> ArenaManager::getStandardAllocator()
{
   if (d_standard_allocator.isNull()) {
      d_standard_allocator = new StandardArena;
   }
   return(d_standard_allocator);
}

Pointer<Arena>
ArenaManager::getFixedAllocator(const size_t bytes)
{
   if (d_fixed_allocator.isNull()) {
      d_fixed_allocator = new FixedArena(0);
   }
   return(d_fixed_allocator->allocateArena(bytes));
}

void ArenaManager::setFixedAllocator(
   const Pointer<FixedArena>& arena)
{
   d_fixed_allocator = arena;
}

void ArenaManager::setScratchAllocator(
   const Pointer<Arena>& arena)
{
   d_scratch_allocator = arena;
}

void ArenaManager::setStandardAllocator(
   const Pointer<Arena>& arena)
{
   d_standard_allocator = arena;
}


}
}
