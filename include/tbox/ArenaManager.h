//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/ArenaManager.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Arena memory management singleton class
//

#ifndef included_tbox_ArenaManager
#define included_tbox_ArenaManager

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_FixedArena
#include "tbox/FixedArena.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class ArenaManager is a Singleton class that acts as an abstract
 * factory for the arena memory management classes.  It is a single point
 * of control for memory allocation.  Class ArenaManager supports three
 * forms of memory management:  standard allocators, scratch allocators, and
 * fixed allocators.  Standard memory allocators are the ones to be used in
 * most cases; by default, the standard memory allocators use new and delete.
 * The scratch allocators are for temporary scratch space.  The fixed size
 * memory allocators are arenas of a fixed size from which many memory
 * requests will be made.
 *
 * See the Design Patterns book by Gamma {\em et al.} for more information
 * about the singleton pattern.
 *
 * @see tbox::Arena
 * @see tbox::FixedArena
 */

class ArenaManager
{
public:
   /**
    * Return a pointer to the single instance of the memory arena manager.
    * All access to the ArenaManager object is through getManager.
    * For example, to get a pointer to the standard allocator arena, we
    * must use ArenaManager::getManager()->getStandardAllocator().
    *
    * Note that when the manager is accessed for the first time, the
    * Singleton instance is registered with the ShutdownRegistry
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    */
   static ArenaManager *getManager();

   /**
    * Deallocate the ArenaManager instance.  It is not necessary
    * to call freeManager() at program termination, since it is automatically
    * called by the ShutdownRegistry class.
    */
   static void freeManager();

   /**
    * Allocate a fixed size arena of the given size.  This arena may then
    * be used in multiple memory allocation calls, although only one call
    * to new will have been made.
    */
   virtual Pointer<Arena> getFixedAllocator(const size_t bytes);

   /**
    * Return a pointer to scratch memory allocation arena.  This arena is
    * to be used for short-lived memory objects that should be allocated
    * from a scratch space.
    */
   virtual Pointer<Arena> getScratchAllocator();

   /**
    * Return a pointer to the standard memory allocation arena.  Use this
    * arena for most standard memory allocation needs.
    */
   virtual Pointer<Arena> getStandardAllocator();

   /**
    * Define a new fixed memory arena allocator.
    */
   virtual void setFixedAllocator(const Pointer<FixedArena>& arena);

   /**
    * Define a new scratch memory arena allocator.
    */
   virtual void setScratchAllocator(const Pointer<Arena>& arena);

   /**
    * Define a new standard memory arena allocator.
    */
   virtual void setStandardAllocator(const Pointer<Arena>& arena);

protected:
   /**
    * The constructor for ArenaManager is protected.  Consistent
    * with the definition of a Singleton class, only subclasses have
    * access to the constructor for the class.
    */
   ArenaManager();

   /**
    * The destructor for ArenaManager is protected.  See the
    * comments for the constructor.
    */
   virtual ~ArenaManager();

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    */
   void registerSingletonSubclassInstance(
      ArenaManager* subclass_instance);

private:
   static ArenaManager *s_arena_instance;
   static bool s_registered_callback;

   Pointer<FixedArena> d_fixed_allocator;
   Pointer<Arena> d_scratch_allocator;
   Pointer<Arena> d_standard_allocator;
};

}
}

#endif
