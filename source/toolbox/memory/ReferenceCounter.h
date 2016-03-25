//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/ReferenceCounter.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Reference counting class for Array and Pointer
//

#ifndef included_tbox_ReferenceCounter
#define included_tbox_ReferenceCounter

#include "SAMRAI_config.h"
#ifndef included_stddef
#define included_stddef
#include <stddef.h>
#endif


namespace SAMRAI {
   namespace tbox {


class Arena;

/**
 * Class ReferenceCounter manages the shared reference counter and
 * arena resources used by Pointer and Array.  It uses a local
 * free pool of objects to speed memory allocation and deallocation.  The
 * locally cached free pool can be freed by calling freeCachedCopies().
 *
 * {\b Do not subclass!}  Changing the size of a ReferenceCounter
 * object will cause my simple memory allocation mechanism to break in
 * horrible and disgusting ways.
 *
 * @see tbox::Array
 * @see tbox::Pointer
 */

class ReferenceCounter
{
public:
   /**
    * Create a ReferenceCounter with an unmanaged memory arena.
    * The number of references is set to one.
    */
   ReferenceCounter();

   /**
    * Create a ReferenceCounter with a managed memory arena.  Argument
    * newArena is the managed memory arena and arenaCounter is the reference
    * counter for that arena.  The number of references is set to one.
    */
   ReferenceCounter(Arena *newArena,
                         ReferenceCounter *arenaCounter);

   /**
    * Destructor for ReferenceCounter.  The destructor releases
    * the managed memory arena if its count has gone to zero.
    */
   ~ReferenceCounter();

   /**
    * Get the managed memory arena (or NULL if none exists).
    */
   Arena *getArena();

   /**
    * Decrement the number of references.  True is returned if the
    * reference count has gone to zero; false otherwise.
    */
   bool deleteReference();

   /**
    * Increment the number of references.
    */
   void addReference();

   /**
    * Release the memory for all currently cached ReferenceCounter
    * copies.  This function may be called at any time.  In general, it
    * should not be necessary to call freeCachedCopies(), since it is
    * called via the ShutdownRegistry mechanism.
    */
   static void freeCachedCopies();

   /**
    * Class-specific operator new.  Data is allocated off of an
    * internal free list to speed memory allocation.
    */
   void *operator new(size_t bytes);

   /**
    * Class-specific operator delete.  Freed data is returned to
    * an internal free list for re-use by operator new.
    */
   void operator delete(void *what);

private:
   ReferenceCounter(const ReferenceCounter&);	// not implemented
   void operator=(const ReferenceCounter&);	// not implemented

   int d_references;
   Arena *d_arena;
   ReferenceCounter *d_counter;

   static ReferenceCounter *s_free_list;
   static bool s_registered_callback;

   static bool s_shutdown;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/ReferenceCounter.I"
#endif
#endif
