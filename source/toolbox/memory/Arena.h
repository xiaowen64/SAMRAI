//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/Arena.h $
// Package:     SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Virtual base class for Arena memory management
//

#ifndef included_tbox_Arena
#define included_tbox_Arena

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_stddef
#define included_stddef
#include <stddef.h>
#endif
#ifndef included_new
#define included_new
#include <new>
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class Arena is an abstract base class for specialized dynamic
 * memory management.  The ``arena'' is a user-controlled portion of memory
 * that is accessed via a special ``new'' command.  Memory allocated in this
 * manner must be returned to the arena using the arena free() function;
 * C++ does not (unfortunately) understand a corresponding form of delete.
 * However, memory management classes such as Pointer and Array
 * can alleviate this burden.
 *
 * @see tbox::Array
 * @see tbox::Pointer
 */

class Arena : public DescribedClass
{
public:
   /**
    * The constructor for Arena.
    */
   Arena();

   /**
    * The virtual destructor for Arena.
    */
   virtual ~Arena();

   /**
    * Abstract virtual function to allocate memory from the arena.
    */
   virtual void *alloc(const size_t bytes) = 0;

   /**
    * Abstract virtual function to return memory to the arena memory pool.
    */
   virtual void free(void *p) = 0;

   /**
    * Static arena function to compute alignment for memory allocation.
    * Data allocations less than the alignment size are rounded up to the
    * next multiple of the allocation size.  All data allocations are
    * aligned on 16 byte boundaries.  Thus, a memory allocation of only
    * 9 bytes will actually return a 16 byte chunk of memory.
    */
   static size_t align(const size_t bytes);

private:
   enum { ArenaAllocationAlignment = 16 };

   Arena(const Arena&);       // not implemented
   void operator=(const Arena&);   // not implemented

};


}
}

#ifdef DEBUG_NO_INLINE
#define inline
#endif

inline void *operator new(size_t bytes, SAMRAI::tbox::Arena *arena);

#ifdef LACKS_NEW_PLACEMENT_OPERATOR
inline void *operator new(size_t bytes, void *ptr);
#endif
#ifdef DEBUG_NO_INLINE
#undef inline
#endif

#ifndef DEBUG_NO_INLINE
#include "tbox/Arena.I"
#endif
#endif
