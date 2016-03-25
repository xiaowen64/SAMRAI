//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/ScratchArena.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Arena memory manager for scratch space
//

#ifndef included_tbox_ScratchArena
#define included_tbox_ScratchArena

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class ScratchArena is intended for temporary memory allocation
 * from a scratch arena.  However, ScratchArena currently uses the
 * standard C++ new and delete operators for memory allocation and
 * deallocation.
 *
 * @see tbox::Arena
 */

class ScratchArena : public Arena
{
public:
   /**
    * The constructor for the scratch memory arena.
    */
   ScratchArena();

   /**
    * The virtual destructor for the scratch memory arena.
    */
   virtual ~ScratchArena();

   /**
    * Allocate memory from the scratch arena.  The current implementation
    * of ScratchArena simply uses the standard C++ new operator.
    */
   virtual void *alloc(const size_t bytes);

   /**
    * Return memory to the scratch arena pool.  The current implementation
    * of ScratchArena simply uses the standard C++ delete operator.
    */
   virtual void free(void *p);

private:
   ScratchArena(const ScratchArena&);	// not implemented
   void operator=(const ScratchArena&);	// not implemented

};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/ScratchArena.I"
#endif
#endif
