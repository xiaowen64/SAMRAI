//
// File:	ScratchArena.h
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
