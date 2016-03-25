//
// File:	StandardArena.h
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Arena memory manager for standard allocation requests
//

#ifndef included_tbox_StandardArena
#define included_tbox_StandardArena

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class StandardArena is intended for standard memory allocation
 * requests.  It currently uses the standard C++ new and delete operators
 * for memory allocation and deallocation.
 *
 * @see tbox::Arena
 */

class StandardArena : public Arena
{
public:
   /**
    * The constructor for the standard memory arena.
    */
   StandardArena();

   /**
    * The virtual destructor for the standard memory arena.
    */
   virtual ~StandardArena();

   /**
    * Allocate memory from the standard memory arena.  The current
    * implementation of StandardArena uses the standard C++ new
    * operator.
    */
   virtual void *alloc(const size_t bytes);

   /**
    * Return memory to the standard arena pool.  The current
    * implementation of StandardArena uses the standard
    * C++ delete operator.
    */
   virtual void free(void *p);

private:
   StandardArena(const StandardArena&);	// not implemented
   void operator=(const StandardArena&);		// not implemneted
};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/StandardArena.I"
#endif
#endif
