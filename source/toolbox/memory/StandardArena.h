//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/StandardArena.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Arena memory manager for standard allocation requests
//

#ifndef included_tbox_StandardArena
#define included_tbox_StandardArena

#include "SAMRAI_config.h"
#include "tbox/Arena.h"

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
