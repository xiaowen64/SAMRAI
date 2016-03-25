//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/ReferenceCounter.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Reference counting class for Array and Pointer
//

#include "tbox/ReferenceCounter.h"
#include "tbox/Arena.h"
#include "tbox/ShutdownRegistry.h"

#include <stdlib.h>

#ifdef DEBUG_NO_INLINE
#include "tbox/ReferenceCounter.I"
#endif

namespace SAMRAI {
   namespace tbox {

ReferenceCounter *ReferenceCounter::s_free_list = NULL;
bool ReferenceCounter::s_registered_callback = false;
bool ReferenceCounter::s_shutdown = false;

void *ReferenceCounter::operator new(size_t bytes)
{
   if (s_free_list) {
      ReferenceCounter *node = s_free_list;
      s_free_list = s_free_list->d_counter;
      return(node);
   } else {
      return(::operator new(bytes));
   }
}

void ReferenceCounter::operator delete(void *what)
{
   ReferenceCounter *node = (ReferenceCounter *) what;
   node->d_counter = s_free_list;
   s_free_list = node;

   if (!s_registered_callback) {
      ShutdownRegistry::registerShutdownRoutine(freeCachedCopies,
		     ShutdownRegistry::priorityReferenceCounter);
      s_registered_callback = true;
   }
}

void ReferenceCounter::freeCachedCopies()
{
   while (s_free_list) {
      void *byebye = s_free_list;
      s_free_list = s_free_list->d_counter;
      ::operator delete(byebye);
   }

   s_registered_callback = false;

   s_shutdown = true;
}

ReferenceCounter::ReferenceCounter(Arena *arena,
                                             ReferenceCounter *counter)
{
   if (!arena) {
      d_references = 1;
      d_arena      = (Arena *) NULL;
      d_counter    = (ReferenceCounter *) NULL;
   } else {
      d_references = 1;
      d_arena      = arena;
      d_counter    = counter;
      if (d_counter) d_counter->d_references++;
   }
}

ReferenceCounter::~ReferenceCounter()
{
   if ((d_counter) && (--d_counter->d_references == 0)) {
      if (d_arena) {
         if (d_counter->d_arena) {
            d_arena->~Arena();
            d_counter->d_arena->free(d_arena);
         } else {
            delete d_arena;
         }
      }
      delete d_counter;
   }
}

}
}
