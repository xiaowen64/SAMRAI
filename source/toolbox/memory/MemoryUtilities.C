//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/MemoryUtilities.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Routines for tracking memory use in SAMRAI.
//

#include "tbox/MemoryUtilities.h"

#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include "tbox/IOStream.h"


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#ifdef HAVE_TAU
#if (PROFILING_ON || TRACING_ON)
#include <Profile/Profiler.h>
/* Register an "event" with Tau to track memory usage. */
TAU_PROFILE_STMT(TauUserEvent ue("memory use"))
#endif
#endif


namespace SAMRAI {
   namespace tbox {

double MemoryUtilities::s_max_memory = 0.;

/*
*************************************************************************
*                                                                       *
* Prints memory usage to specified output stream.  Each time this       *
* method is called, it prints in the format:                            *
*                                                                       *
*    253.0MB (265334688) in 615 allocs, 253.9MB reserved (871952 unused)*
*                                                                       *
* where                                                                 *
*                                                                       *
*    253.0MB is how many megabytes your current allocation has malloced.* 
*    2653346688 is the precise size (in bytes) of your current alloc.   *
*    615 is the number of items allocated with malloc.                  *
*    253.9MB is the current memory reserved by the system for mallocs.  *
*    871952 is the bytes currently not used in this reserved memory.    *
*                                                                       *
*************************************************************************
*/
void MemoryUtilities::printMemoryInfo(std::ostream& os) 
{
   /*
    * NOTE: This was taken directly from John Gyllenhal...
    */
   
   /* Get malloc info structure */
   struct mallinfo my_mallinfo = mallinfo();
   
   /* Get total memory reserved by the system for malloc currently*/
   double reserved_mem = my_mallinfo.arena;
   
   /* Get all the memory currently allocated to user by malloc, etc. */
   double used_mem = my_mallinfo.hblkhd + my_mallinfo.usmblks +
      my_mallinfo.uordblks;
   
   /* Get memory not currently allocated to user but malloc controls */
   double free_mem = my_mallinfo.fsmblks + my_mallinfo.fordblks;
   
   /* Get number of items currently allocated */
   double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks;

   /* Record high-water mark for memory used. */
   s_max_memory = MathUtilities<double>::Max(s_max_memory, used_mem);   

   /* Print out concise malloc info line */
   os << used_mem/(1024.0*1024.0) << "MB ("
      << used_mem << ") in "
      << number_allocated << " allocs, "
      << reserved_mem/(1024.0*1024.0) << "MB reserved ("
      << free_mem << " unused)" << std::endl;
}

/*
*************************************************************************
*                                                                       *
* Records memory usage to user-defined event in TAU.  Note that if TAU  *
* is not included, this method does nothing.                            *
*                                                                       *
*************************************************************************
*/

void MemoryUtilities::recordMemoryInfo(double time) 
{

#ifdef HAVE_TAU
   /*
    * Access information from mallinfo
    */   
   struct mallinfo my_mallinfo = mallinfo();
   
   /* Get total memory reserved by the system for malloc currently*/
   double reserved_mem = my_mallinfo.arena;
   
   /* Get all the memory currently allocated to user by malloc, etc. */
   double used_mem = my_mallinfo.hblkhd + my_mallinfo.usmblks +
      my_mallinfo.uordblks;
   
   /* Get memory not currently allocated to user but malloc controls */
   double free_mem = my_mallinfo.fsmblks + my_mallinfo.fordblks;
   
   /* Get number of items currently allocated */
   double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks;


   /* These vars are unused now but we may use them in the future */
   NULL_USE(reserved_mem);
   NULL_USE(free_mem);
   NULL_USE(number_allocated);

   /* Record high-water mark for memory used. */
   s_max_memory = MathUtilities<double>::Max(s_max_memory, used_mem);   

   /*
    * Record "used_mem" in MB to tau event.
    */
   TAU_PROFILE_STMT(ue.TriggerEvent(used_mem/(1024.0*1024.0)));
   
#endif
}

/*
*************************************************************************
*                                                                       *
* Prints maximum memory used (i.e. high-water mark).  The max is        *
* determined each time the "printMemoryInfo" or "recordMemoryInfo"      *
* functions are called.                                                 *
*                                                                       *
*************************************************************************
*/
void MemoryUtilities::printMaxMemory(std::ostream& os) 
{
   /*
    * Step through all nodes (>0) and send max memory to processor 0,
    * which subsequently writes it out.
    */
   int maxmem = 0;
   int len = 1;
   for (int p = 0; p < SAMRAI_MPI::getNodes(); p++) {
      if (SAMRAI_MPI::getRank() == p) {
         maxmem = (int)s_max_memory;
         SAMRAI_MPI::send(&maxmem, len, 0, false);
      }
      if (SAMRAI_MPI::getRank() == 0) {
         SAMRAI_MPI::recv(&maxmem, len, p, false);
      }
      os << "Maximum memory used on processor " << p
         << ": " << maxmem/(1024.*1024.) << " MB" << std::endl;
   }

}


}
}
