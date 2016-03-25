//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/timers/Clock.h $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Simple utility class for interfacing with system clock
//

#ifndef included_tbox_Clock
#define included_tbox_Clock

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_systime
#ifdef _MSC_VER
#include <time.h>
#else
#include <sys/times.h>
#endif
#endif
#ifndef included_unistd
#ifndef _MSC_VER
#include <unistd.h>
#endif
#endif



namespace SAMRAI {
   namespace tbox {

/**
 * Class Clock serves as a single point of access for system clock 
 * information.  System and user time are computed via the POSIX compliant 
 * times() function.  This is described on p. 137, Lewine, POSIX programmers 
 * guide, 1992.  The methods and structs used in this utility are defined 
 * in <sys/times.h>.  Start and end times are stored as variables of type 
 * clock_t.  A clock_t value can be converted to seconds by dividing by 
 * CLK_TCK (which is defined in <sys/times.h>).  Different systems may use 
 * different CLK_TCK.  Time is accessed by calling the times() function which 
 * takes as an argument a reference to an object of type struct tms.  This 
 * object will record the system and user time (obj.tms_utime \& 
 * obj.tms_stime) and will return the time since the system was started.
 *
 * The return value from the call to times() can be used to compute elapsed
 * wallclock time.  Alternatively, one can use MPI_Wtime() if MPI libraries 
 * are included.  Two methods are defined for accessing system time - one that
 * has a clock_t struct argument for wallclock time (the non-MPI case) and 
 * one that has a double argument to record the value of MPI_Wtime().  
 *
 * Computing user/system/wallclock time with the times() function is performed
 * as follows:
 * \verbatim
 *    struct tms buffer;
 *    clock_t wtime_start = times(&buffer);
 *    clock_t stime_start = buffer.tms_stime;
 *    clock_t utime_start = buffer.tms_utime;
 *     (do some computational work)
 *    clock_t wtime_stop  = times(&buffer);
 *    clock_t stime_stop  = buffer.tms_stime;
 *    clock_t utime_stop  = buffer.tms_utime;
 *    double wall_time   = double(wtime_stop-wtime_start)/double(CLK_TCK);
 *    double user_time   = double(utime_stop-utime_start)/double(CLK_TCK);
 *    double sys_time    = double(stime_stop-stime_start)/double(CLK_TCK);
 * \endverbatim
 *
 */

struct Clock
{
   /**
    * Initialize system clock.  Argument must be in the "clock_t" format 
    * which is a standard POSIX struct provided on most systems in the 
    * <sys/times.h> include file. On Microsoft systems, it is provided in 
    * <time.h>.
    */
   static void initialize(clock_t& clock);

   /**
    * Initialize system clock, where clock is in double format.
    */
   static void initialize(double& clock);

   /**
    * Timestamp clocks for user, system, and wallclock times.  
    */
   static void timestamp(clock_t& user, clock_t& sys, clock_t& wall);

   /**
    * Timestamp user, system, and walltime clocks.  Wallclock argument is in 
    * double format since it will access wallclock times from MPI_Wtime() 
    * function.  
    */
   static void timestamp(clock_t& user, clock_t& sys, double& wall);

   /**
    * Returns clock cycle for the system.   
    */
   static double getClockCycle();

private:
#ifndef _MSC_VER
   static struct tms s_tms_buffer;
#endif
   static clock_t s_null_clock_t;
};

}
}

#endif





