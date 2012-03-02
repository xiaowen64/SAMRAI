/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple utility class for interfacing with system clock
 *
 ************************************************************************/

#include "SAMRAI/tbox/SAMRAI_MPI.h"

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Initialize clock.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Clock::initialize(
   clock_t& clock)
{
#ifdef HAVE_SYS_TIMES_H
   clock = times(&s_tms_buffer);
#endif
}

SAMRAI_INLINE_KEYWORD
void
Clock::initialize(
   double& clock)
{
   clock = 0.0;
}

/*
 *************************************************************************
 *
 * Timestamp the provided structures with current system clock readings.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
Clock::timestamp(
   clock_t& user,
   clock_t& sys,
   double& wall)
{
#ifdef HAVE_SYS_TIMES_H
   s_null_clock_t = times(&s_tms_buffer);
   wall = SAMRAI_MPI::Wtime();
   sys = s_tms_buffer.tms_stime;
   user = s_tms_buffer.tms_utime;
#endif
}

/*
 *************************************************************************
 *
 * Get the clock cycle used by the system (time is then computed
 * as measured_time/clock_cycle)
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double
Clock::getClockCycle()
{
#ifdef _POSIX_VERSION
   double clock_cycle = double(sysconf(_SC_CLK_TCK));
#else
   double clock_cycle = 1.0;
#endif
   return clock_cycle;
}

}
}
