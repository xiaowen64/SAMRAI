//
// File:        Timer.I
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Timing wrapper for various routines.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


/*
 *************************************************************************** 
 *                                                                         *
 * Exclusive start and stop routines for timers.                           *
 *                                                                         *
 * For wallclock time: If we have MPI, we use MPI_Wtime to set the         *
 *                     start/stop point.  If we don't have MPI but do      *
 *                     have access to timer utilities in sys/times.h,      *
 *                     we use the time() utility to set the start/start    *
 *                     point.  If we have neither, we set the wallclock    *
 *                     start/stop time to zero.                            *
 *                                                                         *
 * For user time:      If we have access to timer utilities in sys/times.h,*
 *                     we use the times() utility to compute user and      *
 *                     system start/stop point (passing in the tms struct).*
 *                     If we don't have these utilities, we simply set the *
 *                     user and start/stop times to zero.                  *
 *                                                                         *
 * The timer manager manipulates the exclusive time information;  i.e.     *
 * when the timer's startExclusive and stopExclusive routines are called.  *
 *                                                                         *
 *************************************************************************** 
 */

inline void Timer::startExclusive()
{
   if (d_is_active) {

      Clock::timestamp(d_user_start_exclusive,
                            d_system_start_exclusive,
                            d_wallclock_start_exclusive);

   }
}

inline void Timer::stopExclusive()
{
   if (d_is_active) {
      Clock::timestamp(d_user_stop_exclusive,
                            d_system_stop_exclusive,
                            d_wallclock_stop_exclusive);

      d_wallclock_exclusive +=
         double(d_wallclock_stop_exclusive - d_wallclock_start_exclusive);
      d_user_exclusive += 
         double(d_user_stop_exclusive - d_user_start_exclusive);
      d_system_exclusive +=
         double(d_system_stop_exclusive - d_system_start_exclusive);
   }
}

/*
 *************************************************************************** 
 *                                                                         *
 * Simple utility routines to manipulate timers.                           *
 *                                                                         *
 ***************************************************************************
 */

inline const string &Timer::getName() const
{
   return(d_name);
}

inline int Timer::getIdentifier() const
{
   return(d_identifier);
}

inline bool Timer::isActive() const
{
   return(d_is_active);
}

inline bool Timer::isRunning() const
{
   return(d_is_running);
}

inline double Timer::getTotalUserTime() const
{
   return(d_user_total/Clock::getClockCycle());
}

inline double Timer::getTotalSystemTime() const
{
   return(d_system_total/Clock::getClockCycle());
}

inline double Timer::getTotalWallclockTime() const
{
#ifndef HAVE_MPI
   double clock_cycle = Clock::getClockCycle();
#else 
   double clock_cycle = 1.;
#endif
   return(d_wallclock_total/clock_cycle);
}

inline double Timer::getMaxWallclockTime() const
{
#ifndef HAVE_MPI
   double clock_cycle = Clock::getClockCycle();
#else 
   double clock_cycle = 1.;
#endif
   return(d_max_wallclock/clock_cycle);
}

inline double Timer::getExclusiveUserTime() const
{
   return(d_user_exclusive/Clock::getClockCycle());
}

inline double Timer::getExclusiveSystemTime() const
{
   return(d_system_exclusive/Clock::getClockCycle());
}

inline double Timer::getExclusiveWallclockTime() const
{
#ifndef HAVE_MPI
   double clock_cycle = Clock::getClockCycle();
#else 
   double clock_cycle = 1.;
#endif
   return(d_wallclock_exclusive/clock_cycle);
}

inline void Timer::setConcurrentTimer(const int id)
{
   d_concurrent_timers[id] = true;
}

inline bool* Timer::getConcurrentTimerVector() const
{
   return(d_concurrent_timers);
}

inline int Timer::getNumberAccesses() const
{
   return(d_accesses);
}

/*
 ***************************************************************************
 *                                                                         *
 * Private utility routines to manipulate timers in database.              *
 *                                                                         *
 ***************************************************************************
 */

inline void Timer::setInactive()
{
   d_is_active = false;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

