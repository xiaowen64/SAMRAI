/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Timing wrapper for various routines.
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 ***************************************************************************
 *
 * Exclusive start and stop routines for timers.
 *
 * For wallclock time: If we have MPI, we use MPI_Wtime to set the
 *                     start/stop point.  If we don't have MPI but do
 *                     have access to timer utilities in sys/times.h,
 *                     we use the time() utility to set the start/start
 *                     point.  If we have neither, we set the wallclock
 *                     start/stop time to zero.
 *
 * For user time:      If we have access to timer utilities in sys/times.h,
 *                     we use the times() utility to compute user and
 *                     system start/stop point (passing in the tms struct).
 *                     If we don't have these utilities, we simply set the
 *                     user and start/stop times to zero.
 *
 * The timer manager manipulates the exclusive time information;  i.e.
 * when the timer's startExclusive and stopExclusive routines are called.
 *
 ***************************************************************************
 */
#ifndef ENABLE_SAMRAI_TIMERS
#ifdef __INTEL_COMPILER
#pragma warning (disable:869)
#pragma warning (disable:1101)
#endif // __INTEL_COMPILER
#endif // ENABLE_SAMRAI_TIMERS

SAMRAI_INLINE_KEYWORD void Timer::startExclusive()
{
#ifdef ENABLE_SAMRAI_TIMERS
   if (d_is_active) {

      Clock::timestamp(d_user_start_exclusive,
         d_system_start_exclusive,
         d_wallclock_start_exclusive);

   }
#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD void Timer::stopExclusive()
{
#ifdef ENABLE_SAMRAI_TIMERS
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
#endif // ENABLE_SAMRAI_TIMERS
}

/*
 ***************************************************************************
 *
 * Simple utility routines to manipulate timers.
 *
 ***************************************************************************
 */

SAMRAI_INLINE_KEYWORD const std::string& Timer::getName() const
{
   return d_name;
}

SAMRAI_INLINE_KEYWORD bool Timer::isActive() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_is_active;

#else
   return false;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD bool Timer::isRunning() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_is_running;

#else
   return false;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getTotalUserTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_user_total / Clock::getClockCycle();

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getTotalSystemTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_system_total / Clock::getClockCycle();

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getTotalWallclockTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_wallclock_total;

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getMaxWallclockTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_max_wallclock;

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getExclusiveUserTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_user_exclusive / Clock::getClockCycle();

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getExclusiveSystemTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_system_exclusive / Clock::getClockCycle();

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD double Timer::getExclusiveWallclockTime() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_wallclock_exclusive;

#else
   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD void Timer::addConcurrentTimer(
   const Timer& timer)
{
#ifdef ENABLE_SAMRAI_TIMERS
   if (!isConcurrentTimer(timer)) {
      d_concurrent_timers.push_back(&timer);
   }
#endif // ENABLE_SAMRAI_TIMERS
}

SAMRAI_INLINE_KEYWORD bool Timer::isConcurrentTimer(
   const Timer& timer) const
{
#ifdef ENABLE_SAMRAI_TIMERS
   for (std::vector<const Timer *>::const_iterator i = d_concurrent_timers.begin();
        i != d_concurrent_timers.end();
        ++i) {
      if (*i == &timer) {
         return true;
      }
   }

#endif // ENABLE_SAMRAI_TIMERS
   return false;
}

SAMRAI_INLINE_KEYWORD int Timer::getNumberAccesses() const
{
#ifdef ENABLE_SAMRAI_TIMERS
   return d_accesses;

#else
   return 0;

#endif // ENABLE_SAMRAI_TIMERS
}

/*
 ***************************************************************************
 *
 * Private utility routines to manipulate timers in database.
 *
 ***************************************************************************
 */

SAMRAI_INLINE_KEYWORD void Timer::setActive(
   bool is_active)
{
#ifdef ENABLE_SAMRAI_TIMERS
   d_is_active = is_active;
#endif // ENABLE_SAMRAI_TIMERS
}

}
}
