//
// File:        Timer.C
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Timer class to track elapsed time in portions of a program. 
//

#include "tbox/Timer.h"

#include "tbox/MPI.h"
#include "tbox/IOStream.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"


#ifdef HAVE_VAMPIR
extern "C" {
#include "VT.h"
}
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#define TBOX_TIMER_VERSION (1)

#ifdef DEBUG_NO_INLINE
#include "tbox/Timer.I"
#endif

namespace SAMRAI {
   namespace tbox {

#ifndef TIMER_MAX_CONCURRENT_TIMERS
#define TIMER_MAX_CONCURRENT_TIMERS (128)
#endif

/*
*************************************************************************
*                                                                       *
* The constructor sets the timer name and initializes timer state.      *
*                                                                       *
*************************************************************************
*/

Timer::Timer(const string& name,
                       const int id)
{
   d_name = name;
   d_identifier = id;
   d_is_running = false; 
   d_is_active = true; 
   d_accesses = 0;
   d_concurrent_timers = new bool[TIMER_MAX_CONCURRENT_TIMERS];

#ifdef HAVE_VAMPIR
   string::size_type position;
 
   // parse timers name down to method
   position = name.find("::");
   string class_method = name.substr(position+2);
   position = class_method.find("::");
   string class_name = class_method.substr(0,position);
   string method = class_method.substr(position+2);
 
   // convert strings to char* type
   char* char_method = new char[method.length()+1];
   method.copy(char_method,string::npos);
   char_method[method.length()] = 0;
 
   char* char_class = new char[class_name.length()+1];
   class_name.copy(char_class,string::npos);
   char_class[class_name.length()] = 0;
 
   VT_symdef(id,char_class,char_method);
#endif

#ifdef HAVE_TAU
   /*
    * Create a Tau "timer" to track time.  
    */
   TAU_MAPPING_TIMER_CREATE(tautimer, name, " ", 
                           TAU_USER2, "SAMRAI_DEFAULT");
#endif

   Clock::initialize(d_user_start_exclusive);
   Clock::initialize(d_user_stop_exclusive);
   Clock::initialize(d_system_start_exclusive);
   Clock::initialize(d_system_stop_exclusive);
   Clock::initialize(d_wallclock_start_exclusive);
   Clock::initialize(d_wallclock_stop_exclusive);

   reset();
}

Timer::~Timer()
{
   delete [] d_concurrent_timers;
}

/*
 ***************************************************************************
 *                                                                         *
 * Start and stop routines for timers.                                     *
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
 * Note that the stop routine increments the elapsed time information.     *
 * Also, the timer manager manipulates the exclusive time information      *
 * the timers when start and stop are called.                              *
 *                                                                         *
 ***************************************************************************
 */

void Timer::start()
{
   d_accesses++;

   if (d_is_active) {

      Clock::timestamp(d_user_start_total, 
                            d_system_start_total, 
                            d_wallclock_start_total);

      d_is_running = true;

#ifdef HAVE_VAMPIR
      VT_begin(d_identifier);
#endif

#ifdef HAVE_TAU
      /*
       * Start the TAU timer.  The "tid" is used for threaded systems
       * so it generally won't apply for us.  The profiler accesses
       * the timer (given the "tautimer" for this timer object) which
       * returns "t", and then starts "t".
       */
      int tid = RtsLayer::myThread();
      TAU_MAPPING_PROFILE_TIMER(t, tautimer, tid);
      TAU_MAPPING_PROFILE_START(t, tid);
#endif

      TimerManager::getManager()->startTime(this);

   }
}

void Timer::stop()
{
   if (d_is_active) {

      TimerManager::getManager()->stopTime(this);

#ifdef HAVE_VAMPIR
      VT_end(d_identifier);
#endif

#ifdef HAVE_TAU
      TAU_MAPPING_PROFILE_STOP(RtsLayer::myThread());
#endif

      Clock::timestamp(d_user_stop_total, 
                            d_system_stop_total, 
                            d_wallclock_stop_total);

      d_is_running = false;

      d_wallclock_total +=
         double(d_wallclock_stop_total - d_wallclock_start_total);
      d_user_total += double(d_user_stop_total - d_user_start_total);
      d_system_total += double(d_system_stop_total - d_system_start_total);

   }

}

void Timer::reset()
{
   d_user_total = 0.;
   d_system_total = 0.;
   d_wallclock_total = 0.;

   d_user_exclusive = 0.;
   d_system_exclusive = 0.;
   d_wallclock_exclusive = 0.;

   d_max_wallclock = 0.;

   for (int i = 0; i < TIMER_MAX_CONCURRENT_TIMERS; i++) {
      d_concurrent_timers[i] = false;
   }
}

/*
 ***************************************************************************
 *                                                                         *
 * Compute the load balance efficiency based the wallclock time on each    *
 * processor, using the formula:                                           *
 *                                                                         *
 *      eff = (sum(time summed across processors)/#processors) /           *
 *             max(time across all processors)                             *
 *                                                                         *
 * This formula corresponds to that used to compute load balance           *
 * efficiency based on the processor distribution of the the number of     *
 * cells (i.e. in BalanceUtilities::computeLoadBalanceEfficiency).         *
 *                                                                         *
 ***************************************************************************
 */
double Timer::computeLoadBalanceEfficiency() 
{
   double wall_time = d_wallclock_total;
   double sum = MPI::sumReduction(wall_time);
   computeMaxWallclock();
   int nprocs = MPI::getNodes();
   double eff = 100.;
   if (d_max_wallclock > 0.) {
      eff = 100.*(sum/(double)nprocs)/d_max_wallclock;
   }
   return eff;
}

void Timer::computeMaxWallclock()  
{
   double wall_time = d_wallclock_total;
   d_max_wallclock = MPI::maxReduction(wall_time);   
}

void Timer::putToDatabase(
   Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif
   db->putInteger("TBOX_TIMER_VERSION",
                   TBOX_TIMER_VERSION);

   db->putString("d_name",d_name);

   db->putDouble("d_user_total",d_user_total);
   db->putDouble("d_system_total", d_system_total);
   db->putDouble("d_wallclock_total", d_wallclock_total);

   db->putDouble("d_user_exclusive",d_user_exclusive);
   db->putDouble("d_system_exclusive", d_system_exclusive);
   db->putDouble("d_wallclock_exclusive", d_wallclock_exclusive);
}

void Timer::getFromRestart(
   Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif
   int ver = db->getInteger("TBOX_TIMER_VERSION");
   if (ver != TBOX_TIMER_VERSION) {
      TBOX_ERROR("Restart file version different than class version.");
   }

   d_name = db->getString("d_name");

   d_user_total = db->getDouble("d_user_total");
   d_system_total = db->getDouble("d_system_total");
   d_wallclock_total = db->getDouble("d_wallclock_total");

   d_user_exclusive = db->getDouble("d_user_exclusive");
   d_system_exclusive = db->getDouble("d_system_exclusive");
   d_wallclock_exclusive = db->getDouble("d_wallclock_exclusive");
}



}
}
