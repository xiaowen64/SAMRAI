//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/timers/Clock.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Accesses system times. 
//

#include "tbox/Clock.h"

#include <stdlib.h>
#include "tbox/SAMRAI_MPI.h"


namespace SAMRAI {
   namespace tbox {

#ifndef _MSC_VER
struct tms   Clock::s_tms_buffer;
#endif
clock_t      Clock::s_null_clock_t;

/*
*************************************************************************
*                                                                       *
* Initialize clock.                                                     *
*                                                                       *
*************************************************************************
*/
void Clock::initialize(clock_t& clock)
{
#ifndef _MSC_VER
   clock = times(&s_tms_buffer);
#endif
}

void Clock::initialize(double& clock)
{
   clock = 0.0;
}

/*
*************************************************************************
*                                                                       *
* Timestamp the provided structures with current system clock readings. *
*                                                                       *
*************************************************************************
*/
void Clock::timestamp(clock_t& user, clock_t& sys, clock_t& wall)
{
#ifndef _MSC_VER
   wall = times(&s_tms_buffer);
   sys  = s_tms_buffer.tms_stime;
   user = s_tms_buffer.tms_utime;
#endif
}

void Clock::timestamp(clock_t& user, clock_t& sys, double& wall)
{
#ifndef _MSC_VER
   s_null_clock_t = times(&s_tms_buffer);
#ifdef HAVE_MPI
   wall = MPI_Wtime();
#else
   wall = 0.0;
#endif
   sys  = s_tms_buffer.tms_stime;
   user = s_tms_buffer.tms_utime;
#endif
}

/*
*************************************************************************
*                                                                       *
* Get the clock cycle used by the system (time is then computed         *
* as measured_time/clock_cycle)                                         *
*                                                                       *
*************************************************************************
*/
double Clock::getClockCycle()
{
#ifdef _MSC_VER
   double clock_cycle = 1.0;
#else
   double clock_cycle = double(sysconf(_SC_CLK_TCK));
#endif
   return(clock_cycle);
}


}
}


