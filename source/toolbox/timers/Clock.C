//
// File:        Clock.C
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Accesses system times. 
//

#include "tbox/Clock.h"

#include <stdlib.h>
#include "tbox/MPI.h"


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
   clock = 0.;
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
   wall = 0.;
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
   double clock_cycle = 1;
#else
   double clock_cycle = double(sysconf(_SC_CLK_TCK));
#endif
   return(clock_cycle);
}


}
}


