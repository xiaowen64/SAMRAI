//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/timers/Foo.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Routine to time some routines in the dummy class Foo.
//
#include "Foo.h"
#include "tbox/Pointer.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"

#define LOOP_MAX (10000000)

Foo::Foo()
{
   d_depth = 0;
   d_max_depth = 0;
}

Foo::~Foo() 
{
} 


void Foo::timerOff( )
{
   tbox::Pointer<tbox::Timer> timer = tbox::TimerManager::getManager()->
      getTimer("dummy::SomeClassName::shouldBeTurnedOff");
   timer->start();


   timer->stop();
}

void Foo::timerOn()
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::timerOn()");
   timer->start();

   timer->stop();
}


void Foo::zero(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::zero()");
   if (depth > 0) {
      timer->start();
      one(depth);
      timer->stop();
   }
}

void Foo::one(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::one()");
   if (depth > 1) {
      timer->start();
      two(depth);
      timer->stop();
   }
}

void Foo::two(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::two()");
   if (depth > 2) {
      timer->start();
      three(depth);
      timer->stop();
   }
}

void Foo::three(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::three()");
   if (depth > 3) {
      timer->start();
      four(depth);
      timer->stop();
   }
}

void Foo::four(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::four()");
   if (depth > 4) {
      timer->start();
      five(depth);
      timer->stop();
   }
}

void Foo::five(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::five()");
   if (depth > 5) {
      timer->start();
      six(depth);
      timer->stop();
   }
}

void Foo::six(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::six()");
   if (depth > 6) {
      timer->start();
      seven(depth);
      timer->stop();
   }
}

void Foo::seven(int depth)
{
   tbox::Pointer<tbox::Timer> timer =
      tbox::TimerManager::getManager()->
      getTimer("apps::Foo::seven()");
   if (depth > 7) {
      TBOX_ERROR("Seven levels is maximum implemented in Foo."
		 << "\n please reset exclusive_time_level to something <= 7.");
   }
}

void Foo::startAndStop(string& name) 
{
   tbox::Pointer<tbox::Timer> timer = tbox::TimerManager::getManager()->getTimer(name);
   timer->start();

   timer->stop();
}

void Foo::setMaxDepth(int max_depth)
{
   d_max_depth = max_depth;
}

void Foo::start(string& name) 
{
   d_depth++;

   tbox::Pointer<tbox::Timer> timer = tbox::TimerManager::getManager()->getTimer(name);
   timer->start();

#ifndef LACKS_SSTREAM
   ostringstream osst;
   osst << d_depth;
   string out = osst.str();
   if (d_depth < d_max_depth) {
      start(out);
   }
   
   stop(out);
#endif
}

void Foo::stop(string& name) 
{
   tbox::Pointer<tbox::Timer> timer = tbox::TimerManager::getManager()->getTimer(name);
   timer->stop();
}
