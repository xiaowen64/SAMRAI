//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/timers/Foo.h $
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Simple example to demonstrate input/restart of patch data.
//
 
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include <string>
using namespace std;
#define included_String
#include "tbox/Serializable.h"

using namespace SAMRAI;

class Foo  
{
public:
   Foo();
   ~Foo();
 
   void timerOff();
   void timerOn();
   void zero(int depth);
   void one(int depth);
   void two(int depth);
   void three(int depth);
   void four(int depth);
   void five(int depth);
   void six(int depth);
   void seven(int depth);
   void startAndStop(string& name);
   void setMaxDepth(int max_depth);
   void start(string& name);
   void stop(string& name);

private:

   int d_depth;
   int d_max_depth;
   
};
