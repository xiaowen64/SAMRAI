//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/timers/Foo.h $
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Simple example to demonstrate input/restart of patch data.
//
 
#include "SAMRAI_config.h"

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
