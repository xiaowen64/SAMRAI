//
// File:        Foo.h
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
