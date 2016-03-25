//
// File:        $URL$
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Test program to demonstrate/test logging useing default appenders.
//

#include "SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "tbox/SAMRAI_MPI.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/Utilities.h"
#include "tbox/Logger.h"

#include <string>
#include <iostream>
using namespace std;

using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("defaultloggertest.log");      

      /* 
       * Write a test warning message.
       */
      TBOX_WARNING("Test warning");

      /* 
       * Write a test debug message. Shouldn't see this since
       * Debug messages are off by default.
       */
      TBOX_DEBUG("Test debug1 : should not show up");
      
      tbox::Logger::getInstance() -> setDebug(true);

      /* 
       * Write a test debug message. Should see this 
       * one since we have turned on debug messages.
       */
      TBOX_DEBUG("Test debug2 : should show up");
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
   return(fail_count);
}


