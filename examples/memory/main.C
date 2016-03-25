//
// File:        main.C
// Package:     SAMRAI application
// Copyright:   (c) 1997-2002 The Regents of the University of California
// Revision:    $Revision: 7 $
// Modified:    $Date: 2004-11-30 13:18:17 -0800 (Tue, 30 Nov 2004) $
// Description: Example program to demonstrate timers.
//

#include "SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "tbox/SAMRAIManager.h"
#include "tbox/MemoryUtilities.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"

using namespace SAMRAI;

#define MAX_TEST 500
#define BYTES_DOUBLE 8 

int main( int argc, char *argv[] )
{
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   tbox::PIO::logAllNodes("Timer.log");      

#ifdef HAVE_TAU
//   TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);
//   TAU_PROFILE_SET_NODE(tbox::MPI::getRank());
#endif

   tbox::pout << "\n\nAllocating memory in 1MB chunks until we run out...\n" << endl;

   for (int chunk = 1; chunk < MAX_TEST; chunk++) {
      
      int doubles_in_one_mb = 1024*1024/BYTES_DOUBLE;
 
      double* array = new double[chunk*doubles_in_one_mb];
      if (!array) {
         tbox::pout << "\nRan out of memory!!" << endl;
         break;
      }

      /* 
       * Touch entries to assure they are really allocated.
       */
      for (int i = 0; i < chunk*doubles_in_one_mb; i++) {
         array[i] = 0.;
      }
      
         
      tbox::pout << "Successfully allocated " << chunk << "MB chunk (" 
           << doubles_in_one_mb*chunk << " doubles)"
           << endl;

      tbox::MemoryUtilities::printMemoryInfo(tbox::pout);
      tbox::MemoryUtilities::recordMemoryInfo();

      delete array;
   }
      

   /*
    * We're done.  Shut down application ...
    */
   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();
   return(0);
}


