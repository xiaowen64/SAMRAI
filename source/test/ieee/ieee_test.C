//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/ieee/ieee_test.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	simple test code for the IEEE exception handlers
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
using namespace std;

#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/IEEE.h"

using namespace SAMRAI;

//  A successful signaling operation will stop the program, making it
//  difficult to test a variety of operations in a single executable.  
//  Uncomment the particular case you want to test and recompile:

//#define SIGOPS_SUN_F0
//#define SIGOPS_SUN_D0
//#define SIGOPS_SUN_F1
//#define SIGOPS_SUN_D1
//#define SIGOPS_F2
//#define SIGOPS_D2
//#define SIGOPS_F3
//#define SIGOPS_D3



int main(int argc, char **argv)
{
   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   tbox::IEEE::setupFloatingPointExceptionHandlers();

   /*
    * Try using variables set to NaN and see if the signal handler catches
    * and stops the program.  Note that the TESTSIGNAL flag is not defined
    * by default as it will stop this executable, signaling to the script
    * there was an error, even though this is the desired behavior.
    */

   int fail_count = 0;

#ifdef SIGOPS_SUN_F0
   float f0 = tbox::MathUtilities<float>::getSignalingNaN();
   tbox::pout << "\nThis operation should stop the executable if signal handling " 
        << "\nis working (under solaris)..." << endl;
   tbox::pout << "f0 = " << f0 << endl;
   tbox::pout << "f0 * 1.0 = " << endl;
   tbox::pout << f0 * 1.0 << endl;
   tbox::pout << "The operation worked so SIGNAL HANDLING IS NOT ENABLED" << endl;
#endif

#ifdef SIGOPS_SUN_D0
   double d0 = tbox::MathUtilities<double>::getSignalingNaN();
   tbox::pout << "\nThis operation should stop the executable if signal handling " 
        << "\nis working (under solaris)..." << endl;
   tbox::pout << "d0 = " << d0 << endl;
   tbox::pout << "d0 * 1.0 = " << endl;
   tbox::pout << d0 * 1.0 << endl;
   tbox::pout << "The operation worked so SIGNAL HANDLING IS NOT ENABLED" << endl;
#endif

#ifdef SIGOPS_SUN_F1
   float f1 = tbox::MathUtilities<float>::getSignalingNaN();
   tbox::pout << "\nThis operation should stop the executable if signal handling " 
        << "\nis working (under solaris)..." << endl;
   tbox::pout << "f1 = " << f1 << endl;
   tbox::pout << "f1 * 1.0 = " << endl;
   tbox::pout << f1 * 1.0 << endl;
   tbox::pout << "The operation worked so SIGNAL HANDLING IS NOT ENABLED" << endl;
#endif

#ifdef SIGOPS_SUN_D1
   double d1 = tbox::MathUtilities<double>::getSignalingNaN();
   tbox::pout << "\nThis operation should stop the executable if signal handling " 
        << "\nis working (under solaris)..." << endl;
   tbox::pout << "d1 = " << d1 << endl;
   tbox::pout << "d1 * 1.0 = " << endl;
   tbox::pout << d1 * 1.0 << endl;
   tbox::pout << "The operation worked so SIGNAL HANDLING IS NOT ENABLED" << endl;
#endif

#ifdef SIGOPS_F2
   float f2 = 1.0;
   tbox::pout << "\nThis operation should stop the executable if signal handling " 
        << "is working..." << endl;
   tbox::pout << "f2 = " << f2 << endl;
   tbox::pout << "f2/0.0 = " << endl;
   tbox::pout << f2/0.0 << endl;
   tbox::pout << "The operation worked so SIGNAL HANDLING IS NOT ENABLED" << endl;
#endif

#ifdef SIGOPS_D2
   double d2 = 1.0;
   tbox::pout << "\nThis operation should stop the executable if signal handling " 
        << "is working..." << endl;
   tbox::pout << "d2 = " << d2 << endl;
   tbox::pout << "d2/0.0 = " << endl;
   tbox::pout << d2/0.0 << endl;
   tbox::pout << "The operation worked so SIGNAL HANDLING IS NOT ENABLED" << endl;
#endif

   /*
    * Check if the isNaN() method is working correctly
    */
#ifdef SIGOPS_F0
   bool is_f0_nan = tbox::MathUtilities<float>::isNaN(f0);
   if (!is_f0_nan) {
     fail_count++;
     tbox::perr << "Test f0 FAILED" << endl;
   } else {
     tbox::plog << "Test f0 passed" << endl;
   }
#endif

#ifdef SIGOPS_D0
   bool is_d0_nan = tbox::MathUtilities<float>::isNaN(d0);
   if (!is_d0_nan) {
     fail_count++;
     tbox::perr << "Test d0 FAILED" << endl;
   } else {
     tbox::plog << "Test d0 passed" << endl;
   }
#endif

#ifdef SIGOPS_F1
   bool is_f1_nan = tbox::IEEE::isNaN(f1);
   if (!is_f1_nan) {
     fail_count++;
     tbox::perr << "Test f1 FAILED" << endl;
   } else {
     tbox::plog << "Test f1 passed" << endl;
   }
#endif

#ifdef SIGOPS_D1
   bool is_d1_nan = tbox::IEEE::isNaN(d1);
   if (!is_d1_nan) {
     fail_count++;
     tbox::perr << "Test d1 FAILED" << endl;
   } else {
     tbox::plog << "Test d1 passed" << endl;
   }
#endif

#ifdef SIGOPS_F2
   bool is_f2_nan = tbox::IEEE::isNaN(f2);
   if (is_f2_nan) {
     fail_count++;
     tbox::perr << "Test f2 FAILED" << endl;
   } else {
     tbox::plog << "Test f2 passed" << endl;
   }
#endif

#ifdef SIGOPS_D2
   bool is_d2_nan = tbox::IEEE::isNaN(d2);
   if (is_d2_nan) {
     fail_count++;
     tbox::perr << "Test d2 FAILED" << endl;
   } else {
     tbox::plog << "Test d2 passed" << endl;
   }
#endif

#ifdef SIGOPS_F3
   //  Try some arrays
   tbox::Array<float> f3(5);
   tbox::IEEE::initializeArrayToSignalingNaN(f3);
   for (int i = 0; i < f3.getSize(); i++) {
      bool is_f3_nan = tbox::IEEE::isNaN(f3[i]);
      if (!is_f3_nan) {
         fail_count++;
         tbox::perr << "Test f3 FAILED; i = " << i << endl;
      } else {
         tbox::plog << "Test f3 passed" << endl;
      }
   }
#endif

#ifdef SIGOPS_D3
   tbox::Array<double> d3(5);
   tbox::IEEE::initializeArrayToSignalingNaN(d3);
   for (int i = 0; i < d3.getSize(); i++) {
      bool is_d3_nan = tbox::IEEE::isNaN(d3[i]);
      if (!is_d3_nan) {
         fail_count++;
         tbox::perr << "Test d3 FAILED; i = " << i << endl;
      } else {
         tbox::plog << "Test d3 passed" << endl;
      }
    }
#endif

   if ( fail_count == 0 ) {
      tbox::pout << "\nPASSED:  ieee" << endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return(fail_count);
}
