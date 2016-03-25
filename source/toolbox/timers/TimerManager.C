//
// File:        TimerManager.C
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 477 $
// Modified:    $Date: 2005-07-07 15:33:13 -0700 (Thu, 07 Jul 2005) $
// Description: Class to manage different timer objects used throughout the 
//              library.
//

#include "tbox/TimerManager.h"

#include "tbox/InputDatabase.h"
#include "tbox/MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/IOStream.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
   namespace tbox {


#ifndef NULL
#define NULL (0)
#endif
#ifndef TIMER_MANAGER_MAX_TIMERS
#define TIMER_MANAGER_MAX_TIMERS (128)
#endif
#ifndef TIMER_UNREG_ACCESS_TIME
#define TIMER_UNREG_ACCESS_TIME (3.8e-7)
#endif
#ifndef TIMER_REG_ACCESS_TIME
#define TIMER_REG_ACCESS_TIME (1.0e-5)
#endif
#ifndef MAIN_TIMER_IDENTIFIER
#define MAIN_TIMER_IDENTIFIER  (-1)
#endif
#ifndef INACTIVE_TIMER_IDENTIFIER
#define INACTIVE_TIMER_IDENTIFIER  (-9999)
#endif

TimerManager* TimerManager::s_timer_manager_instance = 
                   (TimerManager*)NULL;
bool TimerManager::s_registered_callback = false;

/*
*************************************************************************
*                                                                       *
* Static timer manager member functions.                                *
*                                                                       *
*************************************************************************
*/

void TimerManager::createManager(
   Pointer<Database> input_db)
{
   if (!s_timer_manager_instance) {
      s_timer_manager_instance = new TimerManager(input_db);
   }
   if (!s_registered_callback) {
      ShutdownRegistry::registerShutdownRoutine(freeManager,
		     ShutdownRegistry::priorityTimerManger);

      s_registered_callback = true;
   }

   s_timer_manager_instance->d_main_timer->start();
   
}

TimerManager* TimerManager::getManager()
{
   if (!s_timer_manager_instance) {
      createManager(NULL);
   }

   return(s_timer_manager_instance);
}

void TimerManager::freeManager()
{
   if (s_timer_manager_instance) delete s_timer_manager_instance;
   s_timer_manager_instance = ((TimerManager*) NULL);

}

void TimerManager::registerSingletonSubclassInstance(
   TimerManager* subclass_instance)
{
   if (!s_timer_manager_instance) {
      s_timer_manager_instance = subclass_instance;
      if (!s_registered_callback) {
         ShutdownRegistry::registerShutdownRoutine(freeManager,
		     ShutdownRegistry::priorityTimerManger);
         s_registered_callback = true;
      }
   } else {
      TBOX_ERROR("TimerManager internal error...\n"
                 << "Attemptng to set Singleton instance to subclass instance,"
                 << "\n but Singleton instance already set." << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Protected timer manager constructor and destructor.                   *
*                                                                       *
*************************************************************************
*/

TimerManager::TimerManager(Pointer<Database> input_db)
{

   /*
    * Create a timer that measures overall solution time.  If the 
    * application uses Tau, this timer will effectively measure 
    * uninstrumented parts of the library.  Hence, use a different name 
    * for the different cases to avoid confusion in the Tau analysis tool.
    */
#ifdef HAVE_TAU
   d_main_timer = new Timer("UNINSTRUMENTED PARTS", 
                                 MAIN_TIMER_IDENTIFIER);
#else 
   d_main_timer = new Timer("TOTAL RUN TIME", 
                                 MAIN_TIMER_IDENTIFIER);
#endif

   d_timers.resizeArray(TIMER_MANAGER_MAX_TIMERS); 
   d_num_timers = 0;

   d_inactive_timers.resizeArray(TIMER_MANAGER_MAX_TIMERS);
   d_num_inactive_timers = 0;

   d_running_timers = new bool[TIMER_MANAGER_MAX_TIMERS];
   for (int i = 0; i < TIMER_MANAGER_MAX_TIMERS; i++) {
      d_running_timers[i] = false;
   }

   d_exclusive_timer_stack.clearItems();

   d_length_package_names = 0;
   d_length_class_names = 0;
   d_length_class_method_names = 0;
  
   d_print_exclusive      = false;
   d_print_total          = true;
   d_print_processor      = true;
   d_print_max            = false;
   d_print_summed         = false;
   d_print_user           = false;
   d_print_sys            = false;
   d_print_wall           = true;
   d_print_percentage     = true;
   d_print_concurrent     = false;
   d_print_timer_overhead = false;
  
   d_print_threshold = 0.25; 

   getFromInput(input_db);


}


TimerManager::~TimerManager()
{

   d_main_timer->stop();
   d_main_timer.setNull();

   for (int i = 0; i < d_num_timers; i++) {
      d_timers[i].setNull();
   }
   d_timers.setNull();
   d_num_timers = 0;

   for (int j = 0; j < d_num_inactive_timers; j++) {
      d_inactive_timers.setNull();
   }
   d_inactive_timers.setNull();
   d_num_inactive_timers = 0;

   delete [] d_running_timers;

   d_exclusive_timer_stack.clearItems();

   d_package_names.clearItems();
   d_class_names.clearItems();
   d_class_method_names.clearItems();
}

/*
*************************************************************************
*                                                                       *
* Utility functions for creating timers, adding them to the manager     *
* database, and checking whether a particular timer exists.             *
*                                                                       *
*    o checkTimerExistsInArray is private.  It returns true if a timer  *
*      matching the name string exists in the array and returns false   *
*      otherwise.  If such a timer exists, the pointer is set to that   *
*      timer; otherwise the pointer is null.                            *
*                                                                       *
*    o getTimer returns a timer with the given name.  It will be active *
*      if its name appears in the input file, or if it is added with    *
*      a `true' argument.  Otherwise the timer will be inactive.        *
*                                                                       *
*    o checkTimerExists returns true if a timer matching the name       *
*      string exists and false otherwise.  If such a timer exists,      *
*      the pointer is set to that timer; otherwise the pointer is null. *
*                                                                       *
*************************************************************************
*/

bool TimerManager::checkTimerExistsInArray(
   Pointer<Timer>& timer,
   const string& name,
   const Array< Pointer<Timer> >& timer_array,
   int array_size) const
{
   timer.setNull();

   bool timer_found = false;

   if (!name.empty()) {
      for (int i = 0; i < array_size; i++) {
         if (timer_array[i]->getName() == name) {
            timer_found = true;
            timer = timer_array[i];
            break;
         }
      }
   }

   return(timer_found);
}

Pointer<Timer> TimerManager::getTimer(
   const string& name,
   bool ignore_timer_input)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif

   bool timer_active = true;
   if (!ignore_timer_input) {   
      /*
       * Check if name is in either the d_package_names, d_class_names,
       * or d_class_method_names lists.  Add it if it is.
       */
      timer_active = checkTimerInNameLists(name);
   }

   /*
    * Add the timer to the appropriate array, if necessary.
    */
   Pointer<Timer> timer;
   if (timer_active) {
      if ( !checkTimerExistsInArray(timer,
                                    name,
                                    d_timers,
                                    d_num_timers) ) {
         if (d_num_timers == TIMER_MANAGER_MAX_TIMERS) {
            TBOX_ERROR("TimerManager::getTimer error ..."
               << "\n   Max timers exceeded with timer " << name << endl);
         }
         timer = new Timer(name, d_num_timers);
         d_timers[d_num_timers] = timer;
         d_num_timers++;
      }
   } else {
       if ( !checkTimerExistsInArray(timer,
                                     name,
                                     d_inactive_timers,
                                     d_num_inactive_timers) ) {
          if (d_num_inactive_timers == TIMER_MANAGER_MAX_TIMERS) {
             TBOX_ERROR("TimerManager::getTimer error ..."
                << "\n   Max timers exceeded with timer " << name << endl);
          }
          timer = new Timer(name, INACTIVE_TIMER_IDENTIFIER);
          timer->setInactive(); 
          d_inactive_timers[d_num_inactive_timers] = timer;
	  d_num_inactive_timers++;
       }
   }
   return(timer);
}

bool TimerManager::checkTimerExists(
   Pointer<Timer>& timer,
   const string& name) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif

   bool timer_found = checkTimerExistsInArray(timer,
                                              name,
                                              d_timers,
                                              d_num_timers);

   if (!timer_found) {
      timer_found = checkTimerExistsInArray(timer,
                                            name,
                                            d_inactive_timers,
                                            d_num_inactive_timers);
   }

   return(timer_found);
}

/*
*************************************************************************
*                                                                       *
* Utility functions to check whether timer is running and to reset      *
* all active timers.                                                    *
*                                                                       *
*************************************************************************
*/

bool TimerManager::checkTimerRunning(
   const string& name) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif

   bool is_running = false;

   Pointer<Timer> timer;
   if (checkTimerExistsInArray(timer, name, d_timers, d_num_timers)) {
      is_running = timer->isRunning();
   }

   return(is_running);
}


void TimerManager::resetAllTimers()
{
   d_main_timer->stop();
   d_main_timer->reset();
   d_main_timer->start();

   for (int i = 0; i < d_num_timers; i++) {
      d_timers[i]->reset();
   }
   for (int j = 0; j < d_num_inactive_timers; j++) {
      d_inactive_timers[j]->reset();
   }
}

/*
*************************************************************************
*                                                                       *
* Protected start and stop routines for exclusive timer management.     *
*                                                                       *
*************************************************************************
*/

void TimerManager::startTime(Timer* timer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(timer == (Timer*)NULL));
#endif
   int id = timer->getIdentifier();
   if (id >= 0) d_running_timers[id] = true;

   if (d_print_exclusive) {
      if (!d_exclusive_timer_stack.isEmpty()) {
         ((Timer*) d_exclusive_timer_stack.getFirstItem())->
            stopExclusive();
      }
      Timer* stack_timer = timer;
      d_exclusive_timer_stack.addItem(stack_timer);
      stack_timer->startExclusive();
   }

   if (d_print_concurrent) {
      for (int i = 0; i < d_num_timers; i++) {
         if ( (id != i) && d_running_timers[i] ) {
            d_timers[i]->setConcurrentTimer(id); 
         } 
      }
   }

}

void TimerManager::stopTime(Timer* timer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!(timer == (Timer*)NULL));
#endif
   int id = timer->getIdentifier();
   if (id >= 0) d_running_timers[id] = false;

   if (d_print_exclusive) {
      timer->stopExclusive();
      if (!d_exclusive_timer_stack.isEmpty()) {
         d_exclusive_timer_stack.removeFirstItem();
         if (!d_exclusive_timer_stack.isEmpty()) {
            ((Timer*) d_exclusive_timer_stack.getFirstItem())->
               startExclusive();
         }
      }
   }
}


/*
*************************************************************************
*                                                                       *
* Parser to see if Timer has been registered or not.                    *
*                                                                       *
*************************************************************************
*/


bool TimerManager::checkTimerInNameLists(
   const string& copy) 
{
   string name = copy;
  /*
   * string::size_type is generally an int, but it may depend on vendor's
   * implementation of string class.  Just to be safe, we use the definition
   * from the string class (string::size_type).
   */
   string::size_type string_length, list_entry_length, position;

  /*
   * Specification of whether we will use the timer after comparing to
   * package, class, and class::method lists.
   */
   bool will_use_timer = false;

  /*
   * Start by evaluating the Timer's package.  If the list of packages
   * has length zero, we can skip this step. We can also skip this step 
   * if the timer does not have two "::" in its name.  
   *
   * The manager will assume the following with regard to the specified timer
   * name:  The name has...
   *  1) Two "::" - it is of the form Package::Class::method.  A user can 
   *     turn on the timer by specifying the timers Package, Class, or 
   *     Class::method combination in the input file.  This is the most
   *     versatile. 
   *  2) One "::" - it has the form Class::method.  The timer is assumed not
   *     to be in a package.  A user can turn on the timer by specifying its
   *     Class or class::method combination.
   *  3) No "::" - it has the form Class.  The timer is can be turned on by
   *     entering the class name in the input file.  
   */
   if ( d_length_package_names > 0 ) {
      /*
       * See how many instances of "::" there are in the name.  If there
       * are at least two, the timer has a package which might be in the
       * package list.  
       */
      int occurrences = 0; 
      position = name.find("::");
      if (position < name.size()) {
         occurrences++;
         string substring = name.substr(position+2);
         position = substring.find("::");
         if (position < substring.size()) {
            occurrences++;
         }
      }

      if ( occurrences >= 2 ) {
        /*
         * The timer may be in the package list.  Parse to get its package
         * name and compare to the list entries.
         */
         position = name.find("::");
         string package = name.substr(0,position);

         /*
          * Iterate through package list and see if the timer's package
          * name is there. 
          */
         bool package_exists = false;
         string_length = package.size();
         for (List<string>::Iterator i(d_package_names); i; i++) {
            list_entry_length = i().size();
            if (string_length == list_entry_length) {
               package_exists = ( i() == package );
            }
            if (package_exists) break;
         }
         will_use_timer = package_exists;
      }
   }

   if (!will_use_timer) {

      /*
       * The timer's package has already been compared to the package list, 
       * so we can parse it off to ease further evaluations of the name.  
       */
      int occurrences = 0; 
      position = name.find("::");
      if (position < name.size()) {
         occurrences++;
         string substring = name.substr(position+2);
         position = substring.find("::");
         if (position < substring.size()) {
            occurrences++;
         }
      }
      if ( occurrences >= 2 ) {
         position = name.find("::");
         name = name.substr(position+2);
      }

      /*
       * See if Timer's class is in d_class_names.  If the list of classes
       * has length zero, we can skip this step.  
       */
      if (d_length_class_names > 0) {

         /*
          * If a "::" exists in the timer name, parse off the part before
          * before it.  This is the class name.  If no "::" exists, assume
          * the timer name is a class name (don't do any parsing) and compare
          * it directly to the class list.
          */
         position = name.find("::");
         string class_name;
         if (position < name.size()) {
            class_name = name.substr(0,position);
         } else {
            class_name = name;
         }

         /*
          * Is class name dimensional?  
          */
         string_length = class_name.size();
         string dim = class_name.substr(string_length-1,1);
         bool is_dimensional = false;
         string nondim_class_name;
         if ( dim == "1" || dim == "2" || dim == "3") { 
            is_dimensional = true; 
            nondim_class_name = class_name.substr(0,string_length-1);
         } else {
            nondim_class_name = class_name;
         }

         /*
          * See if class name is in class list.  Accelerate this by comparing
          * the length of the entries.  Do non-dimensional comparison first. 
          */
         string_length = nondim_class_name.size();
         bool class_exists = false;
         for (List<string>::Iterator i(d_class_names); i; i++) {
            list_entry_length = i().size();
            if (string_length == list_entry_length) {
               class_exists = ( i() == nondim_class_name );
            }
            if (class_exists) break;
         }

         /*
          * Now do dimensional comparison.
          */
         string_length = class_name.size();
         if (is_dimensional && !class_exists) {
            for (List<string>::Iterator i(d_class_names); i; i++) {
               list_entry_length = i().size();
               if (string_length == list_entry_length) {
                  class_exists = ( i() == class_name );
               }
               if (class_exists) break;
            }
         }
         will_use_timer = class_exists;

      }

      
      /*
       * See if Timer's class::method name is in d_class_method_names list.  
       *
       * If the list of class_method_names has length zero, we can skip 
       * this step.  Also, if no "::" exists in the timer name then it 
       * cannot be in the timer list so lets avoid the comparison.
       */
      position = name.find("::");
      occurrences = 0;
      if (position < name.size()) {
         occurrences = 1;
      }

      if (!will_use_timer && d_length_class_method_names > 0 && 
          occurrences > 0) {

         /*
          * Parse name before "::" - this is the class name. 
          */
         position = name.find("::");
         string class_name = name.substr(0,position);

         /*
          * Is class name dimensional?  
          */
         string_length = class_name.size();
         string dim = class_name.substr(string_length-1,1);
         bool is_dimensional = false;
         string nondim_name;
         if ( dim == "1" || dim == "2" || dim == "3") { 
            is_dimensional = true; 
            string nondim_class_name = class_name.substr(0,string_length-1);
            string method_name = name.substr(position);
            nondim_name = nondim_class_name;
            nondim_name += method_name;

         } else {
            nondim_name = name;
         }

         /*
          * See if name is in class_method_names list.  Accelerate this by 
          * comparing the length of the entries.  Do non-dimensional 
          * comparison first. 
          */
         bool class_method_exists = false;
         string_length = nondim_name.size();
         for (List<string>::Iterator i(d_class_method_names); i; i++) {
            list_entry_length = i().size();
            if (string_length == list_entry_length) {
               class_method_exists = ( i() == nondim_name );
            }
            if (class_method_exists) break;
         }

         /*
          * Now do dimensional comparison.
          */
         if (is_dimensional && !class_method_exists) {
            string_length = name.size();
            for (List<string>::Iterator i(d_class_method_names); i; i++) {
               list_entry_length = i().size();
               if (string_length == list_entry_length) {
                  class_method_exists = ( i() == name );
               }
               if (class_method_exists) break;
            }
         }

         will_use_timer = class_method_exists;

      }

   }

   return (will_use_timer);

}


/*
*************************************************************************
*                                                                       *
* Print timer information from each processor.                          *
*                                                                       *
*************************************************************************
*/

void TimerManager::print(ostream& os)
{
   /*
    * There are 18 possible timer values that users may wish to look at.
    * (i.e. User/sys/wallclock time, Total or Exclusive, for individual
    * processor, max across all procs, or summed across all procs).
    * This method builds an array that holds these values and outputs
    * them in column format.
    */
   
   /*
    * First, stop the main_timer so we have an accurate measure of 
    * overall run time so far.
    */
   d_main_timer->stop();

   /*
    * If we are doing max or sum operations, make sure timers are 
    * consistent across processors.
    */
   if (d_print_summed || d_print_max) {
     checkConsistencyAcrossProcessors();
   }

   /*
    * Invoke arrays used to hold timer names, timer_values, and 
    * max_processor_ids (i.e. processor holding the maximum time).
    */
   double (*timer_values)[18] = new double[d_num_timers+1][18];
   int (*max_processor_id)[2] = new int[d_num_timers+1][2];
   Array<string> timer_names(d_num_timers+1);

   /*
    * Fill in timer_values and timer_names arrays, based on values of
    * d_print_total, d_print_exclusive,
    * d_print_user, d_print_sys, d_print_wallclock, 
    * d_print_processor, d_print_summed, d_print_max.
    */
   buildTimerArrays(timer_values, max_processor_id, timer_names);


   /* 
    * Now that we have built up the array, select output options.
    *
    * 1) If user has requested any two (or more) of {user, system, and 
    *    walltime} then use the format:
    *    [Exclusive,Total]
    *    [Processor,Summed,Max]
    *       Name     User    System   Wall  (Max Processor)
    *
    * 2) If user chose just one of {User, system, or walltime}, then use the 
    *    format:
    *    [Exclusive,Total]
    *       Name     Processor   Summed   Max  Max Processor
    *
    * 3) If user chose just one of {user, system, or walltime} and just one 
    *    of {processor, summed, max} then use the format:
    *       Name     [Processor,Summed,Max]  Total   Exclusive
    * 
    * If the user wants overhead stats, print those as:
    *   Timer Overhead:
    *     Timer Name    number calls     Estimated overhead
    *
    * If they want output of a concurrent timer tree, print this as:
    *   Concurrent Tree:
    *     Timer Name    names of timers called by it.
    */

 
   /*
    * Determine which case we are doing - #1, #2, or #3
    * #1 case - user requested any two (or more) of [user,system,wallclock]
    * #2 case - user requested one of [user,system,wallclock]
    * #3 case - user requested one of [user,system,wallclock] and one of 
    *           [processor, summed, max]
    */
   bool case1 = false;
   bool case2 = false;
   bool case3 = false;
   if ( (d_print_user && d_print_sys) ||
        (d_print_sys  && d_print_wall) ||
        (d_print_wall && d_print_user) ) {
      case1 = true;
   } else {
      case2 = true;
      if ( (d_print_processor && !d_print_summed && !d_print_max) ||
           (!d_print_processor && d_print_summed && !d_print_max) ||
           (!d_print_processor && !d_print_summed && d_print_max) ) {
         case2 = false;
         case3 = true;
      }
   }

   string table_title;
   Array<string> column_titles(4);
   int column_ids[3] = {0,0,0};
   int j,k;

   /*
    * Now print out the data
    */
   if (case1) {
   
      column_titles[0] = "";
      column_titles[1] = "";
      column_titles[2] = "";
      column_titles[3] = "Proc";

      if (d_print_user) {
         column_titles[0] = "User Time";
      } else {
         column_titles[0] = "";
      }
      if (d_print_sys) {
         column_titles[1] = "Sys Time";
      } else {
         column_titles[1] = "";
      }
      if (d_print_wall) {
         column_titles[2] = "Wall Time";
      } else {
         column_titles[2] = "";
      }

      for (k = 0; k < 2; k++) {

         if ( (k == 0 && d_print_exclusive) ||
              (k == 1 && d_print_total) ) {

            for (j = 0; j < 3; j++) {
             
	       if ( (j == 0 && d_print_processor) ||
                    (j == 1 && d_print_summed) ||
                    (j == 2 && d_print_max) ) {
 
                  if ( j == 0) {
#ifndef LACKS_SSTREAM
                    ostringstream out;
#endif
                    if ( k == 0 ) {
#ifndef LACKS_SSTREAM
                       out << "EXCLUSIVE TIME \nPROCESSOR:" 
                           << MPI::getRank();
                       table_title = out.str();
#else
                       table_title = "EXCLUSIVE TIME \nPROCESSOR:";
#endif
                       column_ids[0] = 0;
                       column_ids[1] = 1;
                       column_ids[2] = 2;
                    } else if ( k == 1 ) {
#ifndef LACKS_SSTREAM
                       out << "TOTAL TIME \nPROCESSOR:" 
                           << MPI::getRank();
                       table_title = out.str();
#else
                       table_title = "TOTAL TIME \nPROCESSOR:";
#endif
                       column_ids[0] = 9;
                       column_ids[1] = 10;
                       column_ids[2] = 11;
                    }
		    printTable(table_title, 
			       column_titles, 
			       timer_names,
			       column_ids,
			       timer_values,
			       os);
		  } else if (j == 1) {

                    if ( k == 0 ) {
                       table_title = 
                         "EXCLUSIVE TIME \nSUMMED ACROSS ALL PROCESSORS";
                       column_ids[0] = 3;
                       column_ids[1] = 4;
                       column_ids[2] = 5;
                    } else if ( k == 1 ) {
                       table_title = 
                         "TOTAL TIME \nSUMMED ACROSS ALL PROCESSORS:";
                       column_ids[0] = 12;
                       column_ids[1] = 13;
                       column_ids[2] = 14;
                    }
		    printTable(table_title, 
			       column_titles, 
			       timer_names,
			       column_ids,
			       timer_values,
			       os);

		  } else if ( j == 2) {

                    int max_array_id = 0; // identifies which of the two
		                          // max_processor_id values to print
                    if ( k == 0 ) {
                       table_title = 
                          "EXCLUSIVE TIME \nMAX ACROSS ALL PROCESSORS";
                       column_ids[0] = 6;
                       column_ids[1] = 7;
                       column_ids[2] = 8;
                       max_array_id = 0;
                    } else if ( k == 1 ) {
                       table_title = 
                          "TOTAL TIME \nMAX ACROSS ALL PROCESSORS";
                       column_ids[0] = 15;
                       column_ids[1] = 16;
                       column_ids[2] = 17;
                       max_array_id = 1;
                    }
#if 0
		    printTable(table_title, 
			       column_titles, 
			       timer_names,
                               &max_processor_id[0][max_array_id],
			       column_ids,
			       timer_values,
			       os);                    
#else
		    printTable(table_title, 
			       column_titles, 
			       timer_names,
                               max_processor_id,
                               max_array_id,
			       column_ids,
			       timer_values,
			       os);                    
#endif
		  }

	       } // if j
	    } // for j
	 } // if k
      } // for k
   } // if case 1

   if (case2) {
         
       for (k = 0; k < 2; k++) {

         if ( (k == 0 && d_print_exclusive) ||
              (k == 1 && d_print_total) ) {

            int max_array_id = 0;
            string table_title_line_1;
            string table_title_line_2;            
            if ( k == 0 ) {
               table_title_line_1 = "EXCLUSVE \n";
               max_array_id = 0;
            } else if ( k == 1 ) {
	       table_title_line_1 = "TOTAL \n";
               max_array_id = 1;
            }
            if (d_print_user) {
	       table_title_line_2 = "USER TIME"; 
            } else if (d_print_sys) {
	       table_title_line_2 = "SYSTEM TIME";
            } else if (d_print_wall) {
	       table_title_line_2 = "WALLCLOCK TIME";
            }
            table_title = table_title_line_1;
            table_title += table_title_line_2;


            column_titles[0] = "";
            column_titles[1] = "";
            column_titles[2] = "";
            column_titles[3] = "";
            if (d_print_processor) {
#ifndef LACKS_SSTREAM
               ostringstream out;
	       out << "Proc: " << MPI::getRank();
               column_titles[0] = out.str();
#else
               column_titles[0] = "Proc: ";
#endif
            } 
            if (d_print_summed) {
               column_titles[1] = "Summed";
            }
            if (d_print_max) {
               column_titles[2] = "Max";
               column_titles[3] = "Proc";
            }

            if (d_print_user) {
              if ( k == 0 ) {
                column_ids[0] = 0;
                column_ids[1] = 3;
                column_ids[2] = 6;
              } else if ( k == 1 ) {
                column_ids[0] = 9;
                column_ids[1] = 12;
                column_ids[2] = 15;
              }
            } else if (d_print_sys) {
              if ( k == 0 ) {
                column_ids[0] = 1;
                column_ids[1] = 4;
                column_ids[2] = 7;
              } else if ( k == 1 ) {
                column_ids[0] = 10;
                column_ids[1] = 13;
                column_ids[2] = 16;
              }
            } else if (d_print_wall) {
              if ( k == 0 ) {
                column_ids[0] = 2;
                column_ids[1] = 5;
                column_ids[2] = 8;
              } else if ( k == 1 ) {
                column_ids[0] = 11;
                column_ids[1] = 14;
                column_ids[2] = 17;
              } 
            }
#if 0
	    printTable(table_title, 
		       column_titles, 
		       timer_names,
                       &max_processor_id[0][max_array_id],
		       column_ids,
		       timer_values,
		       os);     
#else
	    printTable(table_title, 
		       column_titles, 
		       timer_names,
                       max_processor_id,
                       max_array_id,
		       column_ids,
		       timer_values,
		       os);     
#endif
	 } // if k
       } // for k
   } // if case2


   if (case3) {

      if (d_print_exclusive && !d_print_total) {
         column_titles[0] = "Exclusive";
         column_titles[1] = "";
      } else if (!d_print_exclusive && d_print_total) {
         column_titles[0] = "";
         column_titles[1] = "Total";
      } else if (d_print_exclusive && d_print_total) {
         column_titles[0] = "Exclusive";
         column_titles[1] = "Total"; 
      } 
      column_titles[3] = "";
  
      column_ids[2] = 0;
      if ( d_print_user ) {
	if (d_print_processor) {
#ifndef LACKS_SSTREAM
           ostringstream out;
           out << "USER TIME \nPROCESSOR: " << MPI::getRank();
           table_title = out.str();
#else
           table_title = "USER TIME \nPROCESSOR: ";
#endif
           column_ids[0] = 0;
           column_ids[1] = 9;
        } else if (d_print_summed) {
           table_title = "USER TIME \nSUMMED ACROSS ALL PROCESSORS";
           column_ids[0] = 3;
           column_ids[1] = 12;
        } else if (d_print_max) {
           table_title = "USER TIME \nMAX ACROSS ALL PROCESSORS";
           column_ids[0] = 6;
           column_ids[1] = 15;
        }
      } else if ( d_print_sys ) {
	if (d_print_processor) {
#ifndef LACKS_SSTREAM
           ostringstream out;
           out << "SYSTEM TIME \nPROCESSOR: " << MPI::getRank();
           table_title = out.str();
#else
           table_title = "SYSTEM TIME \nPROCESSOR:";
#endif
           column_ids[0] = 1;
           column_ids[1] = 10;
        } else if (d_print_summed) {
           table_title = "SYSTEM TIME \nSUMMED ACROSS ALL PROCESSORS";
           column_ids[0] = 4;
           column_ids[1] = 13;
        } else if (d_print_max) {
           table_title = "SYSTEM TIME \nMAX ACROSS ALL PROCESSORS";
           column_ids[0] = 7;
           column_ids[1] = 16;
	}
      } else if ( d_print_wall ) {
	if (d_print_processor) {
#ifndef LACKS_SSTREAM
           ostringstream out;
           out << "WALLCLOCK TIME \nPROCESSOR: " << MPI::getRank();
           table_title = out.str();
#else
           table_title = "WALLCLOCK TIME \nPROCESSOR: ";
#endif
           column_ids[0] = 2;
           column_ids[1] = 11;
        } else if (d_print_summed) {
           table_title = "WALLCLOCK TIME \nSUMMED ACROSS ALL PROCESSORS";
           column_ids[0] = 5;
           column_ids[1] = 14;
        } else if (d_print_max) {
           table_title = "WALLCLOCK TIME \nMAX ACROSS ALL PROCESSORS";
           column_ids[0] = 8;
           column_ids[1] = 17;
	}
      }
      printTable(table_title, 
	         column_titles, 
	         timer_names,
	         column_ids,
	         timer_values,
	         os); 
   }

   /*
    * Print overhead stats - number of accesses and estimated cost 
    * (estimated cost computed based on the number of accesses and
    * a fixed REG_ACCESS_TIME - defined as a constant; the REG_ACCESS_TIME
    * is taken as an average of access times measured across a few different
    * machines).  
    * Store the number of accesses in max_processor_id[0] and the estimated
    * cost in timer_values[0] and use the printTable method.
    */   
   if (d_print_timer_overhead) {
      printOverhead(timer_names,
                    timer_values,
                    os);
   }

   /*
    * Print tree of concurrent timers.
    */   
   if (d_print_concurrent) {
      printConcurrent(os);
   }


   delete [] timer_values;
   delete [] max_processor_id;
  /*
   * Lastly, restart the main_timer that we stopped at the beginning of
   * this routine
   */
   d_main_timer->start();
}

void TimerManager::printTable(
                   const string& table_title, 
		   const Array<string> column_titles,   
                   const Array<string> timer_names,
		   const int column_ids[],
		   const double timer_values[][18],
		   ostream& os)
{
   string ascii_line1 =  "++++++++++++++++++++++++++++++++++++++++";
   string ascii_line2 =  "++++++++++++++++++++++++++++++++++++++++\n"; 
   string ascii_line  = ascii_line1;
   ascii_line  += ascii_line2;

   /*
    * By default, output in C++ is right justified with the setw() 
    * option e.g. cout << "[" << setw(5) << 1 << "]" will output 
    * [   1].  The line below makes it left justified, so the same line
    * will generate [1   ].  We us left justification because it is
    * more convenient to output columns of tables.
    */
   os.setf(ios::left);

   os << ascii_line << table_title << "\n";

   int n, i;

   /*
    * Determine maximum name length for formatting
    */
   int maxlen = 10;
   for (n = 0; n < d_num_timers+1; n++) { 
      i = int(timer_names[n].size());
      if (i > maxlen) maxlen = i;
   }


   /*
    * Print table header.  If we are only printing the overall program 
    * timer (i.e. d_num_timers = 0) with only d_print_processor, 
    * d_print_total, and d_print_wall options being true (which
    * is the default case if the user doesn't add a "TimerManager" 
    * section to the input file) then don't bother to print header as 
    * it just clutters up the output.  Also, turn off percentages since 
    * this doesn't mean anything with just one timer.
    */
   bool default_case = d_num_timers == 0 && !d_print_exclusive &&
                      !d_print_summed && !d_print_max &&
                      !d_print_user && !d_print_sys;
   if (default_case) {
      d_print_percentage = false;
   } else {
      os << ascii_line
         << setw(maxlen+3) << "Timer Name" << ' ';
      for (i = 0; i < 3; i++) {
         if ( !column_titles[i].empty() ) {
	   os << setw(15) << column_titles[i].c_str() << "  ";
         }
      }
      os << endl;
   }

   /* 
    * Organize timers largest to smallest.  Apply this to the LAST NONZERO
    * column entry for the table by forming an ordering array - ordered_list 
    * - that orders these values.   
    */
   int last_nonzero_column = 0;
   for (i = 0; i < 3; i++) {
      if ( !column_titles[i].empty() ) {
         last_nonzero_column = column_ids[i];
      }
   }
   int *ordered_list = new int[d_num_timers+1];
   buildOrderedList(timer_values,
                    last_nonzero_column, 
                    ordered_list,
                    d_num_timers);

   /*
    * Tack on TOTAL TIME to end of ordered list
    */
   ordered_list[d_num_timers] = d_num_timers;

   /*
    * Now output the rows of the table.
    */ 
   
   for (int k = 0; k < d_num_timers+1; k++) {
      n = ordered_list[k];

      /*
       * Check the print threshold to see if we should print this timer.
       */
      double frac = computePercentageDouble(
                       timer_values[n][last_nonzero_column],
		       timer_values[d_num_timers][last_nonzero_column] );

      if (frac > d_print_threshold) {
	
	 os << setw(maxlen+3) << timer_names[n].c_str() << ' ';

         /*
          * Print column values
          */
	 for (i = 0; i < 3; i++) {

	    /*
	     * Print column values only title is non-null (i.e. not "")
	     */ 
	    if ( !column_titles[i].empty() ) {

	       /*
		* Print percentages if requested. 
		*/ 
	      int j = column_ids[i];

	      if (d_print_percentage) {
#ifndef LACKS_SSTREAM
                 int perc = computePercentageInt(timer_values[n][j],
					      timer_values[d_num_timers][j]);

		 ostringstream out;
		 out << timer_values[n][j] << " (" << perc << "%)";
		 os << setw(15) << out.str().c_str() << "  ";
#else 
                 os << setw(15) << timer_values[n][j] << "  ";
#endif
	      } else {
	         os << setw(15) << timer_values[n][j] << "  ";
	      }


	    } // if title is non-null

	 } // loop over columns

	 os << endl;

      } // if meets d_print_threshold condition

   } // loop over timers

   delete[] ordered_list;

   os << ascii_line  << endl;
   os.setf(ios::right);

}

void TimerManager::printTable(
                   const string& table_title, 
		   const Array<string> column_titles,   
                   const Array<string> timer_names,
                   const int max_processor_id[][2],
                   const int max_array_id,
		   const int column_ids[],
		   const double timer_values[][18],
		   ostream& os)
{
   string ascii_line1 =  "++++++++++++++++++++++++++++++++++++++++"; 
   string ascii_line2 =  "++++++++++++++++++++++++++++++++++++++++\n"; 
   string ascii_line  = ascii_line1;
   ascii_line  += ascii_line2;

   /*
    * Left-justify all output in this method.
    */
   os.setf(ios::left);

   os << ascii_line 
      << table_title << "\n"
      << ascii_line;
   
   int n,i;

   /*
    * Determine maximum name length for formatting
    */   
   int maxlen = 10;
   for (n = 0; n < d_num_timers+1; n++) { 
      i = int(timer_names[n].size());
      if (i > maxlen) maxlen = i;
   }

   /*
    * Print table header
    */
   os << setw(maxlen+3) << "Timer Name" << ' ';
   for (i = 0; i < 4; i++) {
      if ( !column_titles[i].empty() ) {
	os << setw(15) << column_titles[i].c_str() << "  ";
      }
   }
   os << endl;

   /* 
    * Organize timers largest to smallest.  Apply this to the LAST NONZERO
    * column entry for the table by forming an ordering array - ordered_list 
    * - that orders these values.   
    */
   int last_nonzero_column = 0;
   for (i = 0; i < 3; i++) {
      if ( !column_titles[i].empty() ) {
         last_nonzero_column = column_ids[i];
      }
   }
   int *ordered_list = new int[d_num_timers+1];
   buildOrderedList(timer_values,
                    last_nonzero_column, 
                    ordered_list,
                    d_num_timers);

   /*
    * Tack on TOTAL TIME to end of ordered list
    */
   ordered_list[d_num_timers] = d_num_timers;

   /*
    * Now output the rows of the table.
    */ 
   for (int j = 0; j < d_num_timers+1; j++) {
      n = ordered_list[j];

      /*
       * Check the print threshold to see if we should print this timer.
       */
      double frac = computePercentageDouble(
                       timer_values[n][last_nonzero_column],
		       timer_values[d_num_timers][last_nonzero_column] );

      if (frac > d_print_threshold) {

         os << setw(maxlen+3) << timer_names[n].c_str() << ' ';

         /*
          * Print columns.
          */ 
         for (i = 0; i < 4; i++) {

            /*
             * Print column values only title is non-null (i.e. not "")
             */ 
            if ( !column_titles[i].empty() ) {

               /*
                * Print percentages for columns 0-2
                */ 
	       if (i < 3) {
                  int k = column_ids[i];

                  if (d_print_percentage) {
#ifndef LACKS_SSTREAM
                     int perc = computePercentageInt(timer_values[n][k],
                                              timer_values[d_num_timers][k]);
                     ostringstream out;
                     out << timer_values[n][k] << " (" << perc << "%)";
                     os << setw(15) << out.str().c_str() << "  ";
#else
                     os << setw(15) << timer_values[n][k] << "  ";
#endif
                  } else {
                     os << setw(15) << timer_values[n][k] << "  ";
                  }

               } else {

                  /*
                   * Print column 3 - the processor holding processor ID 
                   * with max times (don't do for TOTAL TIME - this is 
                   * meaningless since all processors are synchronized 
                   * before and after this call).
                   */ 
	          if (n < d_num_timers) {
                     os << setw(15) << max_processor_id[n][max_array_id];
                  }

               } // column 3

            }  // if column title is non-null

         } // loop over columns

         os << endl;

      } //  matches d_print_threshold conditions

   } // loop over timers

   delete[] ordered_list;

   os << ascii_line  << endl;
   os.setf(ios::right);
}

void TimerManager::printOverhead(
                   const Array<string> timer_names,
		   const double timer_values[][18],
		   ostream& os)
{
   string ascii_line1 =  "++++++++++++++++++++++++++++++++++++++++"; 
   string ascii_line2 =  "++++++++++++++++++++++++++++++++++++++++\n"; 
   string ascii_line  = ascii_line1;
   ascii_line  += ascii_line2;

   /*
    * Left-justify all output in this method.
    */
   os.setf(ios::left);

   os << ascii_line 
      << "TIMER OVERHEAD STATISTICS \n"
      << ascii_line;
   
   int n,i;

   /*
    * Determine maximum name length for formatting
    */  
   int maxlen = 10;
   for (n = 0; n < d_num_timers; n++) { 
      i = int(timer_names[n].size());
      if (i > maxlen) maxlen = i;      
   }

   /*
    * Print table header
    */
   os << setw(maxlen+3) << "Timer Name"
      << setw(25) << "Number Accesses" << "  "
      << setw(25) << "Estimated Cost"
      << endl;

   /*
    * Compute totals: total number of REGISTERED accesses and total cost.  
    * Total cost includes inactive timer costs.
    */
   int total_inactive_accesses = 0;
   for (i = 0; i < d_num_inactive_timers; i++) {
      total_inactive_accesses += d_inactive_timers[i]->getNumberAccesses();
   }
   double est_cost = double(total_inactive_accesses) * TIMER_UNREG_ACCESS_TIME;
   double total_est_cost = est_cost;

   int total_accesses = 0;
   for (n = 0; n < d_num_timers; n++) {
      total_accesses += d_timers[n]->getNumberAccesses(); 
   }
   est_cost = double(total_accesses)*TIMER_REG_ACCESS_TIME;
  
   /*
    * If we are keeping exclusive or concurrent times, each access costs
    * roughly four times as much.  Make this correction here...
    */
   if (d_print_exclusive || d_print_concurrent) {
      est_cost *= 4.;
   }
   total_est_cost += est_cost; 

   /*
    * Output the rows of the table.  Start first with the inactive timers...
    */ 
   int num_accesses = total_inactive_accesses;
   est_cost = num_accesses * TIMER_UNREG_ACCESS_TIME;
   int perc = computePercentageInt(est_cost, total_est_cost);

   os << setw(maxlen+3) << "inactive timers"
      << setw(25) << num_accesses << "  ";
#ifndef LACKS_SSTREAM
   ostringstream out;
   out << est_cost << " (" << perc << "%)";
   os << setw(25) << out.str().c_str();
#else
   os << setw(25) << est_cost << " (" << perc << "%)";
#endif
   os << endl;

   /* 
    * Now print the rest of the timers.  While we are cycling through them, 
    * add up the total cost and print it at the end...
    */

   for (n = 0; n < d_num_timers; n++) {

      int num_accesses = d_timers[n]->getNumberAccesses();
      est_cost = num_accesses * TIMER_REG_ACCESS_TIME;

      /*
       * If we are keeping exclusive or concurrent times, each access costs
       * roughly four times as much.  Make this correction here...
       */
      if (d_print_exclusive || d_print_concurrent) {
	 est_cost *= 4.;
      }

      perc = computePercentageInt(est_cost, total_est_cost);

      os  << setw(maxlen+3) << timer_names[n].c_str()
          << setw(25) << num_accesses << "  ";
#ifndef LACKS_SSTREAM
      ostringstream out2;
      out2 << est_cost << " (" << perc << "%)";
      os   << setw(25) << out2.str().c_str();
#else
      os << setw(25) << est_cost << " (" << perc << "%)";
#endif
      os << endl;
   } 

   /* 
    * Output the totals.
    */
   os  << setw(maxlen+3) << "TOTAL:"
       << setw(25) << total_accesses << "  "
       << setw(25) << total_est_cost 
       << "\n" << endl; 

   /*
    * Compare the total estimated cost with overall program wallclock time.  
    * If it is a significant percentage (> 5%) print a warning
    */
   double perc_dbl = computePercentageDouble(total_est_cost,
                                             timer_values[d_num_timers][11]);

   os << "Estimated Timer Costs as a percentage of overall Wallclock Time: "
      << perc_dbl << "% \n";
   if (perc_dbl > 5.) {
     os << "WARNING:  TIMERS ARE USING A SIGNIFICANT FRACTION OF RUN TIME"
        << endl;
   }

   os << ascii_line  << endl;
   os.setf(ios::right);
}

void TimerManager::printConcurrent(ostream& os)
{
   string ascii_line1 =  "++++++++++++++++++++++++++++++++++++++++"; 
   string ascii_line2 =  "++++++++++++++++++++++++++++++++++++++++\n"; 
   string ascii_line  = ascii_line1;
   ascii_line  += ascii_line2;

   os << ascii_line 
      << "CONCURRENT TIMERS\n"
      << ascii_line;
   
   int n,i;

   /*
    * Determine maximum name length for formatting
    */  
   int maxlen = 10;
   for (n = 0; n < d_num_timers; n++) { 
      i = int((d_timers[n]->getName()).size());
      if (i > maxlen) maxlen = i;
      
   }

   /*
    * Print table header
    */
   os << setw(maxlen+3) << "Timer Name"
      << setw(25) << "Nested Timers" 
      << endl;

   /*
    * Output the rows of the table.
    */ 
   
   for (n = 0; n < d_num_timers; n++) {

      os << setw(maxlen+3) << d_timers[n]->getName().c_str();

      bool* bool_vec = d_timers[n]->getConcurrentTimerVector();
      int count = 0;
      for (i = 0; i < TIMER_MANAGER_MAX_TIMERS; i++) {
         if (bool_vec[i]) count++;
      }
      if (count == 0) {
	os << setw(25) << "none " << endl;
      } else {       
        /*
         * Format it like:    Timer Name      Concurrent Timer #1
         *                                    Concurrent Timer #2
         *                                    ... 
         * Use "count" variable defined above to identify the first
         * line or subsequent lines.
         */
        count = 0;
        for (int j = 0; j < TIMER_MANAGER_MAX_TIMERS; j++) {
           if (bool_vec[j]) {
              if (count == 0) {
                 os << setw(25) << d_timers[j]->getName().c_str() << endl;
              } else {
                 os << setw(maxlen+3) << " " 
                    << d_timers[j]->getName().c_str() << endl;
              }
              count++;
           }
        }
      }

   }
   os << ascii_line  << endl;
}


void TimerManager::checkConsistencyAcrossProcessors()
{
   /*
    * Due to the difficulty of comparing strings using MPI calls,
    * we do a rough consistency check of
    * 1. the number of timers and
    * 2. the length of each timer name.
    *
    * Steps:
    * 1. Do global reductions to get the max number of timers
    *    and the max lengths of each timer name.
    * 2. Issue a warning if the number of timers is inconsistent.
    *    This inconsistency would be found on all processes
    *    except those with the biggest number of timers.
    * 3. Issue a warning for each individual timer if
    *    its name length is less than the max length of
    *    all timers at the same index in the timer manager.
    *    Even if the number of timers are consistent, this
    *    would find wrong timer orderings or inconsistent
    *    timer names, unless the errors are for timer names
    *    with identical lengths.
    * 4. Go global reductions to get the number of inconsistencies
    *    of other processes.  Turn off printing of sum and max
    *    if any processes has inconsistencies.
    *
    * In the future, we may want to convert the strings into
    * their MD5 signatures and compare those as integers.
    */

   int max_num_timers = MPI::maxReduction( d_num_timers );

   Array<int> max_timer_lengths(max_num_timers);
   Array<int> rank_of_max(max_num_timers);

   for ( int i=0; i<max_num_timers; ++i ) {
      max_timer_lengths[i] =
         i < d_num_timers ? int(d_timers[i]->getName().size()) : 0;
   }

   MPI::maxReduction( max_timer_lengths.getPointer(),
                      max_num_timers,
                      rank_of_max.getPointer() );

   int inconsistency_count = 0;

   if ( max_num_timers > d_num_timers ) {
      TBOX_WARNING("Timer selections across processors were determined to be"
                   << "\ninconsistent.  This processor has only "
                   << d_num_timers << " while some has " << max_num_timers
                   << ".\nThe consistency check"
                   <<"\nwill continue for this process, but checking only\n"
                   << d_num_timers <<" timers."
                   << "\nIt is not possible to print global"
                   << "\nsummed or max timer information.");
      ++inconsistency_count;
   }

   for (int i = 0; i < d_num_timers; i++) {
      if ( max_timer_lengths[i] != int(d_timers[i]->getName().size()) ) {
         TBOX_WARNING("Timer[" << i << "]: " << d_timers[i]->getName() 
                      << "\nis not consistent across all processors."
                      << "\nOther timer[" << i << "] has up to "
                      << max_timer_lengths[i] << " characters in their names."
                      << "\nIt is not possible to print global"
                      << "\nsummed or max timer information."
                      );
         ++inconsistency_count;
      }
   }

   int max_inconsistency_count = MPI::maxReduction( inconsistency_count );
   if ( max_inconsistency_count > 0 ) {
      d_print_summed = false; 
      d_print_max = false;
      if ( inconsistency_count == 0 ) {
         TBOX_WARNING("Though this process found no timer inconsistencies,"
                      <<"\nother processes did.  It is not possible to print"
                      <<"\nglobal summed or max timer information.");
      }
   }

  /*
   * NOTE:  It might be nice to someday add the capability to remove the
   * inconsistent timers and print the max/summed values of the 
   * consistent ones.   Unfortunately, this is tough to implement.  If it
   * were just a matter of comparing timer names across processors it would be
   * easy. But with MPI, only ints and doubles can be exchanged across 
   * processors so it is difficult to make string comparisons.
   * It is possible to compare the MD5 sum of the strings,
   * but that may make SAMRAI dependent on the MD5 library.
   */
}



void TimerManager::buildTimerArrays(
                   double timer_values[][18], 
                   int max_processor_id[][2],
                   Array<string> timer_names)
{
  /*
   * timer_values - 2D array dimensioned [d_num_timers][18]
   *     For each timer, there are 18 potential values which may be of 
   *     interest.  This array stores them if they are requested.
   * max_processor_id - 2D array dimensioned [d_num_timers][2]
   *     Holds the value of the processor that used the maximum amount 
   *     of time.  [0] is for exclusive time, while [1] is for total time.
   */
   int i,j,k,n;

   /*
    * Initialize arrays
    */
   for (n = 0; n < d_num_timers+1; n++) {
      timer_names[n] = "";
      max_processor_id[n][0] = 0;
      max_processor_id[n][1] = 0;
      for (i = 0; i < 18; i++) {
         timer_values[n][i] = 0.;
      }
   }    

   /*
    * Build arrays.
    */
   for (n = 0; n < d_num_timers; n++) { 
      timer_names[n] = d_timers[n]->getName();
      
      /*
       *  Build timer_values[n][m] array:
       *    m = 0 :  processor exclusive user time
       *    m = 1 :  processor exclusive sys time
       *    m = 2 :  processor exclusive wall time
       *    m = 3 :  summed exclusive user time
       *    m = 4 :  summed exclusive sys time
       *    m = 5 :  summed exclusive wall time
       *    m = 6 :  max exclusive user time
       *    m = 7 :  max exclusive sys time
       *    m = 8 :  max exclusive wall time
       *    m = 9 :  processor total user time
       *    m = 10 :  processor total sys time
       *    m = 11 :  processor total wall time
       *    m = 12 :  summed total user time
       *    m = 13 :  summed total sys time
       *    m = 14 :  summed total wall time
       *    m = 15 :  max total user time
       *    m = 16 :  max total sys time
       *    m = 17 :  max total wall time
       */

      for (k = 0; k < 2; k++) {
	for (j = 0; j < 3; j++) {
             
            if ( (k == 0 && d_print_exclusive) ||
                 (k == 1 && d_print_total) ) {

	       if ( (j == 0 && d_print_processor) ||
                    (j == 1 && d_print_summed) ||
                    (j == 2 && d_print_max) ) {

		    if (k == 0 && j == 0) {
		      if (d_print_user) {
                         timer_values[n][0] =  
                                     d_timers[n]->getExclusiveUserTime(); 
                      }
		      if (d_print_sys) {
                         timer_values[n][1] =  
                                     d_timers[n]->getExclusiveSystemTime(); 
                      }
		      if (d_print_wall) {
                         timer_values[n][2] =  
                            d_timers[n]->getExclusiveWallclockTime(); 
                      }
                    } else if (k == 0 && j == 1) {
		      if (d_print_user) {
                         double user_time = 
                            d_timers[n]->getExclusiveUserTime();
                         timer_values[n][3] = 
                            MPI::sumReduction(user_time);
                      }
		      if (d_print_sys) {
                         double sys_time = 
                            d_timers[n]->getExclusiveSystemTime();
                         timer_values[n][4] =   
                            MPI::sumReduction(sys_time);
                      }
		      if (d_print_wall) {
			 double wall_time = 
			    d_timers[n]->getExclusiveWallclockTime();
                         timer_values[n][5] =   
			    MPI::sumReduction(wall_time);
                      }
                    } else if (k == 0 && j == 2) {
		      if (d_print_user) {
                         double user_time = 
                            d_timers[n]->getExclusiveUserTime();
                         timer_values[n][6] =   
			    MPI::maxReduction(user_time);
                      }
		      if (d_print_sys) {
                         double sys_time = 
                            d_timers[n]->getExclusiveSystemTime();
                         timer_values[n][7] =   
			    MPI::maxReduction(sys_time);
                      }
		      if (d_print_wall) {
                         int max_id;
			 double wall_time = 
			    d_timers[n]->getExclusiveWallclockTime();
                         timer_values[n][8] =
			    MPI::maxReduction(wall_time,&max_id);
                         max_processor_id[n][0] = max_id;
                      }

		   } else if (k == 1 && j == 0) {
		      if (d_print_user) {
                         timer_values[n][9] =  
                                     d_timers[n]->getTotalUserTime(); 
                      }
		      if (d_print_sys) {
                         timer_values[n][10] =  
                                     d_timers[n]->getTotalSystemTime(); 
                      }
		      if (d_print_wall) {
                         timer_values[n][11] =  
                                     d_timers[n]->getTotalWallclockTime(); 
                      }
                    } else if (k == 1 && j == 1) {
		      if (d_print_user) {
                         double user_time = 
                            d_timers[n]->getTotalUserTime();
                         timer_values[n][12] =   
			    MPI::sumReduction(user_time);
                      }
		      if (d_print_sys) {
                         double sys_time = 
                            d_timers[n]->getTotalSystemTime();
                         timer_values[n][13] =   
			    MPI::sumReduction(sys_time);
                      }
		      if (d_print_wall) {
			 double wall_time = 
			    d_timers[n]->getTotalWallclockTime();
                         timer_values[n][14] =   
			    MPI::sumReduction(wall_time);
                      }
                    } else if (k == 1 && j == 2) {
		      if (d_print_user) {
                         double user_time = 
                            d_timers[n]->getTotalUserTime();
                         timer_values[n][15] =   
			    MPI::maxReduction(user_time);
                      }
		      if (d_print_sys) {
                         double sys_time = 
                            d_timers[n]->getTotalSystemTime();
                         timer_values[n][16] =   
			    MPI::maxReduction(sys_time);
                      }
		      if (d_print_wall) {
                         int max_id;
			 double wall_time = 
			    d_timers[n]->getTotalWallclockTime();
                         timer_values[n][17] =
			    MPI::maxReduction(wall_time,&max_id);
                         max_processor_id[n][1] = max_id;
                      }

                    }

 
	         } // if j
	      } // if k

	} // loop over j
      } // loop over k 

   } // loop over n

   /*
    * Store main_timer data in timer_values[d_num_timers][] location.  Max
    * time and exclusive time are not determined since these don't really 
    * mean anything for an overall measurement of run time.
    */
    timer_names[d_num_timers] = "TOTAL RUN TIME:";
    if (d_print_user) {
       double main_time = d_main_timer->getTotalUserTime();
       timer_values[d_num_timers][0] = main_time;
       timer_values[d_num_timers][3] = MPI::sumReduction(main_time);
       timer_values[d_num_timers][6] = main_time;
       timer_values[d_num_timers][9] = main_time;
       timer_values[d_num_timers][12] = MPI::sumReduction(main_time);
       timer_values[d_num_timers][15] = main_time;
    }
    if (d_print_sys) {
       double main_time = d_main_timer->getTotalSystemTime();
       timer_values[d_num_timers][1] = main_time;
       timer_values[d_num_timers][4] = MPI::sumReduction(main_time);
       timer_values[d_num_timers][7] = main_time;
       timer_values[d_num_timers][10] = main_time;
       timer_values[d_num_timers][13] = MPI::sumReduction(main_time);
       timer_values[d_num_timers][16] = main_time;
    }
    if (d_print_wall) {
       double main_time = d_main_timer->getTotalWallclockTime();
       timer_values[d_num_timers][2] = main_time;
       timer_values[d_num_timers][5] = MPI::sumReduction(main_time);
       timer_values[d_num_timers][8] = main_time;
       timer_values[d_num_timers][11] = main_time;
       timer_values[d_num_timers][14] = MPI::sumReduction(main_time);
       timer_values[d_num_timers][17] = main_time;
    }

}



/*
*************************************************************************
*                                                                       *
* Build ordered_list which specifies order of timers - max to min.      *
*                                                                       *
*************************************************************************
*/
void TimerManager::buildOrderedList(
      const double timer_values[][18], 
      const int column,
      int index[],
      const int array_size)
{
   /*
    * initialize the arrays
    */
   Array<double> timer_vals;
   timer_vals.resizeArray(array_size);
   for (int i = 0; i < array_size; i++) {
      index[i] = i;
      timer_vals[i] = timer_values[i][column];
   }

   /*
    * Do a quicksort on timer_values array to build index array 
    * ordered_list.
    */
   quicksort(timer_vals, index, 0, array_size-1);
  
}

/*
*************************************************************************
*                                                                       *
* Sort array a largest to smallest using quicksort algorithm.           *
*                                                                       *
*************************************************************************
*/
void TimerManager::quicksort(
      const Array<double>& a, int index[], int lo, int hi) 
{
   if (hi <= lo) return;

   /*
    * Put a[i] into position for i between lo and hi
    * (i.e. pivot point)
    */
   int i = lo - 1;
   int j = hi;
   double v = a[index[hi]];
   for (;;) {
      while (a[index[++i]] > v) 
	 NULL_STATEMENT;
      while (v > a[index[--j]]) {
         if (j == lo) break;
      }
      if (i >= j) break;

      // exchange i, j indices
      int temp = index[i];
      index[i] = index[j];
      index[j] = temp;
   }
   // exchange i, hi indices
   int temp = index[i];
   index[i] = index[hi];
   index[hi] = temp;

   quicksort( a, index, lo, i-1 );
   quicksort( a, index, i+1, hi );
}

/*
*************************************************************************
*                                                                       *
* Operation performed many times throughout print routines. We have     *
* to have some safety checks to avoid divide-by-zero errors in some     *
* cases.  Thus, I just made it a function.                              *
*                                                                       *
*************************************************************************
*/
int TimerManager::computePercentageInt( const double frac, 
					     const double tot) 
{
   /*
    *  Put a cap on the percentage at 1000.  If tot = 0, this if
    *  test should catch it.
    */
   int perc = 0;
   if (tot > 0.1*frac) {
      perc = int(frac/tot*100.);
   } else {
      perc = 1000; 
   }
   return (perc);
}

double TimerManager::computePercentageDouble( const double frac, 
						   const double tot) 
{
   /*
    *  Put a cap on the percentage at 1000.  If tot = 0, this if
    *  test should catch it.
    */
   double perc = 0;
   if (tot > 0.1*frac) {
      perc = frac/tot*100.;
   } else {
      perc = 1000; 
   }
   return (perc);
}


/*
*************************************************************************
*                                                                       *
* Private member function for processing input data.                    *
*                                                                       *
*************************************************************************
*/

void TimerManager::getFromInput(
   Pointer<Database> input_db)
{

   if (!input_db.isNull()) {

      if (input_db->keyExists("print_exclusive")) {
         d_print_exclusive = input_db->getBool("print_exclusive");
      } else {
         d_print_exclusive = input_db->getBoolWithDefault("print_exclusive",
                                                          d_print_exclusive);
      }

      if (input_db->keyExists("print_total")) {
         d_print_total = input_db->getBool("print_total");
      } else {
	d_print_total = input_db->getBoolWithDefault("print_total",
                                              d_print_total);
      }

      if (input_db->keyExists("print_processor")) {
	d_print_processor = input_db->getBool("print_processor");
      } else {
	d_print_processor = input_db->getBoolWithDefault("print_processor",
                                                  d_print_processor);
      }

      if (input_db->keyExists("print_max")) {
	d_print_max = input_db->getBool("print_max");
      } else {
	d_print_max = input_db->getBoolWithDefault("print_max",
                                            d_print_max);
      }

      if (input_db->keyExists("print_summed")) {
	d_print_summed = input_db->getBool("print_summed");
      } else {
	d_print_summed = input_db->getBoolWithDefault("print_summed",
                                               d_print_summed);
      }

      if (input_db->keyExists("print_user")) {
	d_print_user = input_db->getBool("print_user");
      } else {
	d_print_user = input_db->getBoolWithDefault("print_user",
                                             d_print_user);
      }

      if (input_db->keyExists("print_sys")) {
	d_print_sys = input_db->getBool("print_sys");
      } else {
	d_print_sys = input_db->getBoolWithDefault("print_sys",
                                            d_print_sys);
      }

      if (input_db->keyExists("print_wall")) {
	d_print_wall = input_db->getBool("print_wall");
      } else {
	d_print_wall = input_db->getBoolWithDefault("print_wall",
                                             d_print_wall);
      }

      if (input_db->keyExists("print_percentage")) {
	d_print_percentage = input_db->getBool("print_percentage");
      } else {
	d_print_percentage = input_db->getBoolWithDefault("print_percentage",
                                             d_print_percentage);
      }

      if (input_db->keyExists("print_concurrent")) {
	d_print_concurrent = input_db->getBool("print_concurrent");
      } else {
	d_print_concurrent = input_db->getBoolWithDefault("print_concurrent",
                                                  d_print_concurrent);
      }

      if (input_db->keyExists("print_timer_overhead")) {
	d_print_timer_overhead = input_db->getBool("print_timer_overhead");
      } else {
	d_print_timer_overhead = input_db->
                               getBoolWithDefault("print_timer_overhead",
                                                  d_print_timer_overhead);
      }

      if (input_db->keyExists("print_threshold")) {
	d_print_threshold =
         input_db->getDouble("print_threshold");
      } else {
	d_print_threshold =
	  input_db->getDoubleWithDefault("print_threshold",
                                         d_print_threshold);
      }

      Array<string> timer_list;
      if (input_db->keyExists("timer_list")) {
	timer_list = input_db->getStringArray("timer_list");
      }

      /*
       *  Step thru the input list and call addTimerToNameLists to add
       *  the input file entry to the d_package_names,
       *  d_class_names, and d_class_method_names lists.
       */
      for (int i = 0; i < timer_list.getSize(); i++) {
	string entry = timer_list[i];
	addTimerToNameLists(entry);
      }
      d_length_package_names      = d_package_names.getNumberItems();
      d_length_class_names        = d_class_names.getNumberItems();
      d_length_class_method_names = d_class_method_names.getNumberItems();

   }

}

void TimerManager::addTimerToNameLists(
   const string& name)
{
   /*
    * Evaluate whether the name is a package, class, or class::method 
    * combination.  This parser supports inputs of the form:
    *
    *    *::*::*         - ALL timers added
    *    Package::*::*   - "Package" added to package list.
    *    Class           - "Class" added to class list.  
    *    *::Class        - "Class" added to class list.  
    *    Class::*        - "Class" added to class list. 
    *    *::Class::*     - "Class" added to class list. 
    *    Package::Class::method  - "Class::method" put to class_method list
    *    Class::method   - "Class::method" put to class_method list                                         
    */

   string::size_type position, string_length;

   /*
    *  Step thru the input list and form the d_package_names,
    *  d_class_names, and d_class_method_names lists.
    */

   if (!name.empty()) {    // Nested if #1

      string entry = name;

      /*
       *  Once we have determined whether the entry is a package,
       *  class, or class::method, use this bool to jump to the next 
       *  loop entry.
       */
      bool determined_entry = false;

      /*
       * Check if its a wildcard entry - "*::*::*".  If so, add all package
       * names to the package name list.
       */
      position = entry.find("*::*::*");  // if not found, "position" runs off
                                         // end of entry so pos > entry.size()
      if (position < entry.size()) {
         d_package_names.addItem("algs");
         d_package_names.addItem("apps");
         d_package_names.addItem("appu");
         d_package_names.addItem("geom");
         d_package_names.addItem("hier");
         d_package_names.addItem("math");
         d_package_names.addItem("mesh");
         d_package_names.addItem("pdat");
         d_package_names.addItem("solv");
         d_package_names.addItem("tbox");
         d_package_names.addItem("xfer");
         determined_entry = true;
      }

      /*
       * Is it a package?  Look for "::*::*" string.  If its there, 
       * parse it off and add the package to the package list. 
       */
      if (!determined_entry) {    
         position = entry.find("::*::*");
         if (position < entry.size()) {
            entry = entry.substr(0,position);
            d_package_names.addItem(entry);
            determined_entry = true;
         }
      }
      

      if (!determined_entry) {    // Nested if #2

         /*
          * Is it a class?  If it doesn't have any "::", it must be a class.
          */
         position = entry.find("::");
         if (position > entry.size()) {
            d_class_names.addItem(entry);
            determined_entry = true;
         }
         if (!determined_entry) {    // Nested if #3

            /*
             * At this point, we know the entry has a "::" but wasn't 
             * identified as a package.  There are several options that 
	     * might make Foo a class entry:
             *  1) Foo::*
             *  2) *::Foo::*
             *  3) Package::Foo::*
             *  4) *::Foo
             * Parse these as follows:  First, look for existence of "::*" 
             * at the end of the entry.  This will identify the first 3 
             * options.  Next look for existence of "*::" at front of the 
             * string.  This will identify the fourth choice.
             *
             * Check for first three options... 
             */
            string_length = entry.size();
            string substring = entry.substr(string_length-3,string_length);
            if (substring == "::*") {
               entry = entry.substr(0,string_length-3);

               /*
                * If a preceeding "::" exists at the front of the entry 
                * (i.e. option 2 and 3), parse off anything before it.
                */
               position = entry.find("::");
               if (position < entry.size()) {
                  entry = entry.substr(position+2);
               }
               d_class_names.addItem(entry);
               determined_entry = true;
            }

            if (!determined_entry) {    // Nested if #4

            /*
             * Check for option 4.  The entry has a preceeding *::. Do not
             * accept case where there is a second "::" followed by anything
             * but "*", since this is a class::method combination.
             *
             */
               substring = entry.substr(0,3);
               if (substring == "*::") {
                  entry = entry.substr(3);
                  position = entry.find("::");

                  /*
                   * There is no second "::".  Accept the entry as a class.
                   */
                  if (position > entry.size()) {
                     d_class_names.addItem(entry);
                     determined_entry = true;
                  } else {

                  /*
                   * There is a second "::".  See if it is followed by a 
		   * "*".  If so, parse off the "::*" and accept entry as 
                   * a class.  If not, let it be determined below to be a 
                   * class::method entry.
                   */
                     string_length = entry.size();
                     substring = entry.substr(string_length-1,string_length);
                     if (substring == "*") {
                        entry = entry.substr(0,string_length-3);
                        d_class_names.addItem(entry);
                        determined_entry = true;
                     }
                  }  
               }

               if (!determined_entry) {    // Nested if #5

               /*
                * The entry has not been identified as either a package or 
                * a class.  It must be a class::method combination. There 
                * are three options for entering class::method combinations:
                *  1) Package::Foo::method
                *  2) *::Foo::method
                *  3) Foo::method
                * We only want to maintain "Foo::method" in the package 
                * list.  Check first if there are two "::" in the entry.  
                * If there are, parse of whatever is in front of the 
                * first "::".  If not, just use the entry as is.
                */
                  position = entry.find("::");
                  if (position < entry.size()) {
                     substring = entry.substr(position+2);
                     position = substring.find("::");
                     if (position < substring.size()) {

                        /*
                         * There *are* two "::" so entry must contain a 
                         * package.  Parse it off. 
                         */ 
                        position = entry.find("::");
                        entry = entry.substr(position+2);
                     }
                  }
                  d_class_method_names.addItem(entry);

               } // Nested if #5    

            } // Nested if #4   

         } // Nested if #3

      } // Nested if #2

   } // Nested if #1
}


}
}
