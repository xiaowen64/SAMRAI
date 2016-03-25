//
// File:        main.C
// Description: Test program to demonstrate/test timers.
//

#include "SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h"
#include "Foo.h"
#include "tbox/InputManager.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include <string>
using namespace std;
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#ifndef NULL
#define NULL (0)
#endif

using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   tbox::PIO::logAllNodes("Timer.log");      

   string input_filename;
   string restart_dirname;
   int    restore_num = 0;

   bool is_from_restart = false;
 
   if ( (argc != 2) && (argc != 4) ) {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
           << "<restart dir> <restore number> [options]\n"
           << "  options:\n"
           << "  none at this time"
           << endl;
        tbox::MPI::abort();
        return (-1);
   } else {
      input_filename = argv[1];
      if (argc == 4) {
         restart_dirname = argv[2];
         restore_num = atoi(argv[3]);
 
         is_from_restart = true;
      }
   }

#ifndef HAVE_HDF5
   is_from_restart = false;
#endif

   int i;

   /*
    * Create an input database "input_db" and parse input file (specified
    * on the command line.
    */
   tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
   tbox::InputManager::getManager()->parseInputFile(input_filename,input_db);

   /*
    * Retrieve "Main" section of the input database.  Read ntimes, 
    * which is the number of times the functions are called, and 
    * depth, which is the depth of the exclusive timer tree.
    */
   tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

   int ntimes = 1;
   if (main_db->keyExists("ntimes")) {
      ntimes = main_db->getInteger("ntimes");
   }

   int exclusive_tree_depth = 1;
   if (main_db->keyExists("exclusive_tree_depth")) {
      exclusive_tree_depth = main_db->getInteger("exclusive_tree_depth");
   }

   /*
    * Open the restart file and read information from file into the
    * restart database.
    */
   if (is_from_restart) {
      tbox::RestartManager::getManager()->
         openRestartFile(restart_dirname, restore_num, tbox::MPI::getNodes());
   }

   tbox::Pointer<tbox::Database> restart_db =
        tbox::RestartManager::getManager()->getRootDatabase();

   /*
    * Create timer manager, reading timer list from input.
    */
   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   /*
    * Add a timer "manually" (that is, not thru the input file).
    */
   static tbox::Pointer<tbox::Timer> timer = tbox::TimerManager::getManager()->
      getTimer("apps::main::main");
   timer->start();

   /*
    * We no longer need the restart file so lets close it. 
    */
   tbox::RestartManager::getManager()->closeRestartFile();

   /*
    * Class Foo contains the functions we want to call.
    */
   Foo* foo = new Foo();

   /*
    * Check time to call function with timer name that is NOT
    * registered.  That is, time a NULL timer call.
    */
   static tbox::Pointer<tbox::Timer> timer_off = tbox::TimerManager::getManager()->
      getTimer("apps::main::timer_off");
   timer_off->start();
   for (i=0; i < ntimes; i++) {
      foo->timerOff();
   }
   timer_off->stop();

   /*
    * Check time to call function with timer name that IS
    * registered.  
    */
   static tbox::Pointer<tbox::Timer> timer_on = tbox::TimerManager::getManager()->
      getTimer("apps::main::timer_on");
   tbox::Pointer<tbox::Timer> dummy_timer = tbox::TimerManager::getManager()->
      getTimer("apps::Foo::timerOn()");
   timer_on->start();
   for (i=0; i < ntimes; i++) {
      foo->timerOn();
   }
   timer_on->stop();

   /*
    * Time to call tree-based set of exclusive timers. i.e.
    * Foo->zero() calls Foo->one(), which calls Foo->two(), ... 
    * and so forth until we reach specified "exclusive_tree_depth.
    */  
   static tbox::Pointer<tbox::Timer> timer_excl = tbox::TimerManager::getManager()->
      getTimer("apps::main::exclusive_timer");
   timer_excl->start();
   for (i=0; i < ntimes; i++) {
      foo->zero(exclusive_tree_depth);
   }
   timer_excl->stop();
   timer->stop();

   double eff = timer->computeLoadBalanceEfficiency();
   tbox::pout << "Load Balance eff: " << eff << "%" << endl;


//   tbox::TimerManager::getManager()->print(tbox::plog);

   tbox::TimerManager::getManager()->resetAllTimers();

   timer->start();

   /*
    * Check time to call function with timer name that is NOT
    * registered.  That is, time a NULL timer call.
    */
   timer_off->start();
   for (i=0; i < ntimes; i++) {
      foo->timerOff();
   }
   timer_off->stop();


   /*
    * Check time to call function with timer name that IS
    * registered.  
    */
   timer_on->start();
   for (i=0; i < ntimes; i++) {
      foo->timerOn();
   }
   timer_on->stop();

   /*
    * Time to call tree-based set of exclusive timers. i.e.
    * Foo->zero() calls Foo->one(), which calls Foo->two(), ... 
    * and so forth until we reach specified "exclusive_tree_depth.
    */  
   timer_excl->start();
   for (i=0; i < ntimes; i++) {
      foo->zero(exclusive_tree_depth);
   }
   timer_excl->stop();
   timer->stop();

   tbox::TimerManager::getManager()->print(tbox::plog);


   /*
    * We're done.  Write the restart file.  
    */
   string restart_write_dirname = "restart";
#ifdef HAVE_HDF5
   int timestep = 0;
   tbox::RestartManager::getManager()->writeRestartFile(
               restart_write_dirname,
               timestep);
#endif

   delete foo;

   tbox::pout << "\nPASSED:  timertest" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();
   return(0);
}


