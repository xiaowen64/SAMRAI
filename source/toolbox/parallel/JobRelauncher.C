/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/parallel/JobRelauncher.C $
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Job relauncher handling multiple relaunchable jobs.
 */

#include "tbox/JobRelauncher.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
   namespace tbox {


Pointer<Timer> JobRelauncher::t_compute;
Pointer<Timer> JobRelauncher::t_commwait;
Pointer<Timer> JobRelauncher::t_commwait_any;
Pointer<Timer> JobRelauncher::t_commwait_some;
Pointer<Timer> JobRelauncher::t_commwait_sync;


/*
*******************************************************************
*******************************************************************
*/
JobRelauncher::JobRelauncher()
{

   if ( t_compute.isNull() ) {
      TBOX_ASSERT( t_commwait_any.isNull() );
      TBOX_ASSERT( t_commwait_some.isNull() );
      TBOX_ASSERT( t_commwait_sync.isNull() );
      t_compute = TimerManager::getManager()->
         getTimer("tbox::JobRelauncher::compute");
      t_commwait = TimerManager::getManager()->
         getTimer("tbox::JobRelauncher::commwait");
      t_commwait_any = TimerManager::getManager()->
         getTimer("tbox::JobRelauncher::commwait_any");
      t_commwait_some = TimerManager::getManager()->
         getTimer("tbox::JobRelauncher::commwait_some");
      t_commwait_sync = TimerManager::getManager()->
         getTimer("tbox::JobRelauncher::commwait_sync");
      ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }

   return;
}


/*
*******************************************************************
*                                                                 *
* This method delegates running the algorithm to any of the       *
* methods named runAlgorithm_...(), which try to execute          *
* the BR algorithm with a specific implementation.                *
* The rest of this comment block refers to those methods.         *
* See also comments for the specific version runAlgorithm_...().  *
*                                                                 *
*******************************************************************
*/
void JobRelauncher::runAlgorithm( RelaunchableJob *initial_job )
{

   switch (d_algo_advance_mode) {
   case ROUND_ROBIN:
      runAlgorithm_round_robin( initial_job );
      break;
   case ADVANCE_ANY:
      runAlgorithm_advance_any( initial_job );
      break;
   case ADVANCE_SOME:
      runAlgorithm_advance_some( initial_job );
      break;
   case SYNCHRONOUS:
      runAlgorithm_synchronous( initial_job );
      break;
   }

   return;
}


/*
*********************************************************************
*                                                                   *
* This method sets up the relaunch queue with the initial job.      *
* It goes back and forth between jobs waiting                       *
* for communication and jobs waiting for non-communication          *
* until there is no more jobs.  Communication waits are             *
* performed using d_comm_stage.advanceAny().                          *
*                                                                   *
*********************************************************************
*/
void JobRelauncher::runAlgorithm_advance_any( RelaunchableJob *initial_job )
{
   RelaunchableJob* job_for_relaunch;
   d_relaunch_queue.appendItem(initial_job);

   do {

      t_compute->start();
      /*
       * Relaunch any execution node that needs it.
       *
       * Nodes that are waiting for its own communication to
       * complete are NOT in the relaunch queue because these
       * nodes are continued by the d_comm_stage.
       */
      while ( ! d_relaunch_queue.isEmpty() ) {
         job_for_relaunch = d_relaunch_queue.getFirstItem();
         d_relaunch_queue.removeFirstItem();
         job_for_relaunch->continueJob();
      }
      t_compute->stop();

      /*
       * There is no more node needing relaunching.
       * All current nodes are running but paused to wait
       * for communication.  In this next step, use
       * d_comm_stage to find the first running node that
       * has completed communication and continue it.
       * The continuing node may put more nodes on the
       * d_relaunch_queue, so go back and check the d_relaunch_queue
       * afterwards.
       */
      t_commwait->start();
      t_commwait_any->start();
      job_for_relaunch = d_comm_stage.advanceAny();
      t_commwait_any->stop();
      t_commwait->stop();
      if ( job_for_relaunch != NULL ) {
         t_compute->start();
         job_for_relaunch->continueJob();
         t_compute->stop();
      }

   } while ( ! d_relaunch_queue.isEmpty() || job_for_relaunch != NULL );

   return;
}


/*
***********************************************************************
*                                                                     *
* This method sets up the relaunch queue with the initial job.        *
* It goes back and forth between jobs waiting                         *
* for communication and jobs waiting for non-communication            *
* until there is no more jobs.  Communication waits are               *
* performed using d_comm_stage.advanceSome().                           *
*                                                                     *
***********************************************************************
*/
void JobRelauncher::runAlgorithm_advance_some( RelaunchableJob *initial_job )
{
   int n_comm_group_completed = 0;
   d_relaunch_queue.appendItem(initial_job);

   do {

      t_compute->start();
      while ( ! d_relaunch_queue.isEmpty() ) {
         RelaunchableJob *job_for_relaunch =
            d_relaunch_queue.getFirstItem();
         d_relaunch_queue.removeFirstItem();
         job_for_relaunch->continueJob();
      }
      t_compute->stop();

      Array<RelaunchableJob*> job_for_relaunch;
      t_commwait->start();
      t_commwait_some->start();
      n_comm_group_completed = d_comm_stage.advanceSome( job_for_relaunch );
      t_commwait_some->stop();
      t_commwait->stop();

      t_compute->start();
      for ( int n=0; n<n_comm_group_completed; ++n ) {
         if ( job_for_relaunch[n] != NULL ) {
            job_for_relaunch[n]->continueJob();
         }
      }
      t_compute->stop();

   } while ( ! d_relaunch_queue.isEmpty() || n_comm_group_completed > 0 );

   return;
}


/*
***********************************************************************
*                                                                     *
* This method sets up the relaunch queue with the initial job.        *
* Then it cycles through the relaunch queue until the queue           *
* is empty.  Jobs that did not finish are put back on the             *
* relaunch queue to be relaunched later.                              *
*                                                                     *
***********************************************************************
*/
void JobRelauncher::runAlgorithm_round_robin( RelaunchableJob *initial_job )
{
   d_relaunch_queue.appendItem(initial_job);

   do {

      RelaunchableJob *job_for_relaunch = d_relaunch_queue.getFirstItem();
      d_relaunch_queue.removeFirstItem();
      job_for_relaunch->continueJob();

      if ( job_for_relaunch->getJobState() ==
           RelaunchableJob::COMMUNICATION_WAIT ) {
         d_relaunch_queue.appendItem(job_for_relaunch);
      }

   } while ( ! d_relaunch_queue.isEmpty() );

   return;
}


/*
***********************************************************************
*                                                                     *
* This method sets up the relaunch queue with the initial job.        *
* It goes back and forth between jobs waiting                         *
* for communication and jobs waiting for non-communication            *
* until there is no more jobs.  Jobs waiting for communication        *
* is advanced past the communication, leading to a synchronous        *
* algorithm.                                                          *
*                                                                     *
***********************************************************************
*/
void JobRelauncher::runAlgorithm_synchronous( RelaunchableJob *initial_job )
{
   d_relaunch_queue.appendItem(initial_job);

   do {

      RelaunchableJob *job_for_relaunch = d_relaunch_queue.getFirstItem();

      d_relaunch_queue.removeFirstItem();

      if ( job_for_relaunch->getJobState() ==
           RelaunchableJob::NONCOMMUNICATION_WAIT ) {
         t_compute->start();
         job_for_relaunch->continueJob();
         t_compute->stop();
      }

      while ( job_for_relaunch->getJobState() ==
              RelaunchableJob::COMMUNICATION_WAIT ) {
         t_commwait->start();
         t_commwait_sync->start();
         job_for_relaunch->getCommunicationGroup()->waitOperation();
         t_commwait_sync->stop();
         t_commwait->stop();
         t_compute->start();
         job_for_relaunch->continueJob();
         t_compute->stop();
      }

   } while ( ! d_relaunch_queue.isEmpty() );

   return;
}






void JobRelauncher::setAlgorithmAdvanceMode(
   const std::string &algo_advance_mode )
{
   if ( algo_advance_mode == "ROUND_ROBIN" ) {
      d_algo_advance_mode = ROUND_ROBIN;
   }
   else if ( algo_advance_mode == "ADVANCE_ANY" ) {
      d_algo_advance_mode = ADVANCE_ANY;
   }
   else if ( algo_advance_mode == "ADVANCE_SOME" ) {
      d_algo_advance_mode = ADVANCE_SOME;
   }
   else if ( algo_advance_mode == "SYNCHRONOUS" ) {
      d_algo_advance_mode = SYNCHRONOUS;
   }
   else {
      TBOX_ERROR("No such algorithm choice: " << algo_advance_mode
                 <<"\n");
   }
   return;
}





AsyncCommStage &JobRelauncher::getCommStage()
{
   return d_comm_stage;
}

List<RelaunchableJob*> &JobRelauncher::getRelaunchQueue()
{
   return d_relaunch_queue;
}


/*
***************************************************************************
*                                                                         *
* Release static timers.  To be called by shutdown registry to make sure  *
* memory for timers does not leak.                                        *
*                                                                         *
***************************************************************************
*/
void JobRelauncher::freeTimers()
{
   t_compute.setNull();
   t_commwait.setNull();
   t_commwait_any.setNull();
   t_commwait_some.setNull();
   t_commwait_sync.setNull();
   return;
}


}
}
