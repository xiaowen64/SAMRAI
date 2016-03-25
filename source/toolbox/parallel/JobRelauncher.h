/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/parallel/JobRelauncher.h $
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Job relauncher handling multiple relaunchable jobs.
 */

#ifndef included_tbox_JobRelauncher
#define included_tbox_JobRelauncher

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_AsyncCommStage
#include "tbox/AsyncCommStage.h"
#endif

#ifndef included_tbox_RelaunchableJob
#include "tbox/RelaunchableJob.h"
#endif

#ifndef included_tbox_List
#include "tbox/List.h"
#endif

#ifndef included_String
#include <string>
#define included_String
#endif

#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
   namespace tbox {


/*!
 * @brief Manages an algorithm consisting of multiple relaunchable jobs.
 *
 * The jobs, defined by the RelaunchableJob base class,
 * may be paused to wait for communication or non-communication
 * to finish and then be restarted.
 */
class JobRelauncher
{

public:

   JobRelauncher();


   /*!
    * @brief Begin the algorithm.
    *
    * The algorithm begins with the job given.  It may
    * add more jobs to the relaunch queue as it progresses.
    * See getRelaunchQueue().
    *
    * @param initial_job The initial job in the algorithm.
    */
   void runAlgorithm( RelaunchableJob *initial_job );


   /*!
    * @brief Set the mode for advancing the asynchronous implementation.
    *
    * Choices are:
    * - "SYNCHRONOUS" -> wait for each communication stage to complete
    *   before moving on, thus resulting in synchronous execution.
    * - "ROUND_ROBIN" -> check for completed communication stages in
    *   round-robin fashion instead of waiting for a specific one.
    * - "ADVANCE_ANY" -> advance an execution node through its
    *   communication stage by using AsyncCommStage::advanceAny().
    * - "ADVANCE_SOME" -> advance an execution node through its
    *   communication stage by using AsyncCommStage::advanceSome().
    *
    * The default is "ADVANCE_SOME".  This generally is fastest,
    * but the other modes may be better for debugging.
    *
    * Asynchronous modes are NOT guaranteed to compute the output
    * graph nodes in any particular order.  The order depends on
    * the ordering of message completion, which is not deterministic.
    * If you require consistent outputs, we suggest you have a scheme
    * for reordering your output.
    */
   void setAlgorithmAdvanceMode( const std::string &algo_advance_mode );


   /*!
    * @brief Return the communication stage managing jobs that
    * are waiting for communication to finish.
    */
   AsyncCommStage &getCommStage();

   /*!
    * @brief Return the queue of jobs waiting for non-communications
    * work.
    */
   List<RelaunchableJob*> &getRelaunchQueue();



private:


   /*!
    * @brief Method for advancing the algorithm.
    *
    * Each corresponds to a choice permitted by setAlgorithmAdvanceMode().
    */
   enum AlgoAdvanceMode { ROUND_ROBIN,
                          ADVANCE_ANY,
                          ADVANCE_SOME,
                          SYNCHRONOUS };


   /*!
    * @brief Version of runAlgorithm() that waits for each
    * communication stage to complete before moving on,
    * thus resulting in synchronous execution.
    */
   void runAlgorithm_synchronous( RelaunchableJob *initial_job );

   /*!
    * @brief Version of runAlgorithm() checking communication
    * stages in round-robin fashion instead of using
    * AsyncCommStage.
    */
   void runAlgorithm_round_robin( RelaunchableJob *initial_job );

   /*!
    * @brief Version of runAlgorithm() that advances an
    * execution node through its communication stage by
    * using AsyncCommStage::advanceAny().
    */
   void runAlgorithm_advance_any( RelaunchableJob *initial_job );

   /*!
    * @brief Version of runAlgorithm() that advances an
    * execution node through its communication stage by
    * using AsyncCommStage::advanceSome().
    */
   void runAlgorithm_advance_some( RelaunchableJob *initial_job );


   /*!
    * Free static timers.
    *
    * To be called by shutdown registry to make sure
    * memory for timers does not leak.
    */
   static void freeTimers();


   /*!
    * @brief Queue on which to append jobs to be
    * launched or relaunched.
    *
    * What is placed on the queue depends on the version of
    * runAlgorithm_...() used.
    */
   List<RelaunchableJob*> d_relaunch_queue;

   /*!
    * @brief Stage handling multiple asynchronous communication groups.
    */
   AsyncCommStage d_comm_stage;

   AlgoAdvanceMode d_algo_advance_mode;

   /*!
    * @brief Compute timer.
    */
   static Pointer<Timer> t_compute;
   /*!
    * @brief Communication-wait timer for all modes.
    */
   static Pointer<Timer> t_commwait;
   /*!
    * @brief Communication-wait timer for ADVANCE_ANY mode.
    */
   static Pointer<Timer> t_commwait_any;
   /*!
    * @brief Communication-wait timer for ADVANCE_SOME mode.
    */
   static Pointer<Timer> t_commwait_some;
   /*!
    * @brief Communication-wait timer for SYNCHRONOUS mode.
    */
   static Pointer<Timer> t_commwait_sync;
};

}
}

#endif  // included_tbox_JobRelauncher
