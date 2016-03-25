/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/parallel/AsyncCommStage.C $
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Support for coordinating multiple asynchronous communications
 */

#include "tbox/AsyncCommStage.h"


#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include STL_SSTREAM_HEADER_FILE

using namespace std;
namespace SAMRAI {
   namespace tbox {



#ifdef HAVE_MPI

Pointer<Timer> AsyncCommStage::t_wait_any;
Pointer<Timer> AsyncCommStage::t_wait_some;


AsyncCommStage::AsyncCommStage()
   :
   d_n_group(0),
   d_n_req(0),
   d_group(0, true),
   d_req(0, true),
   d_req_to_group(0, true),
   d_group_to_req(0, true)
{
   if ( t_wait_any.isNull() ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( t_wait_some.isNull() );
#endif
      t_wait_any = TimerManager::getManager()->
         getTimer("tbox::AsyncCommStage::mpi_wait_any");
      t_wait_some = TimerManager::getManager()->
         getTimer("tbox::AsyncCommStage::mpi_wait_some");
      ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }
   return;
}





/*
****************************************************************
*                                                              *
* Make sure the remaining groups are not in the middle         *
* of a communication.  Deallocate allocated groups.            *
*                                                              *
****************************************************************
*/
AsyncCommStage::~AsyncCommStage()
{
   int igroup;
   for ( igroup=0; igroup<d_n_group; ++igroup ) {
      if ( d_group[igroup] != NULL ) {
         // Found an undeallocated group.  Deallocate it.
         if ( ! d_group[igroup]->isDone() ) {
            TBOX_ERROR("Destructing a stage while some groups\n"
                       <<"have pending communication leads to\n"
                       <<"abandoned MPI messages.  Group " << igroup << "\n"
                       <<"is not yet done.");
         }
         delete d_group[igroup];
         d_group[igroup] = NULL;
      }
   }
   return;
}





/*
****************************************************************
*                                                              *
* Allocate a new communication group.  Make sure there         *
* is enough space in the arrays, reallocating space if         *
* needed.  Allocate the new group and set the internal         *
* arrays for it.                                               *
*                                                              *
****************************************************************
*/
AsyncCommGroup *AsyncCommStage::allocateCommGroup(
   int nchild,
   RelaunchableJob *handle )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( nchild < 1 ) {
      TBOX_ERROR("Group nchild must be at least one.\n");
   }
#endif

   /*
    * Allocate space at the end of the current arrays
    * for the group and its needed requests.
    */

   const int new_n_group = d_n_group + 1;
   if ( new_n_group > d_group.size() ) {
      // Additional storage needed for group.
      int new_size = d_n_group == 0 ? 1 : 2*d_n_group;
      d_group.resizeArray(new_size);
      d_group_to_req.resizeArray(new_size);
   }

   d_group[d_n_group] = new StagedGroup( nchild,
                                         this,
                                         d_n_group,
                                         handle );
   d_group_to_req[d_n_group] = d_n_req;

   const int new_n_req = d_n_req + nchild;
   if ( new_n_req > d_req.size() ) {
      int new_size = new_n_req > 2*d_n_req ? new_n_req : 2*d_n_req;
      d_req.resizeArray( new_size );
      d_req_to_group.resizeArray( new_size );
#ifdef DEBUG_CHECK_ASSERTIONS
      for ( int i=new_n_req; i<d_req.size(); ++i ) {
         d_req[i] = MPI_REQUEST_NULL;
         d_req_to_group[i] = -1;
      }
#endif
   }
   for ( int i=d_n_req; i<d_n_req+nchild; ++i ) {
      d_req[i] = MPI_REQUEST_NULL;
      d_req_to_group[i] = d_n_group;
   }

   ++d_n_group;
   d_n_req += nchild;

   return d_group[d_n_group-1];
}





/*
****************************************************************
*                                                              *
* Implementations for AsyncCommStage::StagedGroup.             *
*                                                              *
****************************************************************
*/




AsyncCommStage::StagedGroup::StagedGroup(
   const int nchild,
   AsyncCommStage *stage,
   const int index,
   RelaunchableJob *handle )
   : AsyncCommGroup(nchild),
     d_stage(stage),
     d_index(index),
     d_handle(handle)
{
   return;
}




AsyncCommStage::StagedGroup::~StagedGroup()
{
   if ( !isDone() ) {
      TBOX_ERROR("Cannot deallocate a group with pending communications.\n"
                 <<"It would corrupt message passing algorithms.\n");
   }
   d_stage->destageGroup(this);
   d_stage = NULL;
   return;
}




SAMRAI_MPI::request *AsyncCommStage::StagedGroup::getRequestPointer() const
{
   if ( d_stage == NULL ) {
      TBOX_ERROR("AssyncCommStage::StagedGroup::getRequestPointer():\n"
                 <<"Empty stage encountered!\n"
                 <<"This probably means that the stage that allocated\n"
                 <<"your AsyncCommGroup has been deallocated.\n"
                 <<"(and your AssyncCommGroup deallocated too!).\n"
                 <<"It is an error to deallocate a stage and still\n"
                 <<"use the group it allocated.\n");

   }
   return d_stage->lookupRequestPointer(d_index);
}





/*
******************************************************************
*                                                                *
* Advance one or more communication groups by using MPI_Waitsome *
* to complete one or more requests.                              *
*                                                                *
* Get one or more completed requests and check their             *
* communication groups to see if any group finished              *
* a communication.  If at least one group finished               *
* its communication, return the relaunchable jobs                *
* registered with those that finished.  If no group              *
* has finished, repeat until at lease one has.                   *
*                                                                *
******************************************************************
*/
int AsyncCommStage::advanceSome(
   Array<RelaunchableJob*> &completed )
{
   if ( d_n_group == 0 ) {
      // Short cut for an empty stage.
      return 0;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for ( int i=d_n_req; i<d_req.size(); ++i ) {
      if( d_req[i] != MPI_REQUEST_NULL )
         TBOX_WARNING("non-null request above d_n_reg.");
   }
#endif

   Array<int> index( d_n_req );
   Array<SAMRAI_MPI::status> mpi_stat( d_n_req );

   if ( completed.size() < d_n_group ) {
      completed.setNull();
      completed.resizeArray(d_n_group);
   }

   int n_group_completed = 0;
   int n_req_completed = 0;

   do {

      int errf;
      t_wait_some->start();
      errf = MPI_Waitsome( d_n_req,
                           d_req.getPointer(),
                           &n_req_completed,
                           index.getPointer(),
                           mpi_stat.getPointer() );
      t_wait_some->stop();
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( n_req_completed <= 0 ) {
         /*
           Undocumented feature of some MPI_Waitsome implementations:
           MPI_Waitsome sets n_req_completed to a negative number
           if all the input requests are MPI_REQUEST_NULL.
         */
         for ( int i=0; i<d_n_req; ++i ) {
            if( d_req[i] != MPI_REQUEST_NULL ) {
               TBOX_ERROR("Library error in AsyncCommStage::advanceSome:\n"
                          <<"errf = " << errf << '\n'
                          <<"MPI_SUCCESS = " << MPI_SUCCESS << '\n'
                          <<"MPI_ERR_IN_STATUS = " << MPI_ERR_IN_STATUS << '\n'
                          <<"MPI_REQUEST_NULL = " << MPI_REQUEST_NULL << '\n'
                          <<"d_n_req = " << d_n_req << '\n'
                          <<"d_req.size() = " << d_req.size() << '\n'
                          <<"n_req_completed = " << n_req_completed << '\n'
                          <<"i = " << i << '\n'
                          );
            }
         }
         for ( int i=d_n_req; i<d_req.size(); ++i ) {
            if( d_req[i] != MPI_REQUEST_NULL )
               TBOX_WARNING("non-null request above d_n_reg.");
         }
      }
      if ( n_req_completed == 0 ) {
         TBOX_ASSERT( isDone() );
      }
#endif
      if ( errf != MPI_SUCCESS ) {
         TBOX_ERROR("Error in MPI_Waitsome call.\n"
                    <<"Error-in-status is "
                    << (errf == MPI_ERR_IN_STATUS)
                    << '\n');
      }

      int i, iout, igroup;
      int n_check_group = 0; // Number of groups to check.
      for ( iout=0; iout<n_req_completed; ++iout ) {
         /*
          * Change index from request index to group index.
          * If the group index is not a duplicate, add it to
          * the list of groups to check (which is actually
          * the same list) and increase n_check_group.
          */
         index[iout] = d_req_to_group[ index[iout] ];
         for ( i=0; i<n_check_group; ++i ) if ( index[i] == index[iout] ) break;
         if ( i == n_check_group ) {
            index[n_check_group++] = index[iout];
         }
      }

      /*
       * Check the groups whose requests completed and count up the
       * groups that completed all their communication tasks.
       */
      for ( igroup=0; igroup<n_check_group; ++igroup ) {
         StagedGroup &group = *d_group[ index[igroup] ];
         bool group_done = group.checkOperation();
         if ( group_done ) {
            completed[n_group_completed++] = group.d_handle;
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT( group.isDone() );
#endif
         }
      }

   } while ( n_req_completed > 0 && n_group_completed == 0 );

   return n_group_completed;
}





/*
****************************************************************
*                                                              *
* Advance one communication groups by using MPI_Waitsome       *
* to complete one requests.                                    *
*                                                              *
* Get one completed request and check its                      *
* communication group to see if the group finished             *
* communication.  If it finished its communication,            *
* return the relaunchable job registered with it.              *
* If not, repeat until one group has completed                 *
* communication.                                               *
*                                                              *
****************************************************************
*/
RelaunchableJob *AsyncCommStage::advanceAny()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   for ( int i=d_n_req; i<d_req.size(); ++i ) {
      if( d_req[i] != MPI_REQUEST_NULL )
         TBOX_WARNING("non-null request above d_n_reg.");
   }
#endif
   int ireq = MPI_UNDEFINED;
   int igroup = -1;
   bool group_done;

   do {

      SAMRAI_MPI::status mpi_stat;
      int errf;
      t_wait_any->start();
      errf = MPI_Waitany( d_n_req,
                          d_req.getPointer(),
                          &ireq,
                          &mpi_stat );
      t_wait_any->stop();
      if ( errf != MPI_SUCCESS ) {
         TBOX_ERROR("Error in MPI_Waitany call.\n"
                    <<"Error-in-status is "
                    << (errf == MPI_ERR_IN_STATUS) << '\n'
                    <<"MPI_ERROR value is " << mpi_stat.MPI_ERROR
                    << '\n');
      }

      if ( ireq == MPI_UNDEFINED ) {
         // All input requests are completed even before waiting.
         break;
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( ireq >= 0 && ireq < d_n_req );
#endif

      igroup = d_req_to_group[ireq];

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( igroup >= 0 && igroup < d_n_group );
      TBOX_ASSERT( d_group[igroup] != NULL );
#endif

      StagedGroup *group = d_group[igroup];
      /*
       * Group index igroup had one request completed.
       * See if all of its requests are now completed.
       * If so, run check() to see if the group is done.
       * (The group may not be done because it may initiate
       * follow-up communication tasks.)
       * Exit the do loop when a group is actually done
       * with all of its communication tasks.
       */
      const int init_req = d_group_to_req[group->d_index];
      const int term_req = init_req + group->getNumberOfChildren();
      int i;
      for ( i=init_req; i<term_req; ++i ) {
         if ( d_req[i] != MPI_REQUEST_NULL ) {
            break;
         }
      }
      if ( i == term_req ) {
         /*
          * Group is done with at least its current communication
          * requests.  Follow-up communication may be launched,
          * so check the return value of checkOperation().
          */
         group_done = d_group[igroup]->checkOperation();
      }
      else {
         group_done = false;
         TBOX_ASSERT( ! d_group[igroup]->isDone() );
      }

   } while ( group_done == false );

   return igroup < 0 ? NULL : d_group[igroup]->d_handle;
}



bool AsyncCommStage::isDone() const
{
   /*
    * The stage is done if there is no outstanding MPI requests.
    */
   int ireq;
   for ( ireq=0; ireq<d_n_req; ++ireq ) {
      if ( d_req[ireq] != MPI_REQUEST_NULL ) break;
   }
   return ( ireq == d_n_req );
}



bool AsyncCommStage::numberOfOutstandingRequests() const
{
   int nreq = 0;
   int ireq;
   for ( ireq=0; ireq<d_n_req; ++ireq ) {
      if ( d_req[ireq] != MPI_REQUEST_NULL ) ++nreq;
   }
   return nreq;
}



bool AsyncCommStage::numberOfOutstandingGroups() const
{
   int ngroup = 0;
   int igroup;
   for ( igroup=0; igroup<d_n_group; ++igroup ) {
      if ( d_group[igroup] != NULL && ! d_group[igroup]->isDone() ) ++ngroup;
   }
   return ngroup;
}





/*
****************************************************************
*                                                              *
* Return the request pointer for a communication group         *
* allocated by this object.                                    *
*                                                              *
****************************************************************
*/
SAMRAI_MPI::request *AsyncCommStage::lookupRequestPointer(
   const int igroup ) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( igroup < d_n_group );
   TBOX_ASSERT( d_group[igroup] != NULL );
#endif
   return &d_req[ d_group_to_req[igroup] ];
}





/*
*******************************************************************
*                                                                 *
* Remove mutual reference between group and stage that allocated it. *
* 1. Confirm that the group is currently in the stage (or else it *
*    is an illegal operation).                                    *
* 2. Remove mutual references.                                    *
* 3. Reduce the stage data size, if possible, by skimming unused  *
*    space off the ends of arrays, if possible.                   *
*                                                                 *
*******************************************************************
*/
void AsyncCommStage::destageGroup( StagedGroup *group )
{
   if ( !group->isDone() ) {
      TBOX_ERROR("Cannot clear a group with pending communications.\n"
                 <<"It would corrupt message passing algorithms.\n");
   }

   int igroup = group->d_index;
   if ( d_group[igroup] != group ) {
      // Group was not allocated by this stage.
      TBOX_ERROR("A ArrayCommStage should not clear a\n"
                 <<"AsyncCommGroup object it did not\n"
                 <<"allocate.");
   }

   group->d_stage = NULL;
   d_group[igroup] = NULL;

   /*
    * Reduce d_n_group and d_n_req as much as possible without
    * shifting arrays.  This is only possible if igroup == d_n_group-1.
    */
   while ( d_n_group > 0 && d_group[d_n_group-1] == NULL ) {
      --d_n_group;
   }

   return;
}




/*
***************************************************************************
*                                                                         *
* Release static timers.  To be called by shutdown registry to make sure  *
* memory for timers does not leak.                                        *
*                                                                         *
***************************************************************************
*/
void AsyncCommStage::freeTimers()
{
   t_wait_any.setNull();
   t_wait_some.setNull();
   return;
}




#else


AsyncCommStage::AsyncCommStage()
   : d_n_group(0)
{
   return;
}


AsyncCommStage::~AsyncCommStage()
{
   return;
}



AsyncCommGroup *AsyncCommStage::allocateCommGroup(
   int nchild,
   RelaunchableJob *handle )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( nchild < 1 ) {
      TBOX_ERROR("Group nchild must be at least one.\n");
   }
#endif

   const int new_n_group = d_n_group + 1;

   if ( new_n_group > d_group.size() ) {
      // Additional storage needed for group.
      int new_size = d_n_group == 0 ? 1 : 2*d_n_group;
      d_group.resizeArray(new_size);
   }

   d_group[d_n_group] = new StagedGroup( nchild,
                                         this,
                                         d_n_group,
                                         handle );

   ++d_n_group;

   return d_group[d_n_group-1];
}





RelaunchableJob *AsyncCommStage::advanceAny()
{
   return NULL;
}





int AsyncCommStage::advanceSome(
   Array<RelaunchableJob*> &completed )
{
   completed.setNull();
   return 0;
}



bool AsyncCommStage::isDone() const
{
   return true;
}



bool AsyncCommStage::numberOfOutstandingRequests() const
{
   return 0;
}



bool AsyncCommStage::numberOfOutstandingGroups() const
{
   return 0;
}





void AsyncCommStage::destageGroup( StagedGroup *group )
{
   if ( !group->isDone() ) {
      TBOX_ERROR("Cannot clear a group with pending communications.\n"
                 <<"It would corrupt message passing algorithms.\n");
   }

   int igroup = group->d_index;
   if ( d_group[igroup] != group ) {
      // Group was not allocated by this stage.
      TBOX_ERROR("A ArrayCommStage should not clear a\n"
                 <<"AsyncCommGroup object it did not\n"
                 <<"allocate.");
   }

   group->d_stage = NULL;
   d_group[igroup] = NULL;

   /*
    * Reduce d_n_group and d_n_req as much as possible without
    * shifting arrays.  This is only possible if igroup == d_n_group-1.
    */
   while ( d_n_group > 0 && d_group[d_n_group-1] == NULL ) {
      --d_n_group;
   }

   return;
}





/*
****************************************************************
*                                                              *
* Implementations for AsyncCommStage::StagedGroup.             *
*                                                              *
****************************************************************
*/




AsyncCommStage::StagedGroup::StagedGroup(
   const int nchild,
   AsyncCommStage *stage,
   const int index,
   RelaunchableJob *handle )
   : AsyncCommGroup(nchild),
     d_stage(stage),
     d_index(index),
     d_handle(handle)
{
   return;
}




AsyncCommStage::StagedGroup::~StagedGroup()
{
   if ( !isDone() ) {
      TBOX_ERROR("Cannot deallocate a group with pending communications.\n"
                 <<"It would corrupt message passing algorithms.\n");
   }
   d_stage->destageGroup(this);
   return;
}




SAMRAI_MPI::request *AsyncCommStage::StagedGroup::getRequestPointer() const
{
   return NULL;
}





#endif



}
}
