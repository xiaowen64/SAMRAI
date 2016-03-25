/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/parallel/AsyncCommGroup.C $
 * Package:     SAMRAI toolbox
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2037 $
 * Modified:    $LastChangedDate: 2008-03-05 15:54:45 -0800 (Wed, 05 Mar 2008) $
 * Description: All-to-one and one-to-all communication using a tree.
 */

#include "tbox/AsyncCommGroup.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/AsyncCommGroup.I"
#endif

#include "tbox/PIO.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/MathUtilities.h"
#include STL_SSTREAM_HEADER_FILE

#ifdef OSTRINGSTREAM_TYPE_IS_BROKEN
#ifdef OSTRSTREAM_TYPE_IS_BROKEN
#error "Neither ostringstream nor ostrstream works"
#else
typedef ostringstream ostrstream
#endif
#endif


using namespace std;
namespace SAMRAI {
   namespace tbox {


#ifdef HAVE_MPI

/*
 * This class uses a non-deterministic algorithm, which can be
 * very hard to debug.  To help debugging, we keep some special
 * debugging code that is activated when AsyncCommGroup_DEBUG_OUTPUT
 * is defined.
 */
// #define AsyncCommGroup_DEBUG_OUTPUT

Pointer<Timer> AsyncCommGroup::t_reduce_data;
Pointer<Timer> AsyncCommGroup::t_wait_all;


/*
***********************************************************************
*                                                                     *
* Construct a simple communication group that does not work           *
* with a communication stage.  All parameters are set to reasonable   *
* defaults or, if appropriate, invalid values.                        *
*                                                                     *
***********************************************************************
*/
AsyncCommGroup::AsyncCommGroup( const int nchild )
   :
   d_nchild(nchild),
   d_idx(-1),
   d_root_idx(-1),
   d_parent_rank(-1),
   d_child_data(new ChildData[nchild]),
   d_branch_size_totl(-1),
   d_base_op(undefined),
   d_next_task_op(none),
   d_external_buf(NULL),
   d_external_size(0),
   d_internal_buf(0),
   d_internal_requests(NULL),
   d_mpi_tag(-1),
   d_mpi_communicator(MPI_COMM_WORLD),
   d_use_blocking_send_to_children(false),
   d_use_blocking_send_to_parent(true)
#ifdef DEBUG_CHECK_ASSERTIONS
   , d_group_ranks(0)
#endif
{
   if ( t_reduce_data.isNull() ) {
      TBOX_ASSERT( t_wait_all.isNull() );
      t_reduce_data = TimerManager::getManager()->
         getTimer("tbox::AsyncCommGroup::reduceData()");
      t_wait_all = TimerManager::getManager()->
         getTimer("tbox::AsyncCommGroup::mpi_wait_all");
      ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }
   return;
}


AsyncCommGroup::~AsyncCommGroup()
{
   if ( ! isDone() ) {
      TBOX_ERROR("Deallocating a group while communication is pending\n"
                 <<"leads to lost messages."
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   }
   delete [] d_child_data;
   d_child_data = NULL;
   if ( d_internal_requests != NULL ) {
      delete [] d_internal_requests;
   }
   return;
}





/*
*********************************************************************
*                                                                   *
* Check whether the current (or last) operation has completed.      *
*                                                                   *
*********************************************************************
*/
bool AsyncCommGroup::checkOperation()
{
   switch (d_base_op) {
   case gather: return checkGather();
   case bcast: return checkBcast();
   case min_reduce:
   case max_reduce:
   case sum_reduce: return checkReduce();
   case undefined:
      TBOX_ERROR("There is no current operation to check.\n"
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   default:
      TBOX_ERROR("Library error: attempt to use an operation that\n"
                 <<"has not been written yet"
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   }
   return true;
}


/*
*********************************************************************
*                                                                   *
* Wait for current communication operation to complete.             *
*                                                                   *
* Wait for all requests to come in and call checkOperation()        *
* until all tasks of the communication operation is complete.       *
*                                                                   *
*********************************************************************
*/
void AsyncCommGroup::waitOperation()
{
   SAMRAI_MPI::request * const req = getRequestPointer();
   SAMRAI_MPI::status *mpi_stat = d_next_task_op == none ?
      (SAMRAI_MPI::status*)NULL : new SAMRAI_MPI::status[d_nchild];

   while ( d_next_task_op != none ) {

      t_wait_all->start();
      int errf = MPI_Waitall( d_nchild,
                              req,
                              mpi_stat );
      t_wait_all->stop();

      if ( errf != MPI_SUCCESS ) {
         TBOX_ERROR("Error in MPI_Waitall call.\n"
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
      }

      checkOperation();

   }

   if ( mpi_stat != NULL ) {
      delete [] mpi_stat;
   }

   return;
}





/*
************************************************************************
*                                                                      *
* Set internal parameters for performing the broadcast                 *
* and call checkBcast to perform the communication.                    *
*                                                                      *
************************************************************************
*/
bool AsyncCommGroup::beginBcast( int *buffer, int size )
{
   if ( d_next_task_op != none ) {
      TBOX_ERROR("Cannot begin communication while another is in progress."
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   checkMPIParams();
#endif
   d_external_buf = buffer;
   d_external_size = size;
   d_base_op = bcast;
   d_next_task_op = recv_start;
   return checkBcast();
}


/*
************************************************************************
*                                                                      *
* Broadcast is an all-to-one operation, so we receive from the         *
* parent process and send to the children processes.                   *
*                                                                      *
************************************************************************
*/
bool AsyncCommGroup::checkBcast()
{
   if ( d_base_op != bcast ) {
      TBOX_ERROR("Cannot check nonexistent broadcast operation."
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   }
   SAMRAI_MPI::request * const req = getRequestPointer();
   int ic;
   int flag = 0;

   switch (d_next_task_op) {

   case none:
      break;

   case recv_start:
      if ( d_parent_rank > -1 ) {
         d_mpi_err = MPI_Irecv( d_external_buf,
                                d_external_size,
                                MPI_INT,
                                d_parent_rank,
                                d_mpi_tag,
                                d_mpi_communicator,
                                &req[0] );
         if ( d_mpi_err != MPI_SUCCESS ) {
            TBOX_ERROR("Error in MPI_Irecv."
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
         }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
         tbox::plog << "tag-" << d_mpi_tag
                    << " expecting " << d_external_size
                    << " from " << d_parent_rank
                    << " in checkBcast"
                    << endl;
#endif
      }

   case recv_check:
      if ( req[0] != MPI_REQUEST_NULL ) {
         resetStatus();
         d_mpi_err = MPI_Test( &req[0], &flag, &d_mpi_status );
         if ( d_mpi_err != MPI_SUCCESS ) {
            TBOX_ERROR("Error in MPI_Test.\n"
                       <<"Error-in-status is "
                       << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                       <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                       << '\n'
                       <<"mpi_communicator = " << d_mpi_communicator
                       << "mpi_tag = " << d_mpi_tag);
         }
         if ( flag == 1 ) {
#ifdef DEBUG_CHECK_ASSERTIONS
            int count = -1;
            d_mpi_err = MPI_Get_count( &d_mpi_status, MPI_INT, &count );
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in MPI_Get_count.\n"
                          <<"Error-in-status is "
                          << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                          <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                          << '\n'
                          <<"mpi_communicator = " << d_mpi_communicator
                          << "mpi_tag = " << d_mpi_tag);
            }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            tbox::plog << "tag-" << d_mpi_tag
                       << " received " << count
                       << " from " << d_mpi_status.MPI_SOURCE
                       << " in checkBcast"
                       << endl;
#endif
            TBOX_ASSERT( count <= d_external_size );
            TBOX_ASSERT( d_mpi_status.MPI_TAG == d_mpi_tag );
            TBOX_ASSERT( d_mpi_status.MPI_SOURCE == d_parent_rank );
            TBOX_ASSERT( req[0] == MPI_REQUEST_NULL );
#endif
         }
         else {
            d_next_task_op = recv_check;
            break;
         }
      }

   case send_start:
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( d_child_data[ic].rank >= 0 ) {
            if ( d_use_blocking_send_to_children ) {
               d_mpi_err = MPI_Send( d_external_buf,
                                     d_external_size,
                                     MPI_INT,
                                     d_child_data[ic].rank,
                                     d_mpi_tag,
                                     d_mpi_communicator );
            }
            else {
               d_mpi_err = MPI_Isend( d_external_buf,
                                      d_external_size,
                                      MPI_INT,
                                      d_child_data[ic].rank,
                                      d_mpi_tag,
                                      d_mpi_communicator,
                                      &req[ic] );
            }
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in send."
                          <<"mpi_communicator = " << d_mpi_communicator
                          << "mpi_tag = " << d_mpi_tag);
            }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            tbox::plog << "tag-" << d_mpi_tag
                       << " sending " << d_external_size
                       << " to " << d_child_data[ic].rank
                       << " in checkBcast"
                       << endl;
#endif
         }
      }

   case send_check:
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( req[ic] != MPI_REQUEST_NULL ) {
            resetStatus();
            d_mpi_err = MPI_Test( &req[ic], &flag, &d_mpi_status );
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in MPI_Test.\n"
                          <<"Error-in-status is "
                          << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                          <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                          << '\n'
                          <<"mpi_communicator = " << d_mpi_communicator
                          << "mpi_tag = " << d_mpi_tag);
            }
#ifdef DEBUG_CHECK_ASSERTIONS
            if ( req[ic] == MPI_REQUEST_NULL ) {
               int count = -1;
               MPI_Get_count( &d_mpi_status, MPI_INT, &count );
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               tbox::plog << "tag-" << d_mpi_tag
                          << " sent unknown size (MPI convention)"
                          << " to " << d_child_data[ic].rank
                          << " in checkBcast"
                          << endl;
#endif
            }
#endif
         }
      }
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( req[ic] != MPI_REQUEST_NULL ) {
            break;
         }
      }
      if ( ic < d_nchild ) {
         d_next_task_op = send_check;
         break;
      }
      d_next_task_op = none;
      break;

   default:
      TBOX_ERROR("checkBcast is incompatible with current state."
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   }

   if ( d_parent_rank == -1 ) {
      TBOX_ASSERT( d_next_task_op != recv_check );
   }

   return d_next_task_op == none;
}






/*
************************************************************************
*                                                                      *
* Allocate enough memory internally to store all descendent data.      *
* Place local process's contribution in the internall buffer.          *
* Call checkGather to obtain data from descendants and send            *
* to parent..                                                          *
*                                                                      *
************************************************************************
*/
bool AsyncCommGroup::beginGather( int *buffer, int size )
{
   /*
    * The internal buffer size is d_branch_size_totl+1 times the
    * message size.  There is one data block for each descendent,
    * plus one for the local contribution:
    *
    * |<------------------------ internal buffer ---------------------->|
    * |                                                                 |
    * |<--msg_size-->|<--msg_size-->| ... |<--msg_size-->|<--msg_size-->|
    * |              |              | ... |              |              |
    * |  recv from   |  recv from   | ... |  recv from   |    local     |
    * | descendant 0 | descendant 1 | ... | dsndt (nb-1) | contribution |
    *
    * (nb = d_branch_size_totl)
    *
    * The message size is the data buffer size, plus one integer
    * describing the index of the process contributing the data.
    */

   if ( d_next_task_op != none ) {
      TBOX_ERROR("Cannot begin communication while another is in progress."
                 <<"mpi_communicator = " << d_mpi_communicator
                 << "mpi_tag = " << d_mpi_tag);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   checkMPIParams();
#endif

   d_external_buf = buffer;
   d_external_size = size;

   /*
    * Allocate enough space for data from this position plus all
    * descendant positions.  Each position contributes its external
    * data plus some data to help sort the final gathered data.
    */
   int per_proc_msg_size = (1+d_external_size);
   d_internal_buf.setNull();
   d_internal_buf.resizeArray( (d_branch_size_totl+1)*per_proc_msg_size );

   /*
    * Add our contribution to the gathered data.
    */
   int *ptr = d_internal_buf.getPointer()
      + (d_branch_size_totl)*per_proc_msg_size;
   *(ptr++) = d_idx;
   int i;
   for ( i=0; i<d_external_size; ++i ) {
      ptr[i] = d_external_buf[i];
   }

   d_base_op = gather;
   d_next_task_op = recv_start;

   return checkGather();
}




/*
************************************************************************
*                                                                      *
* Gather is an all-to-one operation, so we receive from the            *
* children processes and send to the parent process.                   *
*                                                                      *
************************************************************************
*/
bool AsyncCommGroup::checkGather()
{
   if ( d_base_op != gather ) {
      TBOX_ERROR("Cannot check nonexistent gather operation\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }

   SAMRAI_MPI::request * const req = getRequestPointer();
   int per_proc_msg_size = (1+d_external_size);

   int i;
   int ic, older_sibling_size;
   int flag = 0;

   switch (d_next_task_op) {

   case none:
      break;

   case recv_start:
      older_sibling_size = 0;
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( d_child_data[ic].rank >= 0 ) {
            d_mpi_err = MPI_Irecv( d_internal_buf.getPointer() +
                                      per_proc_msg_size*older_sibling_size,
                                   d_child_data[ic].size*per_proc_msg_size,
                                   MPI_INT,
                                   d_child_data[ic].rank,
                                   d_mpi_tag,
                                   d_mpi_communicator,
                                   &req[ic] );
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in MPI_Irecv.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
            }
            older_sibling_size += d_child_data[ic].size;
         }
      }

   case recv_check:

      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( req[ic] != MPI_REQUEST_NULL ) {
            resetStatus();
            d_mpi_err = MPI_Test( &req[ic], &flag, &d_mpi_status );
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in MPI_Test.\n"
                          <<"Error-in-status is "
                          << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                          <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                          << '\n'
                          <<"mpi_communicator = " << d_mpi_communicator << '\n'
                          <<"mpi_tag = " << d_mpi_tag << '\n');
            }
#ifdef DEBUG_CHECK_ASSERTIONS
            if ( flag == 1 ) {
               TBOX_ASSERT( d_mpi_status.MPI_TAG == d_mpi_tag );
               TBOX_ASSERT( d_mpi_status.MPI_SOURCE == d_child_data[ic].rank );
               TBOX_ASSERT( req[ic] == MPI_REQUEST_NULL );
               int count = -1;
               MPI_Get_count( &d_mpi_status, MPI_INT, &count );
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               tbox::plog << "tag-" << d_mpi_tag
                          << " received " << count
                          << " from " << d_mpi_status.MPI_SOURCE
                          << " in checkGather"
                          << endl;
#endif
               if( count > d_child_data[ic].size*per_proc_msg_size ) {
                  TBOX_ERROR("Message size bigger than expected from proc "
                             << d_child_data[ic].rank << "\n"
                             <<"Expect "
                             << d_child_data[ic].size*per_proc_msg_size
                             << "\n"
                             <<"Actual " << count << '\n'
                             <<"mpi_communicator = " << d_mpi_communicator << '\n'
                             <<"mpi_tag = " << d_mpi_tag << '\n');
               }
            }
#endif
         }
      }
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( req[ic] != MPI_REQUEST_NULL ) {
            break;
         }
      }
      if ( ic < d_nchild ) {
         d_next_task_op = recv_check;
         break;
      }

      if ( d_parent_rank < 0 ) {
         /*
          * Transfer the internal buffer into the external buffer,
          * unshuffling data in the process.
          */
         int n;
         for ( n=0; n<d_group_size; ++n ) {
            int *ptr = d_internal_buf.getPointer()
                     + per_proc_msg_size*n;
            const int source_idx = *(ptr++);
            if( source_idx < 0 && source_idx >= d_group_size ) {
               TBOX_ERROR("Gathered data has out of range index "
                          << source_idx << '\n'
                          <<"mpi_communicator = " << d_mpi_communicator << '\n'
                          <<"mpi_tag = " << d_mpi_tag << '\n');
            }
            int *dest_buf = d_external_buf
                          + d_external_size*(source_idx);
            for ( i=0; i<d_external_size; ++i, ++ptr ) {
               dest_buf[i] = *ptr;
            }
         }
      }

   case send_start:
      if ( d_parent_rank >= 0 ) {
         if ( d_use_blocking_send_to_parent ) {
            d_mpi_err = MPI_Send( d_internal_buf.getPointer(),
                                  d_internal_buf.size(),
                                  MPI_INT,
                                  d_parent_rank,
                                  d_mpi_tag,
                                  d_mpi_communicator );
         }
         else {
            d_mpi_err = MPI_Isend( d_internal_buf.getPointer(),
                                   d_internal_buf.size(),
                                   MPI_INT,
                                   d_parent_rank,
                                   d_mpi_tag,
                                   d_mpi_communicator,
                                   &req[0] );
         }
         if ( d_mpi_err != MPI_SUCCESS ) {
            TBOX_ERROR("Error in send.\n"
                       <<"mpi_communicator = " << d_mpi_communicator << '\n'
                       <<"mpi_tag = " << d_mpi_tag << '\n');
         }
      }

   case send_check:
      if ( req[0] != MPI_REQUEST_NULL ) {
         resetStatus();
         d_mpi_err = MPI_Test( &req[0], &flag, &d_mpi_status );
         if ( d_mpi_err != MPI_SUCCESS ) {
            TBOX_ERROR("Error in MPI_Test.\n"
                       <<"Error-in-status is "
                       << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                       <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                       << '\n'
                       <<"mpi_communicator = " << d_mpi_communicator << '\n'
                       <<"mpi_tag = " << d_mpi_tag << '\n');
         }
      }
      if ( req[0] != MPI_REQUEST_NULL ) {
         d_next_task_op = send_check;
      }
      else {
         d_next_task_op = none;
      }
      break;

   default:
      TBOX_ERROR("checkGather is incompatible with current state.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }

   return d_next_task_op == none;
}





/*
**********************************************************************
*                                                                    *
* Set flag indicating a sum reduce and call the generic reduce method *
* to do the actual work.                                             *
*                                                                    *
**********************************************************************
*/
bool AsyncCommGroup::beginSumReduce( int *buffer, int size )
{
   if ( d_next_task_op != none ) {
      TBOX_ERROR("Cannot begin communication while another is in progress.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   checkMPIParams();
#endif
   d_base_op = sum_reduce;
   d_external_buf = buffer;
   d_external_size = size;
   return beginReduce();
}


/*
**********************************************************************
*                                                                    *
* The fact that the current reduction is a sum reduction is recorded *
* in the d_base_op parameter.  We simply use the generic method for  *
* checking reduction.                                                *
*                                                                    *
**********************************************************************
*/
bool AsyncCommGroup::checkSumReduce()
{
   return checkReduce();
}






/*
************************************************************************
*                                                                      *
* Allocate enough memory internally to store all children data.        *
* Place local process's contribution in the internall buffer.          *
* Call checkReduce to obtain data from children, reduce the            *
* data and send result to parent..                                     *
*                                                                      *
************************************************************************
*/
bool AsyncCommGroup::beginReduce()
{
   int msg_size = d_external_size;
   /*
    * For reducing data, nc = number of actual children.  nc <= d_nchild.
    *
    * The internal buffer stores up to nc+1 times the message size:
    *
    * |<------------------------ internal buffer ---------------------->|
    * |                                                                 |
    * |<--msg_size-->|<--msg_size-->| ... |<--msg_size-->|<--msg_size-->|
    * |              |              | ... |              |              |
    * |  recv from   |  recv from   | ... |  recv from   |   send to    |
    * |   child 0    |   child 1    | ... | child (nc-1) |    parent    |
    * |   (if any)   |   (if any)   | ... |   (if any)   |   (if any)   |
    *
    * If a data block is not needed, it is not created and the
    * following blocks shift over.  For non-root processes,
    * the reduced data is placed in the "send to parent" section
    * so it can be passed up the tree.
    * For the root process, reduced data is placed directly
    * into d_external_buf.
    */

   /*
    * Compute the number of actual children, nc.  Note that
    * d_nchild is just the upper limit to the number of actual
    * children.
    */
   const int oldest_pos = toOldest(toPosition(d_idx));
   int limit_pos = toYoungest(toPosition(d_idx)) + 1;
   limit_pos = tbox::MathUtilities<int>::Min(limit_pos, d_group_size);
   const int n_children = (limit_pos > oldest_pos ? limit_pos - oldest_pos : 0);

   d_internal_buf.setNull();
   d_internal_buf.resizeArray( msg_size *
                               ( n_children + (d_parent_rank > -1) ) );

   if ( d_parent_rank > -1 ) {
      int *ptr = d_internal_buf.getPointer() + d_internal_buf.size() - msg_size;
      for ( int i=0; i<d_external_size; ++i ) {
         ptr[i] = d_external_buf[i];
      }
   }

   d_next_task_op = recv_start;

   return checkReduce();
}


/*
************************************************************************
*                                                                      *
* Reduction is an all-to-one operation, so we receive from             *
* the children and send to the parent.  The transfer data              *
* toward the root process is the same for all types of                 *
* reductions.  After each child's data is received, the                *
* specific data reduction arithmetic is done locally.                  *
* The result of the local reduction is sent to the parent.             *
*                                                                      *
************************************************************************
*/
bool AsyncCommGroup::checkReduce()
{
   if ( ! ( d_base_op == max_reduce ||
            d_base_op == min_reduce ||
            d_base_op == sum_reduce ) ) {
      TBOX_ERROR("Cannot check nonexistent reduce operation.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }

   SAMRAI_MPI::request * const req = getRequestPointer();
   int msg_size = d_external_size;

   int ic;
   int flag = 0;
   switch (d_next_task_op) {

   case none:
      break;

   case recv_start:
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( d_child_data[ic].rank >= 0 ) {
            d_mpi_err = MPI_Irecv( d_internal_buf.getPointer() + ic*msg_size,
                                   msg_size,
                                   MPI_INT,
                                   d_child_data[ic].rank,
                                   d_mpi_tag,
                                   d_mpi_communicator,
                                   &req[ic] );
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in MPI_Irecv.\n"
                          <<"mpi_communicator = " << d_mpi_communicator << '\n'
                          <<"mpi_tag = " << d_mpi_tag << '\n');
            }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            tbox::plog << "tag-" << d_mpi_tag
                       << " expecting " << msg_size
                       << " from " << d_child_data[ic].rank
                       << " in checkReduce"
                       << endl;
#endif
         }
      }

   case recv_check:

      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( req[ic] != MPI_REQUEST_NULL ) {
            resetStatus();
            d_mpi_err = MPI_Test( &req[ic], &flag, &d_mpi_status );
            if ( d_mpi_err != MPI_SUCCESS ) {
               TBOX_ERROR("Error in MPI_Test.\n"
                          <<"Error-in-status is "
                          << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                          <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                          << '\n'
                          <<"mpi_communicator = " << d_mpi_communicator << '\n'
                          <<"mpi_tag = " << d_mpi_tag << '\n');
            }
            if ( flag == 1 ) {
#ifdef DEBUG_CHECK_ASSERTIONS
               TBOX_ASSERT( d_mpi_status.MPI_TAG == d_mpi_tag );
               TBOX_ASSERT( d_mpi_status.MPI_SOURCE == d_child_data[ic].rank );
               TBOX_ASSERT( req[ic] == MPI_REQUEST_NULL );
               int count = -1;
               MPI_Get_count( &d_mpi_status, MPI_INT, &count );
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               tbox::plog << "tag-" << d_mpi_tag
                          << " received " << count
                          << " from " << d_mpi_status.MPI_SOURCE
                          << " in checkReduce"
                          << endl;
#endif
               if( count != msg_size ) {
                  TBOX_ERROR("Did not get the expected message size from proc "
                             << d_child_data[ic].rank << "\n"
                             <<"Expect " << msg_size << "\n"
                             <<"Actual " << count << '\n'
                             <<"mpi_communicator = " << d_mpi_communicator << '\n'
                             <<"mpi_tag = " << d_mpi_tag << '\n');
               }
#endif
            }
         } 
      } 
      for ( ic=0; ic<d_nchild; ++ic ) {
         if ( req[ic] != MPI_REQUEST_NULL ) {
            break;
         }
      }
      if ( ic < d_nchild ) {
         d_next_task_op = recv_check;
         break;
      }


      {
         int *local_data = d_parent_rank < 0 ? d_external_buf :
            d_internal_buf.getPointer() + d_internal_buf.size() - msg_size;
         for ( ic=0; ic<d_nchild; ++ic ) {
            if ( d_child_data[ic].rank > -1 ) {
               int *child_data = d_internal_buf.getPointer() + ic*msg_size;
               t_reduce_data->start();
               reduceData( local_data, child_data );
               t_reduce_data->stop();
            }
         }
      }

   case send_start:
      if ( d_parent_rank >= 0 ) {
         int *ptr = d_internal_buf.getPointer() +
            d_internal_buf.size() - msg_size;
         if ( d_use_blocking_send_to_parent ) {
            d_mpi_err = MPI_Send( ptr,
                                  msg_size,
                                  MPI_INT,
                                  d_parent_rank,
                                  d_mpi_tag,
                                  d_mpi_communicator );
         }
         else {
            d_mpi_err = MPI_Isend( ptr,
                                   msg_size,
                                   MPI_INT,
                                   d_parent_rank,
                                   d_mpi_tag,
                                   d_mpi_communicator,
                                   &req[0] );
         }
         if ( d_mpi_err != MPI_SUCCESS ) {
            TBOX_ERROR("Error in send.\n"
                       <<"mpi_communicator = " << d_mpi_communicator << '\n'
                       <<"mpi_tag = " << d_mpi_tag << '\n');
         }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
         tbox::plog << "tag-" << d_mpi_tag
                    << " sending " << msg_size
                    << " to " << d_parent_rank
                    << " in checkReduce"
                    << endl;
#endif
      }

   case send_check:
      if ( req[0] != MPI_REQUEST_NULL ) {
         resetStatus();
         d_mpi_err = MPI_Test( &req[0], &flag, &d_mpi_status );
         if ( d_mpi_err != MPI_SUCCESS ) {
            TBOX_ERROR("Error in MPI_Test.\n"
                       <<"Error-in-status is "
                       << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                       <<"MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                       << '\n'
                       <<"mpi_communicator = " << d_mpi_communicator << '\n'
                       <<"mpi_tag = " << d_mpi_tag << '\n');
         }
      }
      if ( req[0] != MPI_REQUEST_NULL ) {
#ifdef DEBUG_CHECK_ASSERTIONS
         int count = -1;
         MPI_Get_count( &d_mpi_status, MPI_INT, &count );
#ifdef AsyncCommGroup_DEBUG_OUTPUT
         tbox::plog << "tag-" << d_mpi_tag
                    << " sent " << count
                    << " to " << d_parent_rank
                    << " in checkReduce"
                    << endl;
#endif
#endif
         d_next_task_op = send_check;
      }
      else {
         d_next_task_op = none;
      }
      break;

   default:
      TBOX_ERROR("checkReduce is incompatible with current state.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }

   if ( d_parent_rank == -1 ) {
      TBOX_ASSERT( d_next_task_op != send_check );
   }

   /*
    * Receive from parent.
    */
   return d_next_task_op == none;
}


void AsyncCommGroup::setGroupAndRootRank( const Array<int> &group,
                                          const int root_rank )
{
   int i;
   for ( i=0; i<group.size(); ++i ) {
      if ( group[i] == root_rank ) {
         break;
      }
   }
   if ( i == group.size() ) {
      TBOX_ERROR("New root " << root_rank << " is not in the group.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }
   setGroupAndRootIndex( group, i );
   return;
}


void AsyncCommGroup::setGroupAndRootIndex( const Array<int> &group,
                                           const int root_index )
{
   if ( d_next_task_op != none ) {
      TBOX_ERROR("Changing group while a communication is occuring can\n"
                 <<"corrupt data and so is disallowed.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }

   const int rank = SAMRAI_MPI::getRank();

   // Set the index for local and root processes.
   d_group_size = group.size();
   d_idx = -1;
   for ( int i=0; i<d_group_size; ++i ) {
      if ( group[i] == rank ) {
         d_idx = i;
         break;
      }
   }
   d_root_idx = root_index;

#ifdef DEBUG_CHECK_ASSERTIONS
   // Set d_group_ranks and do some sanity checks.
   int comm_size;
   d_mpi_err = MPI_Comm_size( d_mpi_communicator, &comm_size );
   if ( d_mpi_err != MPI_SUCCESS ) {
      TBOX_ERROR("Error in MPI_Comm_size.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }
   if ( comm_size < d_group_size ) {
      TBOX_ERROR("Group size must not be greater than the size of\n"
                 <<"the MPI communicator group.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');

   }
   TBOX_ASSERT( d_group_size > 0 );
   if ( d_idx == -1 ) {
      TBOX_ERROR("The local process MUST be in the communication group.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }
   d_group_ranks.setNull();
   d_group_ranks.resizeArray( d_group_size );
   int dup = 0;
   for ( int i=0; i<d_group_size; ++i ) {
      if ( group[i] < 0 || group[i] >= comm_size ) {
         TBOX_ERROR("Rank " << group[i] << " is not in the current\n"
                    <<"MPI communicator.\n"
                    <<"mpi_communicator = " << d_mpi_communicator << '\n'
                    <<"mpi_tag = " << d_mpi_tag << '\n');
      }
      if ( group[i] == rank ) ++dup;
      d_group_ranks[i] = group[i];
   }
   if ( dup != 1 ) {
      TBOX_ERROR("The local process must appear exactly once in the group.\n"
                 <<"It appeared " << dup << " times.\n"
                 <<"mpi_communicator = " << d_mpi_communicator << '\n'
                 <<"mpi_tag = " << d_mpi_tag << '\n');
   }
#endif

   computeDependentData( group );

   return;
}




/*
***************************************************************************
*                                                                         *
* Compute data that is dependent on the group, root and local processes.  *
*                                                                         *
***************************************************************************
*/
void AsyncCommGroup::computeDependentData( const Array<int> &group )
{
   /*
    * Compute number of descendants in each child branch and in all branches.
    * To find the number of descendants in each branch find the oldest
    * and youngest descendants for each generation.  Add up contribution
    * of all descendant generations.
    */
   d_branch_size_totl = 0;
   int ic;
   for ( ic=0; ic<d_nchild; ++ic ) {
      d_child_data[ic].size = 0;
      int pos_of_child_ic = toChildPosition( toPosition(d_idx), ic );
      int oldest = pos_of_child_ic;
      int yngest = pos_of_child_ic;
      while ( oldest < d_group_size ) {
         d_child_data[ic].size +=
            ( yngest < d_group_size ? yngest : d_group_size-1 ) - oldest + 1;
         oldest = toOldest(oldest);
         if ( yngest < d_group_size ) yngest = toYoungest(yngest);
      }
      d_branch_size_totl += d_child_data[ic].size;
   }

   /*
    * Find ranks of parent and children using specific arithmetic
    * relationships between parent and children positions.
    */
   const int pos = toPosition(d_idx);
   if ( pos > 0 ) {
      d_parent_rank = group[ toIndex((pos-1)/d_nchild) ];
   }
   else d_parent_rank = -1;
   for ( ic=0; ic<d_nchild; ++ic ) {
      const int pos_of_child_ic = toChildPosition( pos, ic );
      if ( pos_of_child_ic < d_group_size ) {
         d_child_data[ic].rank = group[ toIndex(pos_of_child_ic) ];
      }
      else d_child_data[ic].rank = -1;
   }

   return;
}





SAMRAI_MPI::request *AsyncCommGroup::getRequestPointer() const
{
   if ( d_internal_requests == NULL ) {
      d_internal_requests = new SAMRAI_MPI::request [d_nchild];
   }
   return d_internal_requests;
}





void AsyncCommGroup::checkMPIParams()
{
   if ( d_mpi_tag < 0 ) {
      TBOX_ERROR("AsyncCommGroup: Invalid MPI tag value "
                 << d_mpi_tag << "\nUse setMPITag() to set it.");
   }
}




void AsyncCommGroup::logCurrentState( ostream &co ) const
{
   SAMRAI_MPI::request * const req = getRequestPointer();
   co << "State=" << 10*d_base_op + d_next_task_op
      << "  tag=" << d_mpi_tag
      << "  communicator=" << d_mpi_communicator
      << "  extern. buff=" << d_external_buf
      << "  size=" << d_external_size
      << "  parent=" << d_parent_rank
      ;
   int i;
   for ( i=0; i<d_nchild; ++i ) {
      co << "  [" << i << "]=" << d_child_data[i].rank
         << " (" << req[i] << ')';
   }
   co << '\n';
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
void AsyncCommGroup::freeTimers()
{
   t_reduce_data.setNull();
   t_wait_all.setNull();
   return;
}




#else



/*
***************************************************************************
*                                                                         *
* When MPI is not available, this communication class is practically empty. *
* public methods are trivial and private methods do not exist.            *
*                                                                         *
***************************************************************************
*/



AsyncCommGroup::AsyncCommGroup( const int nchild )
   : d_nchild(nchild)
{
   return;
}


AsyncCommGroup::~AsyncCommGroup()
{
   return;
}





bool AsyncCommGroup::checkOperation()
{
   return true;
}


void AsyncCommGroup::waitOperation()
{
   return;
}





bool AsyncCommGroup::beginBcast( int *buffer, int size )
{
   return true;
}


bool AsyncCommGroup::checkBcast()
{
   return true;
}






bool AsyncCommGroup::beginGather( int *buffer, int size )
{
   return true;
}




bool AsyncCommGroup::checkGather()
{
   return true;
}





bool AsyncCommGroup::beginSumReduce( int *buffer, int size )
{
   return true;
}


bool AsyncCommGroup::checkSumReduce()
{
   return true;
}


void AsyncCommGroup::setGroupAndRootRank( const Array<int> &group,
                                          const int root_rank )
{
   return;
}


void AsyncCommGroup::setGroupAndRootIndex( const Array<int> &group,
                                           const int root_index )
{
   return;
}





SAMRAI_MPI::request *AsyncCommGroup::getRequestPointer() const
{
   return NULL;
}





void AsyncCommGroup::checkMPIParams()
{
   return;
}



#endif


}
}
