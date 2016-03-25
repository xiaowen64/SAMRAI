//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/parallel/Schedule.h $
// Package:	SAMRAI communication and data transfer package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2410 $
// Modified:	$LastChangedDate: 2008-10-08 14:30:40 -0700 (Wed, 08 Oct 2008) $
// Description:	Schedule of communication transactions between processors
//

#ifndef included_tbox_Schedule
#define included_tbox_Schedule

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MessageStream.h"
#include "tbox/Transaction.h"

namespace SAMRAI {
   namespace tbox {


/*!
 * @brief Class Schedule is used to construct and execute a set of 
 * data communication transactions.  Each transaction represents some 
 * data dependency and exchange between two processors, or locally involving
 * a single processor.  Once a communication schedule is constructed, transactions 
 * are provided to the schedule, using either the addTransaction() method or 
 * the appendTransaction() method.  The schedule is then executed forcing
 * the communication, either interprocessor or local to occur.  The basic idea 
 * behind the schedule is that it enables the cost of assembling communication
 * dependencies and data transfers over many communication phases.
 * 
 * Note that since the transactions are stored in lists, the "add" and "append" 
 * mimick the semantics of the List class.  That is, addTransaction() will 
 * put the transaction at the head of the list, while appendTransaction() will 
 * put the transaction at the end of the list.  This flexibility is provided for
 * situations where the order of transaction execution matters.  Regardless of 
 * which method is used to assemble the transactions, they will be executed in the
 * order in which they appear in the list.
 *
 * @see tbox::Transaction
 * @see tbox::List
 */

class Schedule : public DescribedClass
{
public:
   /*!
    * Create an empty schedule with no transactions.
    */
   Schedule();

   /*!
    * The destructor deletes the schedule and all associated storage.  Note
    * that the schedule should not be deleted during a communication phase.
    */
   virtual ~Schedule();

   /*!
    * Add a data transaction to the head of the list of transactions already
    * assembled in the schedule.  The transaction must involve the local 
    * processor as either a source or destination or both.  If the transaction 
    * does not include the local processor, then the transaction
    * is not placed on the schedule.
    *
    * @param transaction  Pointer to transaction added to the schedule.
    */
   void addTransaction(const Pointer<Transaction>& transaction);

   /*!
    * Append a data transaction to the tail of the list of transactions already
    * assembled in the schedule.  The transaction must involve the local
    * processor as either a source or destination or both.  If the transaction
    * does not include the local processor, then the transaction
    * is not placed on the schedule.
    *
    * @param transaction  Pointer to transaction appended to the schedule.
    */
   void appendTransaction(const Pointer<Transaction>& transaction);

   /*!
    * Perform the communication described by the schedule.
    */
   void communicate();

   /*!
    * Begin the communication process but do not deliver data to the
    * transaction objects.  Member function <TT>finalizeCommunication()</TT>
    * must be called to finish the message communication.
    */
   void beginCommunication();

   /*!
    * Finish the communication and deliver the messages.  This member
    * function completes communication began by <TT>beginCommunication()</TT>.
    */
   void finalizeCommunication();

   /*!
    * Print class data to the specified output stream.
    */
   void printClassData(std::ostream& stream) const;

   /*!
    * This internal meassage stream structure must be declared public for the 
    * Sun CC compiler.
    */
   struct ScheduleMessageStream {
      int d_bytes_in_stream;
      bool d_must_communicate_byte_size;
      bool d_stream_in_use;
      Pointer<MessageStream> d_stream;
#ifdef HAVE_MPI
      MPI_Request d_request_id;
#endif
   };

private:


   /*!
    * @brief Set up things for the entire class.
    */
   void firstConstructorTasks();

   /*!
     Free static timers.

     To be called by shutdown registry to make sure
     memory for timers does not leak.
   */
   static void freeTimers();

   void calculateSendSizes();
   void calculateReceiveSizes();
   void postMessageReceives();
   void sendMessages();
   void performLocalCopies();
   void processIncomingMessages();
   void deallocateSendBuffers();

   Schedule(const Schedule&);		// not implemented
   void operator=(const Schedule&);	// not implemented

   int d_nnodes;

#ifdef LACKS_NAMESPACE_IN_DECLARE
   Array< ScheduleMessageStream > d_incoming;
   Array< ScheduleMessageStream > d_outgoing;
#else
   Array< Schedule::ScheduleMessageStream > d_incoming;
   Array< Schedule::ScheduleMessageStream > d_outgoing;
#endif

   Array< List< Pointer< Transaction > > > d_send_set;
   Array< List< Pointer< Transaction > > > d_recv_set;


   List< Pointer< Transaction > > d_local_set;

};


}
}

#endif
