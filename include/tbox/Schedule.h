//
// File:	Schedule.h
// Package:	SAMRAI communication and data transfer package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Schedule of communication transactions between processors
//

#ifndef included_tbox_Schedule
#define included_tbox_Schedule

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_MPI
#include "tbox/MPI.h"
#endif
#ifndef included_tbox_MessageStream
#include "tbox/MessageStream.h"
#endif
#ifndef included_tbox_Transaction
#include "tbox/Transaction.h"
#endif

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
    * transaction objects.  Member function {\tt finalizeCommunication()}
    * must be called to finish the message communication.
    */
   void beginCommunication();

   /*!
    * Finish the communication and deliver the messages.  This member
    * function completes communication began by {\tt beginCommunication()}.
    */
   void finalizeCommunication();

   /*!
    * Print class data to the specified output stream.
    */
   void printClassData(ostream& stream) const;

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
