//
// File:	StatTransaction.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Communication transaction structure for statistic data copies
//
 
#ifndef included_tbox_StatTransaction
#define included_tbox_StatTransaction

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Statistic
#include "tbox/Statistic.h"
#endif
#ifndef included_tbox_Transaction
#include "tbox/Transaction.h"
#endif

namespace SAMRAI {
   namespace tbox {


/**
 * A stattistic transaction represents a simple copy communication 
 * transaction between two processors for sending and gathering statistic 
 * information generated on different processors.
 *
 * @see tbox::Schedule
 * @see tbox::Transaction
 */

class StatTransaction : public Transaction
{
public:
   /**
    * Create a transaction for communicating local statistic information.
    * This transaction will be responsible for either: (1) packing a
    * message stream with statistic information if this processor number
    * is the same as the given source processor id, or (2) unpacking 
    * tatistic information from the message stream if this processor number
    * is the same as the given destination processor.  The statistic pointer
    * passed in through the argument list must be non-null and will be used
    * as either the source or destination statistic depending on whether
    * case (1) or (2) applies.  
    *
    * Note that generally this transaction class is used to pass information
    * between two different processors and unexpected behavior may result
    * if the source and destination processors are the same.  Also, note
    * that the copyLocalData() routine has an empty implementation. 
    */
   StatTransaction(Pointer<Statistic> stat,
                        int src_proc_id,
                        int dst_proc_id);

   /**
    * The virtual destructor for the copy transaction releases all
    * memory associated with the transaction.
    */
   virtual ~StatTransaction();

   /**
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communications protocol kicks in and the message size is transmitted
    * between nodes.
    */
   virtual bool canEstimateIncomingMessageSize();

   /**
    * Return the amount of buffer space needed for the incoming message.
    * This routine is only called if the transaction can estimate the
    * size of the incoming message.
    */
   virtual int computeIncomingMessageSize();

   /**
    * Return the buffer space needed for the outgoing message.
    */
   virtual int computeOutgoingMessageSize();

   /**
    * Return the sending processor for the communications transaction.
    */
   virtual int getSourceProcessor();

   /**
    * Return the receiving processor for the communications transaction.
    */
   virtual int getDestinationProcessor();

   /**
    * Pack the transaction data into the message stream.
    */
   virtual void packStream(AbstractStream& stream);

   /**
    * Unpack the transaction data from the message stream.
    */
   virtual void unpackStream(AbstractStream& stream);

   /**
    * Perform the local data copy for the transaction.  This function
    * drops through as it is not needed.
    */
   virtual void copyLocalData();

   /**
    * Print transaction information to given output stream.
    */
   virtual void printClassData(ostream& stream) const;

private:
   StatTransaction(const StatTransaction&);	// not implemented
   void operator=(const StatTransaction&);	// not implemented
   
   Pointer<Statistic> d_stat;
   int d_src_id;
   int d_dst_id;

};

}
}
#endif
