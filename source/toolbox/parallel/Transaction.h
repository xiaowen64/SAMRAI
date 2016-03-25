//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/parallel/Transaction.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Abstract base class for all schedule transactions
//

#ifndef included_tbox_Transaction
#define included_tbox_Transaction

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "tbox/AbstractStream.h"

namespace SAMRAI {
   namespace tbox {


/**
 * Class Transaction describes a single communication between two 
 * processors or a local data copy.  It is an abstract base class for each 
 * data transaction in a communication schedule.
 */

class Transaction : public DescribedClass
{
public:
   /**
    * The constructor for transaction does nothing interesting.
    */
   Transaction();

   /**
    * The virtual destructor for transaction does nothing interesting.
    */
   virtual ~Transaction();

   /**
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communications protocol kicks in and the message size is transmitted
    * between nodes.
    */
   virtual bool canEstimateIncomingMessageSize() = 0;

   /**
    * Return the amount of buffer space needed for the incoming message.
    * This routine is only called if the transaction can estimate the
    * size of the incoming message.
    */
   virtual int computeIncomingMessageSize() = 0;

   /**
    * Return the buffer space needed for the outgoing message.
    */
   virtual int computeOutgoingMessageSize() = 0;

   /**
    * Return the sending processor for the communications transaction.
    */
   virtual int getSourceProcessor() = 0;

   /**
    * Return the receiving processor for the communications transaction.
    */
   virtual int getDestinationProcessor() = 0;

   /**
    * Pack the transaction data into the message stream.
    */
   virtual void packStream(AbstractStream& stream) = 0;

   /**
    * Unpack the transaction data from the message stream.
    */
   virtual void unpackStream(AbstractStream& stream) = 0;

   /**
    * Perform the local data copy for the transaction.
    */
   virtual void copyLocalData() = 0;

   /**
    * Print out transaction information.
    */
   virtual void printClassData(std::ostream& stream) const = 0;

private:
   Transaction(const Transaction&);	// not implemented
   void operator=(const Transaction&);	// not implemented

};


}
}

#endif
