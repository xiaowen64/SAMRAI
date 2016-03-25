//
// File:	StatTransaction.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Communication transaction structure for statistic data copies
//
 
#include "tbox/StatTransaction.h"

namespace SAMRAI {
   namespace tbox {

StatTransaction::StatTransaction(
   Pointer<Statistic> stat,
   int src_proc_id,
   int dst_proc_id)
{
   d_stat = stat;
   d_src_id = src_proc_id;
   d_dst_id = dst_proc_id;
}

StatTransaction::~StatTransaction()
{
}
 
bool StatTransaction::canEstimateIncomingMessageSize()
{
   return(d_stat->canEstimateDataStreamSize());
}

int StatTransaction::computeIncomingMessageSize()
{
   return(d_stat->getDataStreamSize());
}

int StatTransaction::computeOutgoingMessageSize()
{
   return(d_stat->getDataStreamSize());
}

int StatTransaction::getSourceProcessor()
{
   return(d_src_id);
}

int StatTransaction::getDestinationProcessor()
{
   return(d_dst_id);
}

void StatTransaction::packStream(AbstractStream& stream)
{
   d_stat->packStream(stream);
}

void StatTransaction::unpackStream(AbstractStream& stream)
{
   d_stat->unpackStream(stream);
}

void StatTransaction::copyLocalData()
{
   // Nothing to do here!
}

void StatTransaction::printClassData(ostream& stream) const
{
   stream << "Stat Transaction" << endl;
   stream << "   source processor:   " << d_src_id      << endl;
   stream << "   destination processor:   " << d_dst_id      << endl;
   stream << "   stat name:   " << d_stat->getName() << endl;
}

}
}
