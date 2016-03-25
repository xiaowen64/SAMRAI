//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/timers/StatTransaction.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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

void StatTransaction::printClassData(std::ostream& stream) const
{
   stream << "Stat Transaction" << std::endl;
   stream << "   source processor:   " << d_src_id      << std::endl;
   stream << "   destination processor:   " << d_dst_id      << std::endl;
   stream << "   stat name:   " << d_stat->getName() << std::endl;
}

}
}
