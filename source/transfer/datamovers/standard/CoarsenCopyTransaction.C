//
// File:	CoarsenCopyTransaction.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 684 $
// Modified:	$Date: 2005-10-21 14:59:38 -0700 (Fri, 21 Oct 2005) $
// Description:	Communication transaction for data copies during data coarsening
//
 
#ifndef included_xfer_CoarsenCopyTransaction_C
#define included_xfer_CoarsenCopyTransaction_C

#include "CoarsenCopyTransaction.h"

#include "Patch.h"
#include "PatchData.h"
#include "tbox/MPI.h"
#include "CoarsenClasses.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
    namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

/*
*************************************************************************
*                                                                       *
* Initialization, set/unset functions for static array of coarsen items.*
*                                                                       *
*************************************************************************
*/

template<int DIM> const typename CoarsenClasses<DIM>::Data** 
   CoarsenCopyTransaction<DIM>::s_coarsen_items = 
       (const typename CoarsenClasses<DIM>::Data**)NULL;
template<int DIM> int CoarsenCopyTransaction<DIM>::s_num_coarsen_items = 0;

template<int DIM> void CoarsenCopyTransaction<DIM>::setCoarsenItems(
   const typename CoarsenClasses<DIM>::Data** coarsen_items,
   int num_coarsen_items)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(coarsen_items != (const typename CoarsenClasses<DIM>::Data**)NULL);
   assert(num_coarsen_items >= 0);
#endif
   s_coarsen_items = coarsen_items;
   s_num_coarsen_items = num_coarsen_items;
}

template<int DIM> void CoarsenCopyTransaction<DIM>::unsetCoarsenItems()
{
   s_coarsen_items = (const typename CoarsenClasses<DIM>::Data**)NULL;
   s_num_coarsen_items = 0;
}

/*
*************************************************************************
*                                                                       *
* Constructor sets state of transaction.                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>  CoarsenCopyTransaction<DIM>::CoarsenCopyTransaction(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
   int dst_patch,
   int src_patch,
   int coarsen_item_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst_level.isNull());
   assert(!src_level.isNull());
   assert(!overlap.isNull());
   assert(dst_patch >= 0 && dst_patch < dst_level->getNumberOfPatches());
   assert(src_patch >= 0 && src_patch < src_level->getNumberOfPatches());
   assert(coarsen_item_id >= 0);
   // Note: s_num_coarsen_items cannot be used at this point!
#endif

   d_dst_level        = dst_level;
   d_src_level        = src_level;
   d_overlap          = overlap;
   d_dst_patch        = dst_patch;
   d_src_patch        = src_patch;
   d_coarsen_item_id  = coarsen_item_id;
   d_incoming_bytes   = 0;
   d_outgoing_bytes   = 0;
}

template<int DIM>  CoarsenCopyTransaction<DIM>::~CoarsenCopyTransaction()
{
}

/*
*************************************************************************
*                                                                       *
* Functions overridden in tbox::Transaction base class.                  *
*                                                                       *
*************************************************************************
*/
 
template<int DIM> bool CoarsenCopyTransaction<DIM>::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (getSourceProcessor() == tbox::MPI::getRank()) {
      can_estimate = 
         d_src_level->getPatch(d_src_patch)
                    ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                                   d_src)
                    ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate = 
         d_dst_level->getPatch(d_dst_patch)
                    ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                                   d_dst)
                    ->canEstimateStreamSizeFromBox();
   }
   return(can_estimate);
}

template<int DIM> int CoarsenCopyTransaction<DIM>::computeIncomingMessageSize()
{
   d_incoming_bytes = 
      d_dst_level->getPatch(d_dst_patch)
                 ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                                d_dst)
                 ->getDataStreamSize(*d_overlap);
   return(d_incoming_bytes);
}

template<int DIM> int CoarsenCopyTransaction<DIM>::computeOutgoingMessageSize()
{
   d_outgoing_bytes = 
      d_src_level->getPatch(d_src_patch)
                 ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                                d_src)
                 ->getDataStreamSize(*d_overlap);
   return(d_outgoing_bytes);
}

template<int DIM> int CoarsenCopyTransaction<DIM>::getSourceProcessor()
{
   return(d_src_level->getMappingForPatch(d_src_patch));
}

template<int DIM> int CoarsenCopyTransaction<DIM>::getDestinationProcessor()
{
   return(d_dst_level->getMappingForPatch(d_dst_patch));
}

template<int DIM> void CoarsenCopyTransaction<DIM>::packStream(tbox::AbstractStream& stream)
{
   d_src_level->getPatch(d_src_patch)
              ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                             d_src)
              ->packStream(stream, *d_overlap);
}

template<int DIM> void CoarsenCopyTransaction<DIM>::unpackStream(tbox::AbstractStream& stream)
{
   d_dst_level->getPatch(d_dst_patch)
              ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                             d_dst)
              ->unpackStream(stream, *d_overlap);
}

template<int DIM> void CoarsenCopyTransaction<DIM>::copyLocalData()
{
   d_dst_level->getPatch(d_dst_patch)
              ->getPatchData(s_coarsen_items[d_coarsen_item_id]->
                             d_dst)
              ->copy(*d_src_level->getPatch(d_src_patch)
                                 ->getPatchData(
                                     s_coarsen_items[d_coarsen_item_id]->
                                        d_src), *d_overlap);
}

/*
*************************************************************************
*                                                                       *
* Function to print state of transaction.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
void CoarsenCopyTransaction<DIM>::printClassData(ostream& stream) const
{
   stream << "Coarsen Copy Transaction"                            << endl;
   stream << "   coarsen item array:        " 
          << (typename CoarsenClasses<DIM>::Data**)s_coarsen_items << endl;
   stream << "   num coarsen items:      " << s_num_coarsen_items << endl;
   stream << "   destination patch:      " << d_dst_patch       << endl;
   stream << "   source patch:           " << d_src_patch       << endl;
   stream << "   coarsen item id:        " << d_coarsen_item_id << endl;
   stream << "   destination patch data: " 
          << s_coarsen_items[d_coarsen_item_id]->d_dst          << endl;
   stream << "   source patch data:      " 
          << s_coarsen_items[d_coarsen_item_id]->d_src          << endl;
   stream << "   incoming bytes:         " << d_incoming_bytes  << endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes  << endl;
   stream << "   destination level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                     << endl;
   stream << "   source level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                     << endl;
   stream << "   overlap:                "                      << endl;
   d_overlap->print(stream);
}

}
}
#endif
