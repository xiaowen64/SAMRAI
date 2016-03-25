//
// File:	RefineCopyTransaction.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 684 $
// Modified:	$Date: 2005-10-21 14:59:38 -0700 (Fri, 21 Oct 2005) $
// Description:	Communication transaction for data copies during data refining
//
 
#ifndef included_xfer_RefineCopyTransaction_C
#define included_xfer_RefineCopyTransaction_C

#include "RefineCopyTransaction.h"

#include "Patch.h"
#include "PatchData.h"
#include "tbox/MPI.h"
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
* Initialization, set/unset functions for static array of refine items. *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
const typename RefineClasses<DIM>::Data** 
   RefineCopyTransaction<DIM>::s_refine_items = 
      (const typename RefineClasses<DIM>::Data**)NULL;
template<int DIM> int RefineCopyTransaction<DIM>::s_num_refine_items = 0;

template<int DIM> void RefineCopyTransaction<DIM>::setRefineItems(
   const typename RefineClasses<DIM>::Data** refine_items,
   int num_refine_items)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(refine_items != (const typename RefineClasses<DIM>::Data**)NULL);
   assert(num_refine_items >= 0);
#endif
   s_refine_items = refine_items;
   s_num_refine_items = num_refine_items;
}

template<int DIM> void RefineCopyTransaction<DIM>::unsetRefineItems()
{
   s_refine_items = (const typename RefineClasses<DIM>::Data**)NULL;
   s_num_refine_items = 0;
}

/*
*************************************************************************
*                                                                       *
* Constructor sets state of transaction.                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>  RefineCopyTransaction<DIM>::RefineCopyTransaction(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
   int dst_patch,
   int src_patch,
   int refine_item_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!dst_level.isNull());
   assert(!src_level.isNull());
   assert(!overlap.isNull());
   assert(dst_patch >= 0 && dst_patch < dst_level->getNumberOfPatches());
   assert(src_patch >= 0 && src_patch < src_level->getNumberOfPatches());
   assert(refine_item_id >= 0);
   // Note: s_num_refine_items cannot be used at this point!
#endif
   d_dst_level        = dst_level;
   d_src_level        = src_level;
   d_overlap          = overlap;
   d_dst_patch        = dst_patch;
   d_src_patch        = src_patch;
   d_refine_item_id   = refine_item_id;
   d_incoming_bytes   = 0;
   d_outgoing_bytes   = 0;
}

template<int DIM>  RefineCopyTransaction<DIM>::~RefineCopyTransaction()
{
}

/*
*************************************************************************
*                                                                       *
* Functions overridden in tbox::Transaction base class.                 *
*                                                                       *
*************************************************************************
*/
 
template<int DIM> bool RefineCopyTransaction<DIM>::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (getSourceProcessor() == tbox::MPI::getRank()) {
      can_estimate = 
         d_src_level->getPatch(d_src_patch)
                    ->getPatchData(s_refine_items[d_refine_item_id]->
                                   d_src)
                    ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate = 
         d_dst_level->getPatch(d_dst_patch)
                    ->getPatchData(s_refine_items[d_refine_item_id]->
                                   d_scratch)
                    ->canEstimateStreamSizeFromBox();
   }
   return(can_estimate);
}

template<int DIM> int RefineCopyTransaction<DIM>::computeIncomingMessageSize()
{
   d_incoming_bytes = 
      d_dst_level->getPatch(d_dst_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_scratch)
                 ->getDataStreamSize(*d_overlap);
   return(d_incoming_bytes);
}

template<int DIM> int RefineCopyTransaction<DIM>::computeOutgoingMessageSize()
{
   d_outgoing_bytes = 
      d_src_level->getPatch(d_src_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_src)
                 ->getDataStreamSize(*d_overlap);
   return(d_outgoing_bytes);
}

template<int DIM> int RefineCopyTransaction<DIM>::getSourceProcessor()
{
   return(d_src_level->getMappingForPatch(d_src_patch));
}

template<int DIM> int RefineCopyTransaction<DIM>::getDestinationProcessor()
{
   return(d_dst_level->getMappingForPatch(d_dst_patch));
}

template<int DIM> void RefineCopyTransaction<DIM>::packStream(tbox::AbstractStream& stream)
{
   d_src_level->getPatch(d_src_patch)
              ->getPatchData(s_refine_items[d_refine_item_id]->
                             d_src)
              ->packStream(stream, *d_overlap);
}

template<int DIM> void RefineCopyTransaction<DIM>::unpackStream(tbox::AbstractStream& stream)
{
   d_dst_level->getPatch(d_dst_patch)
              ->getPatchData(s_refine_items[d_refine_item_id]->
                             d_scratch)
              ->unpackStream(stream, *d_overlap);
}

template<int DIM> void RefineCopyTransaction<DIM>::copyLocalData()
{
   d_dst_level->getPatch(d_dst_patch)
              ->getPatchData(s_refine_items[d_refine_item_id]->
                             d_scratch)
              ->copy(*d_src_level->getPatch(d_src_patch)
                                 ->getPatchData(
                                     s_refine_items[d_refine_item_id]->
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
void RefineCopyTransaction<DIM>::printClassData(ostream& stream) const
{
   stream << "Refine Copy Transaction"                            << endl;
   stream << "   refine item array:        "
          << (typename RefineClasses<DIM>::Data**)s_refine_items  << endl;
   stream << "   num refine items:       " << s_num_refine_items << endl;
   stream << "   destination patch:      " << d_dst_patch      << endl;
   stream << "   source patch:           " << d_src_patch      << endl;
   stream << "   refine item id:         " << d_refine_item_id << endl;
   stream << "   destination patch data: " 
          << s_refine_items[d_refine_item_id]->d_scratch       << endl;
   stream << "   source patch data:      " 
          << s_refine_items[d_refine_item_id]->d_src           << endl;
   stream << "   incoming bytes:         " << d_incoming_bytes << endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes << endl;
   stream << "   destination level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                    << endl;
   stream << "   source level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                    << endl;
   stream << "   overlap:                "                     << endl;
   d_overlap->print(stream);
}

}
}
#endif
