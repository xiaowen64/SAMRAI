//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/RefineCopyTransaction.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Communication transaction for data copies during data refining
//
 
#ifndef included_xfer_RefineCopyTransaction_C
#define included_xfer_RefineCopyTransaction_C

#include "RefineCopyTransaction.h"

#include "Patch.h"
#include "PatchData.h"
#include "tbox/SAMRAI_MPI.h"

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
   TBOX_ASSERT(refine_items != (const typename RefineClasses<DIM>::Data**)NULL);
   TBOX_ASSERT(num_refine_items >= 0);
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
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(!overlap.isNull());
   TBOX_ASSERT(dst_patch >= 0 && dst_patch < dst_level->getNumberOfPatches());
   TBOX_ASSERT(src_patch >= 0 && src_patch < src_level->getNumberOfPatches());
   TBOX_ASSERT(refine_item_id >= 0);
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
   if (getSourceProcessor() == tbox::SAMRAI_MPI::getRank()) {
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
void RefineCopyTransaction<DIM>::printClassData(std::ostream& stream) const
{
   stream << "Refine Copy Transaction"                            << std::endl;
   stream << "   refine item array:        "
          << (typename RefineClasses<DIM>::Data**)s_refine_items  << std::endl;
   stream << "   num refine items:       " << s_num_refine_items << std::endl;
   stream << "   destination patch:      " << d_dst_patch      << std::endl;
   stream << "   source patch:           " << d_src_patch      << std::endl;
   stream << "   refine item id:         " << d_refine_item_id << std::endl;
   stream << "   destination patch data: " 
          << s_refine_items[d_refine_item_id]->d_scratch       << std::endl;
   stream << "   source patch data:      " 
          << s_refine_items[d_refine_item_id]->d_src           << std::endl;
   stream << "   incoming bytes:         " << d_incoming_bytes << std::endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes << std::endl;
   stream << "   destination level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                    << std::endl;
   stream << "   source level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                    << std::endl;
   stream << "   overlap:                "                     << std::endl;
   d_overlap->print(stream);
}

}
}
#endif
