//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/standard/CoarsenCopyTransaction.C $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Communication transaction for data copies during data coarsening
//
 
#ifndef included_xfer_CoarsenCopyTransaction_C
#define included_xfer_CoarsenCopyTransaction_C

#include "CoarsenCopyTransaction.h"

#include "Patch.h"
#include "PatchData.h"
#include "tbox/SAMRAI_MPI.h"
#include "CoarsenClasses.h"

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
   TBOX_ASSERT(coarsen_items != (const typename CoarsenClasses<DIM>::Data**)NULL);
   TBOX_ASSERT(num_coarsen_items >= 0);
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
   TBOX_ASSERT(!dst_level.isNull());
   TBOX_ASSERT(!src_level.isNull());
   TBOX_ASSERT(!overlap.isNull());
   TBOX_ASSERT(dst_patch >= 0 && dst_patch < dst_level->getNumberOfPatches());
   TBOX_ASSERT(src_patch >= 0 && src_patch < src_level->getNumberOfPatches());
   TBOX_ASSERT(coarsen_item_id >= 0);
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
   if (getSourceProcessor() == tbox::SAMRAI_MPI::getRank()) {
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
void CoarsenCopyTransaction<DIM>::printClassData(std::ostream& stream) const
{
   stream << "Coarsen Copy Transaction"                            << std::endl;
   stream << "   coarsen item array:        " 
          << (typename CoarsenClasses<DIM>::Data**)s_coarsen_items << std::endl;
   stream << "   num coarsen items:      " << s_num_coarsen_items << std::endl;
   stream << "   destination patch:      " << d_dst_patch       << std::endl;
   stream << "   source patch:           " << d_src_patch       << std::endl;
   stream << "   coarsen item id:        " << d_coarsen_item_id << std::endl;
   stream << "   destination patch data: " 
          << s_coarsen_items[d_coarsen_item_id]->d_dst          << std::endl;
   stream << "   source patch data:      " 
          << s_coarsen_items[d_coarsen_item_id]->d_src          << std::endl;
   stream << "   incoming bytes:         " << d_incoming_bytes  << std::endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes  << std::endl;
   stream << "   destination level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                     << std::endl;
   stream << "   source level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level                     << std::endl;
   stream << "   overlap:                "                      << std::endl;
   d_overlap->print(stream);
}

}
}
#endif
