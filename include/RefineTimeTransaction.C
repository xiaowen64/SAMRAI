//
// File:	RefineTimeTransaction.C
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 684 $
// Modified:	$Date: 2005-10-21 14:59:38 -0700 (Fri, 21 Oct 2005) $
// Description:	Communication transaction for time interpolation during data refining
//

#ifndef included_xfer_RefineTimeTransaction_C
#define included_xfer_RefineTimeTransaction_C

#include <typeinfo>
using namespace std;

#include "RefineTimeTransaction.h"

#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "tbox/MPI.h"
#include "tbox/Utilities.h"
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
* Initialization, set/unset functions for static array of refine items  *
* and interpolation time.                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> double RefineTimeTransaction<DIM>::s_time = 0.0;

template<int DIM> 
const typename RefineClasses<DIM>::Data** 
   RefineTimeTransaction<DIM>::s_refine_items = 
      (const typename RefineClasses<DIM>::Data**)NULL;
template<int DIM> int RefineTimeTransaction<DIM>::s_num_refine_items = 0;

template<int DIM> 
void RefineTimeTransaction<DIM>::setTransactionTime(const double time)
{
   s_time = time;
}

template<int DIM> void RefineTimeTransaction<DIM>::setRefineItems(
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

template<int DIM> void RefineTimeTransaction<DIM>::unsetRefineItems()
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

template<int DIM>  RefineTimeTransaction<DIM>::RefineTimeTransaction(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
   int dst_patch,
   int src_patch,
   const hier::Box<DIM>& box,
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
   d_box              = box;
   d_refine_item_id   = refine_item_id;
}

template<int DIM>  RefineTimeTransaction<DIM>::~RefineTimeTransaction()
{
}

/*
*************************************************************************
*                                                                       *
* Functions overridden in tbox::Transaction base class.                 *
*                                                                       *
*************************************************************************
*/
 
template<int DIM> bool RefineTimeTransaction<DIM>::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (getSourceProcessor() == tbox::MPI::getRank()) {
      can_estimate = 
         d_src_level->getPatch(d_src_patch)
                    ->getPatchData(s_refine_items[d_refine_item_id]->
                                   d_src_told)
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

template<int DIM> int RefineTimeTransaction<DIM>::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_level->getPatch(d_dst_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_scratch)
                 ->getDataStreamSize(*d_overlap);
   return(d_incoming_bytes);
}

template<int DIM> int RefineTimeTransaction<DIM>::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_level->getPatch(d_src_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_src_told)
                 ->getDataStreamSize(*d_overlap);
   return(d_outgoing_bytes);
}

template<int DIM> int RefineTimeTransaction<DIM>::getSourceProcessor()
{
   return(d_src_level->getMappingForPatch(d_src_patch));
}

template<int DIM> int RefineTimeTransaction<DIM>::getDestinationProcessor()
{
   return(d_dst_level->getMappingForPatch(d_dst_patch));
}

template<int DIM> void RefineTimeTransaction<DIM>::packStream(tbox::AbstractStream& stream)
{
   const tbox::Pointer< hier::Patch<DIM> >& patch = d_src_level->getPatch(d_src_patch);

   tbox::Pointer< hier::PatchData<DIM> > temp =
      d_src_level->getPatchDescriptor()
                 ->getPatchDataFactory(s_refine_items[d_refine_item_id]->
                                       d_src_told)
                 ->allocate(d_box);
   temp->setTime(s_time);

   timeInterpolate(
      temp,
      patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_told),
      patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_tnew) );

   temp->packStream(stream, *d_overlap);
}

template<int DIM> void RefineTimeTransaction<DIM>::unpackStream(tbox::AbstractStream& stream)
{
   d_dst_level->getPatch(d_dst_patch)
              ->getPatchData(s_refine_items[d_refine_item_id]->
                             d_scratch)
              ->unpackStream(stream, *d_overlap);
}

template<int DIM> void RefineTimeTransaction<DIM>::copyLocalData()
{
   const tbox::Pointer< hier::Patch<DIM> >& patch = d_src_level->getPatch(d_src_patch);

   /*
    * If there is no offset between the source and destination, then 
    * time interpolate directly to the destination patchdata.  Otherwise,
    * time interpolate into a temporary patchdata and copy the result
    * to the destination patchdata.
    */
   if (d_overlap->getSourceOffset() == hier::IntVector<DIM>(0)) {

      timeInterpolate(
         d_dst_level->getPatch(d_dst_patch)->
                getPatchData(s_refine_items[d_refine_item_id]->
                             d_scratch),
         patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_told),
         patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_tnew) );

   } else {

      tbox::Pointer< hier::PatchData<DIM> > temp =
         d_src_level->getPatchDescriptor()
                    ->getPatchDataFactory(s_refine_items[d_refine_item_id]->
                                          d_src_told)
                    ->allocate(d_box);
      temp->setTime(s_time);
   
      timeInterpolate(
         temp,
         patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_told),
         patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_tnew) );
   
      d_dst_level->getPatch(d_dst_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->d_scratch) 
                 ->copy(*temp,*d_overlap);

   }

}

template<int DIM> void RefineTimeTransaction<DIM>::timeInterpolate(
   const tbox::Pointer< hier::PatchData<DIM> >& pd_dst,
   const tbox::Pointer< hier::PatchData<DIM> >& pd_old,
   const tbox::Pointer< hier::PatchData<DIM> >& pd_new)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!pd_old.isNull());
   assert(!pd_dst.isNull());
   assert(tbox::Utilities::deq(pd_dst->getTime(), s_time));
#endif
   if (tbox::Utilities::deq(pd_old->getTime(), s_time)) {
      s_refine_items[d_refine_item_id]->
         d_optime->timeInterpolate(*pd_dst, d_box, *pd_old, *pd_old);
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!pd_new.isNull());
   assert(pd_old->getTime() < s_time);
   assert(pd_new->getTime() >= s_time);
#endif
      s_refine_items[d_refine_item_id]->
         d_optime->timeInterpolate(*pd_dst, d_box, *pd_old, *pd_new);
   }
}

/*
*************************************************************************
*                                                                       *
* Function to print state of transaction.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
void RefineTimeTransaction<DIM>::printClassData(ostream& stream) const
{
   stream << "Refine Time Transaction" << endl;
   stream << "   transaction time:        " << s_time            << endl;
   stream << "   refine item array:        "
          << (typename RefineClasses<DIM>::Data**)s_refine_items << endl;
   stream << "   num refine items:        " << s_num_refine_items << endl;
   stream << "   destination patch:       " << d_dst_patch       << endl;
   stream << "   source patch:            " << d_src_patch       << endl;
   stream << "   time interpolation box:  " << d_box             << endl;
   stream << "   refine item id:          " << d_refine_item_id  << endl;
   stream << "   destination patch data:  " 
          << s_refine_items[d_refine_item_id]->d_scratch         << endl;
   stream << "   source (old) patch data: " 
          << s_refine_items[d_refine_item_id]->d_src_told        << endl;
   stream << "   source (new) patch data: " 
          << s_refine_items[d_refine_item_id]->d_src_tnew        << endl;
   stream << "   time interpolation name: " 
          << typeid(*s_refine_items[d_refine_item_id]->d_optime).name() << endl;
   stream << "   incoming bytes:          " << d_incoming_bytes  << endl;
   stream << "   outgoing bytes:          " << d_outgoing_bytes  << endl;
   stream << "   destination level:           "
          << (hier::PatchLevel<DIM>*)d_src_level                    << endl;
   stream << "   source level:           "
          << (hier::PatchLevel<DIM>*)d_src_level                    << endl;
   stream << "   overlap:                 " << endl;
   d_overlap->print(stream);
}

}
}
#endif
