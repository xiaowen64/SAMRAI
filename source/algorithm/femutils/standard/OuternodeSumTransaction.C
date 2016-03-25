//
// File:        OuternodeSumTransaction.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 684 $
// Modified:    $Date: 2005-10-21 14:59:38 -0700 (Fri, 21 Oct 2005) $
// Description: Communication transaction for summing outernode data
//

#ifndef included_algs_OuternodeSumTransaction_C
#define included_algs_OuternodeSumTransaction_C
 
#include "OuternodeSumTransaction.h"

#include "Patch.h"
#include "PatchData.h"
#include "ArrayDataBasicOps.h"
#include "NodeGeometry.h"
#include "OuternodeData.h"
#include "tbox/MPI.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
    namespace algs {

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

template<int DIM> const typename xfer::RefineClasses<DIM>::Data**
   OuternodeSumTransaction<DIM>::s_refine_items =
      (const typename xfer::RefineClasses<DIM>::Data**)NULL;
template<int DIM> int OuternodeSumTransaction<DIM>::s_num_refine_items = 0;

template<int DIM> void OuternodeSumTransaction<DIM>::setRefineItems(
   const typename xfer::RefineClasses<DIM>::Data** refine_items,
   int num_refine_items)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(refine_items != (const typename xfer::RefineClasses<DIM>::Data**)NULL);
   assert(num_refine_items >= 0);
#endif
   s_refine_items = refine_items;
   s_num_refine_items = num_refine_items;
}

template<int DIM> void OuternodeSumTransaction<DIM>::unsetRefineItems()
{
   s_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;
   s_num_refine_items = 0;
}

/*
*************************************************************************
*                                                                       *
* Constructor sets state of transaction.                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> OuternodeSumTransaction<DIM>::OuternodeSumTransaction(
   tbox::Pointer<hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer<hier::PatchLevel<DIM> > src_level,
   tbox::Pointer<hier::BoxOverlap<DIM> > overlap,
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

template<int DIM> OuternodeSumTransaction<DIM>::~OuternodeSumTransaction()
{
}

/*
*************************************************************************
*                                                                       *
* Functions overridden in tbox::Transaction base class.                 *
*                                                                       *
*************************************************************************
*/
 
template<int DIM> bool 
OuternodeSumTransaction<DIM>::canEstimateIncomingMessageSize()
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

template<int DIM> int 
OuternodeSumTransaction<DIM>::computeIncomingMessageSize()
{
   d_incoming_bytes = 
      d_dst_level->getPatch(d_dst_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_scratch)
                 ->getDataStreamSize(*d_overlap);
   return(d_incoming_bytes);
}

template<int DIM> int 
OuternodeSumTransaction<DIM>::computeOutgoingMessageSize()
{
   d_outgoing_bytes = 
      d_src_level->getPatch(d_src_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_src)
                 ->getDataStreamSize(*d_overlap);
   return(d_outgoing_bytes);
}

template<int DIM> int 
OuternodeSumTransaction<DIM>::getSourceProcessor()
{
   return(d_src_level->getMappingForPatch(d_src_patch));
}

template<int DIM> int 
OuternodeSumTransaction<DIM>::getDestinationProcessor()
{
   return(d_dst_level->getMappingForPatch(d_dst_patch));
}

template<int DIM> void 
OuternodeSumTransaction<DIM>::packStream(tbox::AbstractStream& stream)
{
   d_src_level->getPatch(d_src_patch)
              ->getPatchData(s_refine_items[d_refine_item_id]->
                             d_src)
              ->packStream(stream, *d_overlap);
}

template<int DIM> void 
OuternodeSumTransaction<DIM>::unpackStream(tbox::AbstractStream& stream)
{
   tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_dst_data =
      d_dst_level->getPatch(d_dst_patch)->
         getPatchData(s_refine_items[d_refine_item_id]->d_scratch);
   TBOX_ASSERT(!onode_dst_data.isNull());

   pdat::OuternodeData<DIM,double> onode_tmp_data(
      onode_dst_data->getBox(),
      onode_dst_data->getDepth() );
   onode_tmp_data.fillAll(0.0);
   onode_tmp_data.unpackStream(stream, *d_overlap);

   /*
    * Apply the sum - onode_dst_data += onode_tmp_data
    */
   sumData( *onode_dst_data, onode_tmp_data );
}

template<int DIM> void 
OuternodeSumTransaction<DIM>::copyLocalData()
{
   tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_dst_data =
      d_dst_level->getPatch(d_dst_patch)->
         getPatchData(s_refine_items[d_refine_item_id]->d_scratch);
   TBOX_ASSERT(!onode_dst_data.isNull());
 
   tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_src_data =
      d_src_level->getPatch(d_src_patch)->
         getPatchData(s_refine_items[d_refine_item_id]->d_src);
   TBOX_ASSERT(!onode_src_data.isNull());

   const hier::IntVector<DIM>& src_offset = d_overlap->getSourceOffset();
   hier::Box<DIM> shifted_src_box = onode_src_data->getBox();
   shifted_src_box.shift(src_offset);
   
   pdat::OuternodeData<DIM,double> onode_src_shift(shifted_src_box,
                                               onode_src_data->getDepth() );
   onode_src_shift.fillAll(0.0);
   onode_src_shift.copy(*onode_src_data, *d_overlap);
   
   /*
    * Apply the sum - onode_dst_data += onode_src_shift
    */
   sumData( *onode_dst_data, onode_src_shift);
}

/*
*************************************************************************
*                                                                       *
* Private method to sum data in the intersecting region of src (add)    *
* and dst outernode data objects.                                       *
*                                                                       *
*************************************************************************
*/ 
template<int DIM> void 
OuternodeSumTransaction<DIM>::sumData(
   pdat::OuternodeData<DIM,double> &dst,
   pdat::OuternodeData<DIM,double> &add) const
{
   TBOX_ASSERT( dst.getDepth() == add.getDepth() );

   const hier::Box<DIM> dst_node_box = pdat::NodeGeometry<DIM>::toNodeBox(dst.getBox());
   const hier::Box<DIM> add_node_box = pdat::NodeGeometry<DIM>::toNodeBox(add.getBox());
   const hier::Box<DIM> node_box_int = dst_node_box * add_node_box;
   
   /*
    * Compute overlapping regions between the outernode faces and
    * do an add operation on the corresponding array data.
    */
   math::ArrayDataBasicOps<DIM,double> mathops;

   for ( int dd = 0; dd < DIM; dd++ ) {
      for ( int ds = 0; ds < 2; ds++ ) {

         pdat::ArrayData<DIM,double> &dst_data_array =
            dst.getArrayData(dd, ds);

         if (dst_data_array.isInitialized()) {

            for (int ad = 0; ad < DIM; ad++ ) {
               for (int as = 0; as < 2; as++ ) {
                  
                  pdat::ArrayData<DIM,double> &add_data_array =
                     add.getArrayData(ad, as);
                  if (add_data_array.isInitialized()) {

                     mathops.add( dst_data_array,
                                  dst_data_array,
                                  add_data_array,
                                  node_box_int );
                  }
               }
            }
         }
      }
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
void OuternodeSumTransaction<DIM>::printClassData(ostream& stream) const
{
   stream << "Outernode Sum Transaction"                    << endl;
   stream << "   refine item array:        " 
          << (typename xfer::RefineClasses<DIM>::Data**)s_refine_items       
          << endl;
   stream << "   num refine items:      " << s_num_refine_items << endl;
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
          << (hier::PatchLevel<DIM>*)d_src_level               << endl;
   stream << "   source level:           " 
          << (hier::PatchLevel<DIM>*)d_src_level               << endl;
   stream << "   overlap:                "                     << endl;
   d_overlap->print(stream);
}

}
}
#endif

