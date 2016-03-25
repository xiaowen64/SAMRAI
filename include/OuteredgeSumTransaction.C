//
// File:        OuteredgeSumTransaction.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 693 $
// Modified:    $Date: 2005-10-28 13:39:23 -0700 (Fri, 28 Oct 2005) $
// Description: Communication transaction for summing outeredge data
//
 
#ifndef included_algs_OuteredgeSumTransaction_C
#define included_algs_OuteredgeSumTransaction_C

#include "OuteredgeSumTransaction.h"

#include "Patch.h"
#include "PatchData.h"
#include "ArrayDataBasicOps.h"
#include "EdgeGeometry.h"
#include "OuteredgeData.h"
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
   OuteredgeSumTransaction<DIM>::s_refine_items =
      (const typename xfer::RefineClasses<DIM>::Data**)NULL;
template<int DIM> int OuteredgeSumTransaction<DIM>::s_num_refine_items = 0;

template<int DIM> void OuteredgeSumTransaction<DIM>::setRefineItems(
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

template<int DIM> void OuteredgeSumTransaction<DIM>::unsetRefineItems()
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

template<int DIM> OuteredgeSumTransaction<DIM>::OuteredgeSumTransaction(
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

 
template<int DIM> OuteredgeSumTransaction<DIM>::~OuteredgeSumTransaction()
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
OuteredgeSumTransaction<DIM>::canEstimateIncomingMessageSize()
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
OuteredgeSumTransaction<DIM>::computeIncomingMessageSize()
{
   d_incoming_bytes = 
      d_dst_level->getPatch(d_dst_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_scratch)
                 ->getDataStreamSize(*d_overlap);
   return(d_incoming_bytes);
}

template<int DIM> int 
OuteredgeSumTransaction<DIM>::computeOutgoingMessageSize()
{
   d_outgoing_bytes = 
      d_src_level->getPatch(d_src_patch)
                 ->getPatchData(s_refine_items[d_refine_item_id]->
                                d_src)
                 ->getDataStreamSize(*d_overlap);
   return(d_outgoing_bytes);
}

template<int DIM> int 
OuteredgeSumTransaction<DIM>::getSourceProcessor()
{
   return(d_src_level->getMappingForPatch(d_src_patch));
}

template<int DIM> int 
OuteredgeSumTransaction<DIM>::getDestinationProcessor()
{
   return(d_dst_level->getMappingForPatch(d_dst_patch));
}

template<int DIM> void 
OuteredgeSumTransaction<DIM>::packStream(tbox::AbstractStream& stream)
{
   d_src_level->getPatch(d_src_patch)
              ->getPatchData(s_refine_items[d_refine_item_id]->
                             d_src)
              ->packStream(stream, *d_overlap);
}

template<int DIM> void 
OuteredgeSumTransaction<DIM>::unpackStream(tbox::AbstractStream& stream)
{
   tbox::Pointer<pdat::OuteredgeData<DIM,double> > oedge_dst_data = 
      d_dst_level->getPatch(d_dst_patch)->
         getPatchData(s_refine_items[d_refine_item_id]->d_scratch);
   TBOX_ASSERT(!oedge_dst_data.isNull());

   pdat::OuteredgeData<DIM,double> oedge_tmp_data(
      oedge_dst_data->getBox(),
      oedge_dst_data->getDepth() );
   oedge_tmp_data.fillAll(0.0);
   oedge_tmp_data.unpackStream(stream, *d_overlap);

   /*
    * Apply the sum - oedge_dst_data += oedge_tmp_data
    */
   sumData( *oedge_dst_data, oedge_tmp_data );
}

template<int DIM> void 
OuteredgeSumTransaction<DIM>::copyLocalData()
{
   tbox::Pointer< pdat::OuteredgeData<DIM,double> > oedge_dst_data = 
      d_dst_level->getPatch(d_dst_patch)->
         getPatchData(s_refine_items[d_refine_item_id]->d_scratch);
   TBOX_ASSERT(!oedge_dst_data.isNull());

   tbox::Pointer< pdat::OuteredgeData<DIM,double> > oedge_src_data = 
      d_src_level->getPatch(d_src_patch)->
         getPatchData(s_refine_items[d_refine_item_id]->d_src);
   TBOX_ASSERT(!oedge_src_data.isNull());

   const hier::IntVector<DIM>& src_offset = d_overlap->getSourceOffset();
   hier::Box<DIM> shifted_src_box = oedge_src_data->getBox();
   shifted_src_box.shift(src_offset);
   
   pdat::OuteredgeData<DIM,double> oedge_src_shift(shifted_src_box,
                                                   oedge_src_data->getDepth() );
   oedge_src_shift.fillAll(0.0);
   oedge_src_shift.copy(*oedge_src_data, *d_overlap);

   /*
    * Apply the sum - oedge_dst_data += oedge_src_shift
    */
   sumData( *oedge_dst_data, oedge_src_shift );

}

template<int DIM> void OuteredgeSumTransaction<DIM>::sumData(
   pdat::OuteredgeData<DIM,double> &dst,
   pdat::OuteredgeData<DIM,double> &add ) const
{
   TBOX_ASSERT( dst.getDepth() == add.getDepth() );

   /*
    * Compute overlapping regions between the outeredge faces and
    * do an add operation on the corresponding array data.
    */
   math::ArrayDataBasicOps<DIM,double> mathops;

   // a = axis
   for ( int a = 0; a < DIM; a++ ) {

      const hier::Box<DIM> dst_edge_box = 
         pdat::EdgeGeometry<DIM>::toEdgeBox(dst.getBox(),a);
      const hier::Box<DIM> add_edge_box = 
         pdat::EdgeGeometry<DIM>::toEdgeBox(add.getBox(),a);
      const hier::Box<DIM> edge_box_intsct = 
         dst_edge_box * add_edge_box;

      // df = dest face normal
      for ( int df = 0; df < DIM; df++ ) {
         if (df != a) {
            // ds = dest upper/lower side
            for ( int ds = 0; ds < 2; ds++ ) {
               pdat::ArrayData<DIM,double> &dst_data =
                  dst.getArrayData(a,df,ds);

               // sf = source face normal
               for ( int sf = 0; sf < DIM; sf++ ) {
                  if (sf != a) {
                     // ss = source upper/lower side
                     for ( int ss = 0; ss < 2; ss++ ) {
                        pdat::ArrayData<DIM,double> &add_data =
                           add.getArrayData(a,sf,ss);
               
                        // dst_data += add_data
                        if (dst_data.isInitialized() && 
                            add_data.isInitialized()) {
                           mathops.add( dst_data,
                                        dst_data,
                                        add_data,
                                        edge_box_intsct );
                        }
                     } // loop over ss
                  } // sf != a
               } // loop over sf
                     
            } // loop over df
         } // df != a   
      } // loop over df

   } // loop over a

}

/*
*************************************************************************
*                                                                       *
* Function to print state of transaction.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
OuteredgeSumTransaction<DIM>::printClassData(ostream& stream) const
{
   stream << "Outeredge Sum Transaction"                        << endl;
   stream << "   refine item array:        " 
          << (typename xfer::RefineClasses<DIM>::Data**)s_refine_items       
          << endl;
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
          << (hier::PatchLevel<DIM>*)d_src_level               << endl;
   stream << "   source level:           "
          << (hier::PatchLevel<DIM>*)d_src_level               << endl;
   stream << "   overlap:                "                     << endl;
   d_overlap->print(stream);
}

}
}
#endif
