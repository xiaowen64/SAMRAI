//
// File:        OuteredgeSumTransaction.h
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 684 $
// Modified:    $Date: 2005-10-21 14:59:38 -0700 (Fri, 21 Oct 2005) $
// Description: Communication transaction for summing outeredge data
//
 
#ifndef included_algs_OuteredgeSumTransaction
#define included_algs_OuteredgeSumTransaction

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_hier_BoxOverlap
#include "BoxOverlap.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_pdat_OuteredgeData
#include "OuteredgeData.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Transaction
#include "tbox/Transaction.h"
#endif
#ifndef included_xfer_RefineClasses
#include "RefineClasses.h"
#endif


namespace SAMRAI {
    namespace algs {

/*!
 * @brief Class OuteredgeSumTransaction<DIM> represents a single outeredge data sum 
 * communication transaction between two processors or a local data sum for refine schedules.
 * Note that to there is an implicit hand-shaking between objects of this class and
 * the xfer::RefineSchedule<DIM> object that constructs them.  Following the refine
 * schedule implementation, the source patch data index for a transaction always refers
 * to the source data and the destination patch data index for a transaction is always
 * the scratch data, all as defined in the xfer::RefineClasses<DIM> class.
 *
 * @see xfer::RefineSchedule
 * @see xfer::RefineClasses
 * @see tbox::Schedule
 * @see tbox::Transaction
 */

template<int DIM> class OuteredgeSumTransaction : public tbox::Transaction
{
public:
   /*!
    * Static member function to set the array of refine class data items that
    * is shared by all object instances of this sum transaction class during
    * data transfers.  The array must be set before any transactions are executed.
    * The array is set in the RefineSchedule<DIM> class.
    */
   static void setRefineItems(
      const typename xfer::RefineClasses<DIM>::Data** refine_items,
      int num_refine_items);

   /*!
    * Static member function to unset the array of refine class data items that
    * is shared by all object instances of this sum transaction class during
    * data transfers.  The unset function is used to prevent erroneous execution
    * of different schedules.  The array is unset in the RefineSchedule<DIM> class.
    */
   static void unsetRefineItems();

   /*!
    * Construct a transaction with the specified source and destination
    * levels, patches, and patch data components found in the refine class
    * item with the given id owned by the calling refine schedule.  In 
    * general, this constructor is called by a xfer::RefineSchedule<DIM> 
    * object for each data transaction (specifically summing outeredge 
    * data) that must occur.  This transaction will be responsible for one 
    * of the following: (1) a local data copy and sum, or (2) packing a 
    * message stream with source patch data, or (3) unpacking and summing 
    * destination patch data from a message stream.
    *
    * @param dst_level        tbox::Pointer to destination patch level.
    * @param src_level        tbox::Pointer to source patch level.
    * @param overlap          tbox::Pointer to overlap region between patches.
    * @param dst_patch        Integer index of destination patch in 
    *                         destination 
    *                         patch level. 
    * @param src_patch        Integer index of source patch in source patch 
    *                         level. 
    * @param refine_item_id   Integer id of refine data item owned by refine 
    *                         schedule.
    *
    * When assertion checking is active, an assertion will result if any of the pointer
    * arguments is null, or if any of the integer arguments is invalid (i.e., < 0).
    */
   OuteredgeSumTransaction(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
      int dst_patch,
      int src_patch,
      int refine_item_id);

   /*!
    * The virtual destructor for the copy transaction releases all
    * memory associated with the transaction.
    */
   virtual ~OuteredgeSumTransaction<DIM>();

   /*!
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communication protocol kicks in and the message size is transmitted
    * between edges.
    */
   virtual bool canEstimateIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the incoming 
    * message. This routine is only called if the transaction can estimate the
    * size of the incoming message.  See canEstimateIncomingMessageSize().
    */
   virtual int computeIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the outgoing message.
    */
   virtual int computeOutgoingMessageSize();

   /*!
    * Return the sending processor number for the communications transaction.
    */
   virtual int getSourceProcessor();

   /*!
    * Return the receiving processor number for the communications transaction.
    */
   virtual int getDestinationProcessor();

   /*!
    * Pack the transaction data into the message stream.
    */
   virtual void packStream(tbox::AbstractStream& stream);

   /*!
    * Unpack the transaction data from the message stream.
    */
   virtual void unpackStream(tbox::AbstractStream& stream);

   /*!
    * Perform the local data copy for the transaction.
    */
   virtual void copyLocalData();

   /*!
    * Print out transaction information.
    */
   virtual void printClassData(ostream& stream) const;

private:
   OuteredgeSumTransaction<DIM>(const OuteredgeSumTransaction<DIM>&); // not implemented
   void operator=(const OuteredgeSumTransaction<DIM>&); // not implemented

   static const typename xfer::RefineClasses<DIM>::Data** s_refine_items;
   static int s_num_refine_items;

   void sumData(pdat::OuteredgeData<DIM,double> &dst,
                pdat::OuteredgeData<DIM,double> &add) const;

   tbox::Pointer< hier::PatchLevel<DIM> >    d_dst_level;
   tbox::Pointer< hier::PatchLevel<DIM> >    d_src_level;
   tbox::Pointer< hier::BoxOverlap<DIM> >    d_overlap;
   int                               d_dst_patch;
   int                               d_src_patch;
   int                               d_refine_item_id;
   int                               d_incoming_bytes;
   int                               d_outgoing_bytes;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuteredgeSumTransaction.C"
#endif






