//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/datamovers/standard/RefineTimeTransaction.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Communication transaction for time interpolation during data refining
//
 
#ifndef included_xfer_RefineTimeTransaction
#define included_xfer_RefineTimeTransaction

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_hier_BoxOverlap
#include "BoxOverlap.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
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
#ifndef included_xfer_TimeInterpolateOperator
#include "TimeInterpolateOperator.h"
#endif

namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Class RefineTimeTransaction<DIM> represents a single time interpolation 
 * communication transaction between two processors or a local data copy or refine 
 * schedules.  Note that to there is an implicit hand-shaking between objects of 
 * this class and the RefineSchedule<DIM> object that constructs them.  Following 
 * the refine schedule implementation, the source patch data indices for a time 
 * transaction are always refer to the old and new source data and the destination 
 * patch data index for a time transaction is always the scratch data, all as defined 
 * in the RefineClasses<DIM> class.  This transaction is used by the refine schedule.
 *
 * @see xfer::RefineSchedule
 * @see xfer::RefineClasses
 * @see tbox::Schedule
 * @see tbox::Transaction
 */

template<int DIM> class RefineTimeTransaction : public tbox::Transaction
{
public:
   /*!
    * Static member function to set the array of refine class data items that
    * is shared by all object instances of this time transaction class during 
    * data transfers.  The array must be set before any transactions are executed.  
    * The array is set in the RefineSchedule<DIM> class.
    */
   static void setRefineItems(const typename RefineClasses<DIM>::Data** refine_items,
                              int num_refine_items);

   /*!
    * Static member function to unset the array of refine class data items that
    * is shared by all object instances of this time transaction class during
    * data transfers.  The unset function is used to prevent erroneous execution 
    * of different schedules.  The array is unset in the RefineSchedule<DIM> class.
    */
   static void unsetRefineItems();

   /*!
    * Static member function to set the transaction time that will be shared 
    * by all object instances of this time transaction class during time 
    * interpolation.  This transaction time must be set before any transactions 
    * are executed.
    */
   static void setTransactionTime(const double time);

   /*!
    * Construct a transaction with the specified source and destination levels,
    * patches, and patch data components found in the refine class
    * item with the given id owned by the calling refine schedule.  In general,
    * this constructor is called by a RefineSchedule<DIM> object for each data
    * transaction (involving time interpolation) that must occur.  This transaction 
    * will be responsible for one of the following: (1) performing a local copy, 
    * (2) packing a message stream from the source, or (3) unpacking a message stream 
    * from the destination.  The transaction will perform time interpolation between
    * the source old and new times using the time interpolation operator found in
    * the refine class item.
    * 
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch      Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch      Integer index of source patch in source patch level.
    * @param box            hier::Box region in which to time interpolate.
    * @param refine_item_id   Integer id of refine data item owned by refine schedule.
    *
    * When assertion checking is active, an assertion will result if any of the pointer
    * arguments is null, or if any of the integer arguments are invalid (i.e., < 0);
    */
   RefineTimeTransaction(
      tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
      tbox::Pointer< hier::PatchLevel<DIM> > src_level,
      tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
      int dst_patch,
      int src_patch,
      const hier::Box<DIM>& box,
      int refine_item_id);

   /*!
    * The virtual destructor for time transaction releases all
    * memory associated with the transaction.
    */
   virtual ~RefineTimeTransaction<DIM>();

   /*!
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communications protocol kicks in and the message size is transmitted
    * between nodes.
    */
   virtual bool canEstimateIncomingMessageSize();

   /*!
    * Return the integer amount of buffer space (in bytes) needed for the 
    * incoming message.  This routine is only called if the transaction 
    * can estimate the size of the incoming message.
    */
   virtual int computeIncomingMessageSize();

   /*!
    * Return the integer buffer space needed (in bytes) for the outgoing message.
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
   virtual void printClassData(std::ostream& stream) const;

private:
   RefineTimeTransaction(const RefineTimeTransaction<DIM>&); // not implemented
   void operator=(const RefineTimeTransaction<DIM>&);	// not implemented

   static double s_time;
   static const typename RefineClasses<DIM>::Data** s_refine_items;
   static int s_num_refine_items;

   void timeInterpolate(
      const tbox::Pointer< hier::PatchData<DIM> >& pd_dst,
      const tbox::Pointer< hier::PatchData<DIM> >& pd_old,
      const tbox::Pointer< hier::PatchData<DIM> >& pd_new);

   tbox::Pointer< hier::PatchLevel<DIM> >    d_dst_level;
   tbox::Pointer< hier::PatchLevel<DIM> >    d_src_level;
   tbox::Pointer< hier::BoxOverlap<DIM> >    d_overlap;
   int                               d_dst_patch;
   int                               d_src_patch;
   hier::Box<DIM>                         d_box;
   int                               d_refine_item_id;
   int                               d_incoming_bytes;
   int                               d_outgoing_bytes;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "RefineTimeTransaction.C"
#endif
