/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for summing outernode data
 *
 ************************************************************************/

#ifndef included_algs_OuternodeSumTransaction
#define included_algs_OuternodeSumTransaction

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineClasses.h"

#include <iostream>

namespace SAMRAI {
namespace algs {

/*!
 * @brief Class OuternodeSumTransaction represents a single outernode data sum
 * communication transaction between two processors or a local data sum for
 * refine schedules.  Note that to there is an implicit hand-shaking between
 * objects of this class and the xfer::RefineSchedule object that constructs
 * them.  Following the refine schedule implementation, the source patch data
 * index for a transaction always refers to the source data and the destination
 * patch data index for a transaction is always the scratch data, all as
 * defined in the xfer::RefineClasses class.
 *
 * @see xfer::RefineSchedule
 * @see xfer::RefineClasses
 * @see tbox::Schedule
 * @see tbox::Transaction
 */

class OuternodeSumTransaction:public tbox::Transaction
{
public:
   /*!
    * Static member function to set the array of refine class data items that
    * is shared by all object instances of this sum transaction class during
    * data transfers.  The array must be set before any transactions are
    * executed.  The array is set in the xfer::RefineSchedule class.
    *
    * @pre refine_items != 0);
    * @pre num_refine_items >= 0
    */
   static void
   setRefineItems(
      const xfer::RefineClasses::Data *const* refine_items,
      int num_refine_items)
   {
      TBOX_ASSERT(refine_items != 0);
      TBOX_ASSERT(num_refine_items >= 0);
      s_refine_items = refine_items;
      s_num_refine_items = num_refine_items;
   }

   /*!
    * Static member function to unset the array of refine class data items that
    * is shared by all object instances of this sum transaction class during
    * data transfers.  The unset function is used to prevent erroneous
    * execution of different schedules.  The array is unset in the
    * RefineSchedule class.
    */
   static void
   unsetRefineItems()
   {
      s_refine_items = 0;
      s_num_refine_items = 0;
   }

   /*!
    * Construct a transaction with the specified source and destination
    * levels, patches, and patch data components found in the refine class
    * item with the given id owned by the calling refine schedule.  In
    * general, this constructor is called by a xfer::RefineSchedule
    * object for each data transaction (specifically summing outernode
    * data) that must occur.  This transaction will be responsible for one
    * of the following: (1) a local data copy and sum, or (2) packing a
    * message stream with source patch data, or (3) unpacking and summing
    * destination patch data from a message stream.
    *
    * @param dst_level        boost::shared_ptr to destination patch level.
    * @param src_level        boost::shared_ptr to source patch level.
    * @param overlap          boost::shared_ptr to overlap region between
    *                         patches.
    * @param dst_node         Destination Box in destination patch level.
    * @param src_node         Source Box in source patch level.
    * @param refine_item_id   Integer id of refine data item owned by refine
    *                         schedule.
    *
    * @pre dst_level
    * @pre src_level
    * @pre overlap
    * @pre dst_node.getLocalId() >= 0
    * @pre src_node.getLocalId() >= 0
    * @pre refine_item_id >= 0
    * @pre (dst_level->getDim() == src_level->getDim()) &&
    *      (dst_level->getDim() == dst_node.getDim()) &&
    *      (dst_level->getDim() == src_node.getDim())
    */
   OuternodeSumTransaction(
      const boost::shared_ptr<hier::PatchLevel>& dst_level,
      const boost::shared_ptr<hier::PatchLevel>& src_level,
      const boost::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_node,
      const hier::Box& src_node,
      int refine_item_id);

   /*!
    * The virtual destructor for the copy transaction releases all
    * memory associated with the transaction.
    */
   virtual ~OuternodeSumTransaction();

   /*!
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communication protocol kicks in and the message size is transmitted
    * between nodes.
    */
   virtual bool
   canEstimateIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the incoming
    * message.  This routine is only called if the transaction can estimate the
    * size of the incoming message.  See canEstimateIncomingMessageSize().
    */
   virtual size_t
   computeIncomingMessageSize();

   /*!
    * Return the integer buffer space (in bytes) needed for the outgoing
    * message.
    */
   virtual size_t
   computeOutgoingMessageSize();

   /*!
    * Return the sending processor number for the communications transaction.
    */
   virtual int
   getSourceProcessor();

   /*!
    * Return the receiving processor number for the communications transaction.
    */
   virtual int
   getDestinationProcessor();

   /*!
    * Pack the transaction data into the message stream.
    */
   virtual void
   packStream(
      tbox::MessageStream& stream);

   /*!
    * Unpack the transaction data from the message stream.
    *
    * @pre d_dst_level->getPatch(d_dst_node.getGlobalId())->getPatchData(s_refine_items[d_refine_item_id]->d_scratch) is actually a boost::shared_ptr<pdat::OuternodeData<double> >
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream);

   /*!
    * Perform the local data copy for the transaction.
    *
    * @pre d_dst_level->getPatch(d_dst_node.getGlobalId())->getPatchData(s_refine_items[d_refine_item_id]->d_scratch) is actually a boost::shared_ptr<pdat::OuternodeData<double> >
    * @pre d_src_level->getPatch(d_src_node.getGlobalId())->getPatchData(s_refine_items[d_refine_item_id]->d_src) is actually a boost::shared_ptr<pdat::OuternodeData<double> >
    */
   virtual void
   copyLocalData();

   /*!
    * Print out transaction information.
    */
   virtual void
   printClassData(
      std::ostream& stream) const;

private:
   OuternodeSumTransaction(
      const OuternodeSumTransaction&);                      // not implemented
   void
   operator = (
      const OuternodeSumTransaction&);             // not implemented

   static const xfer::RefineClasses::Data*const* s_refine_items;
   static int s_num_refine_items;

   boost::shared_ptr<hier::PatchLevel> d_dst_level;
   boost::shared_ptr<hier::PatchLevel> d_src_level;
   boost::shared_ptr<hier::BoxOverlap> d_overlap;
   hier::Box d_dst_node;
   hier::Box d_src_node;
   int d_refine_item_id;
   int d_incoming_bytes;
   int d_outgoing_bytes;

};

}
}

#endif
