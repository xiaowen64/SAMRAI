/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TransitLoad by vouchers.
 *
 ************************************************************************/

#ifndef included_mesh_VoucherTransitLoad_C
#define included_mesh_VoucherTransitLoad_C

#include "SAMRAI/mesh/VoucherTransitLoad.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int VoucherTransitLoad::VoucherTransitLoad_DEMANDTAG;
const int VoucherTransitLoad::VoucherTransitLoad_SUPPLYTAG;

const std::string VoucherTransitLoad::s_default_timer_prefix("mesh::VoucherTransitLoad");
std::map<std::string, VoucherTransitLoad::TimerStruct> VoucherTransitLoad::s_static_timers;

tbox::StartupShutdownManager::Handler
VoucherTransitLoad::s_initialize_finalize_handler(
   VoucherTransitLoad::initializeCallback,
   0,
   0,
   VoucherTransitLoad::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);


/*
*************************************************************************
*************************************************************************
*/
VoucherTransitLoad::VoucherTransitLoad( const PartitioningParams &pparams ) :
   TransitLoad(),
   d_voucher_set(),
   d_sumload(0),
   d_pparams(&pparams),
   d_print_steps(false),
   d_print_edge_steps(false)
{
   setTimerPrefix(s_default_timer_prefix);
}


/*
*************************************************************************
*************************************************************************
*/
VoucherTransitLoad::VoucherTransitLoad( const VoucherTransitLoad &other ) :
   TransitLoad(other),
   d_voucher_set(other.d_voucher_set),
   d_sumload(other.d_sumload),
   d_pparams(other.d_pparams),
   d_print_steps(other.d_print_steps),
   d_print_edge_steps(other.d_print_edge_steps),
   d_object_timers(other.d_object_timers)
{
}


/*
*************************************************************************
Initialize sets to a new (empty) container but retains current
supplemental data such as diagnostic parameters.
*************************************************************************
*/
void VoucherTransitLoad::initialize()
{
   d_voucher_set.clear();
   d_sumload = 0.0;
}


/*
*************************************************************************
Initialize sets to a new (empty) container but retains current
supplemental data such as diagnostic parameters.
*************************************************************************
*/
boost::shared_ptr<TransitLoad> VoucherTransitLoad::clone() const
{
   boost::shared_ptr<VoucherTransitLoad> new_object =
      boost::make_shared<VoucherTransitLoad>(*this);
   new_object->setAllowBoxBreaking(getAllowBoxBreaking());
   new_object->initialize();
   return new_object;
}


/*
*************************************************************************
*************************************************************************
*/
void VoucherTransitLoad::insertAll( const hier::BoxContainer &other )
{
   for ( hier::BoxContainer::const_iterator bi=other.begin(); bi!=other.end(); ++bi ) {
      insertCombine( Voucher( LoadType(bi->size()), bi->getOwnerRank() ) );
   }
}


/*
*************************************************************************
*************************************************************************
*/
void VoucherTransitLoad::insertAll( const TransitLoad &other_transit_load )
{
   const VoucherTransitLoad &other = recastTransitLoad(other_transit_load);
   for ( const_iterator si=other.d_voucher_set.begin();
         si!=other.d_voucher_set.end(); ++si ) {
      insertCombine(*si);
   }
}


/*
*************************************************************************
*************************************************************************
*/
size_t VoucherTransitLoad::getNumberOfItems() const
{
   return d_voucher_set.size();
}


/*
*************************************************************************
*************************************************************************
*/
size_t VoucherTransitLoad::getNumberOfOriginatingProcesses() const
{
   return d_voucher_set.size();
}


/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::putToMessageStream( tbox::MessageStream &msg ) const
{
   msg << d_voucher_set.size();
   for (const_iterator ni = d_voucher_set.begin(); ni != d_voucher_set.end(); ++ni) {
      msg << ni->d_issuer_rank << ni->d_load;
   }
}


/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::getFromMessageStream( tbox::MessageStream &msg )
{
   /*
    * As we pull each Voucher out, give it a new id that reflects
    * its new owner.
    */
   size_t num_vouchers = 0;
   msg >> num_vouchers;
   Voucher v;
   for (size_t i = 0; i < num_vouchers; ++i) {
      msg >> v.d_issuer_rank >> v.d_load;
      insert(v);
   }
}


/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level
 * and populate balanced<==>unbalanced).
 *
 * This method does two things:
 * - Voucher redeemers request and receive work vouchers they hold.
 * - Voucher fulfillers receive and fulfill redemption requests.
 * The code is writen to let each process be both redeemers and
 * fulfillers.  Logic should drop through correctly on processes
 * that plays just one role.
 *
 * There are four major steps, organized to overlap communication.
 * 1. Request work for vouchers to be redeemed.
 * 2. Receive redemption requests.
 * 3. Fulfill redemption requests.
 * 4. Receive work for redeemed vouchers.
 */
void
VoucherTransitLoad::assignContentToLocalProcessAndPopulateMaps(
   hier::BoxLevel& balanced_box_level,
   hier::MappingConnector &balanced_to_unbalanced,
   hier::MappingConnector &unbalanced_to_balanced )
{
   d_object_timers->t_assign_content_to_local_process_and_generate_map->start();

   if ( d_print_steps ) {
      tbox::plog << "VoucherTransitLoad::assignContentToLocalProcessAndPopulateMaps: entered." << std::endl;
   }

   const hier::BoxLevel &unbalanced_box_level = unbalanced_to_balanced.getBase();
   const tbox::SAMRAI_MPI &mpi = unbalanced_box_level.getMPI();


   /*
    * unaccounted_work is amount of work we started with, minus what
    * we can still acount for (the part we still hold).  The rest have
    * been sent off in the form of vouchers and we don't know where
    * they ended up.
    */
   LoadType unaccounted_work = LoadType(unbalanced_box_level.getLocalNumberOfCells())
      - yankVoucher(mpi.getRank()).d_load;

   if ( d_print_edge_steps ) {
      tbox::plog << "unaccounted work is " << unaccounted_work << std::endl;
   }


   // 1. Send work demands for vouchers we want to redeem.

   std::map<int,VoucherRedemption> redemptions_to_request;

   const hier::LocalId last_local_id = unbalanced_box_level.getLastLocalId();
   const hier::LocalId local_id_inc(1 + d_voucher_set.size());

   int count = 1;
   for ( const_iterator si=d_voucher_set.begin();
         si!=d_voucher_set.end(); ++si, ++count ) {

      hier::SequentialLocalIdGenerator id_gen( last_local_id + count, local_id_inc );

      redemptions_to_request[si->d_issuer_rank].sendWorkDemand( *si, id_gen, mpi );
   }

   // Set up the reserve for fulfilling incoming redemption requests.
   BoxTransitSet reserve(*d_pparams);
   reserve.insertAll( unbalanced_box_level.getBoxes() );


   // 2. Receive work demands for voucher we generated but can't account for.

   /*
    * If there is unaccounted work, then some process must have our
    * voucher.  Reveive their demand for work until we've accounted
    * for everything.  Don't supply work until all demands are
    * received, because that leads to non-deterministic results.
    */

   std::map<int,VoucherRedemption> redemptions_to_fulfill;

   while ( unaccounted_work > d_pparams->getLoadComparisonTol() ) {

      tbox::SAMRAI_MPI::Status status;
      mpi.Probe( MPI_ANY_SOURCE, VoucherTransitLoad_DEMANDTAG, &status );

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count( &status, MPI_CHAR, &count );

      VoucherRedemption &vr = redemptions_to_fulfill[source];
      vr.recvWorkDemand( source, count, mpi );

      unaccounted_work -= vr.d_voucher.d_load;
      TBOX_ASSERT( unaccounted_work > -d_pparams->getLoadComparisonTol() );

   }


   // 3. Supply work according to received demands.

   for ( std::map<int,VoucherRedemption>::iterator mi=redemptions_to_fulfill.begin();
         mi!=redemptions_to_fulfill.end(); ++mi ) {

      VoucherRedemption &vr = mi->second;
      vr.sendWorkSupply( reserve,
                         unbalanced_to_balanced,
                         balanced_to_unbalanced,
                         *d_pparams );
   }


   // Anything left in reserve is kept locally.
   hier::SequentialLocalIdGenerator id_gen( last_local_id, local_id_inc );
   reserve.reassignOwnership(
      id_gen,
      balanced_box_level.getMPI().getRank() );

   reserve.putInBoxLevel(balanced_box_level);
   reserve.generateLocalBasedMapEdges(
      unbalanced_to_balanced,
      balanced_to_unbalanced);


   // 4. Receive work according to the demands we sent.

   while ( !redemptions_to_request.empty() ) {

      tbox::SAMRAI_MPI::Status status;
      mpi.Probe( MPI_ANY_SOURCE, VoucherTransitLoad_SUPPLYTAG, &status );

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count( &status, MPI_CHAR, &count );

      VoucherRedemption &vr = redemptions_to_request[source];
      vr.recvWorkSupply( count,
                         balanced_box_level,
                         unbalanced_to_balanced,
                         balanced_to_unbalanced,
                         *d_pparams );

      redemptions_to_request.erase(source);
   }


   if ( d_print_steps ) {
      tbox::plog << "VoucherTransitLoad::assignContentToLocalProcessAndPopulateMaps: exiting." << std::endl;
   }

   d_object_timers->t_assign_content_to_local_process_and_generate_map->stop();
}



/*
*************************************************************************
* Start redeeming a voucher by sending the demand for work to the
* voucher issuer.
*************************************************************************
*/
void VoucherTransitLoad::VoucherRedemption::sendWorkDemand(
   const Voucher &voucher,
   const hier::SequentialLocalIdGenerator &id_gen,
   const tbox::SAMRAI_MPI &mpi )
{
   d_voucher = voucher;
   d_id_gen = id_gen;
   d_mpi = mpi;

   d_msg = boost::make_shared<tbox::MessageStream>();
   (*d_msg) << d_id_gen << d_voucher;

   d_mpi.Isend(
      (void*)(d_msg->getBufferStart()),
      static_cast<int>(d_msg->getCurrentSize()),
      MPI_CHAR,
      d_voucher.d_issuer_rank,
      VoucherTransitLoad_DEMANDTAG,
      &d_mpi_request);

   d_msg.reset();
}



/*
*************************************************************************
*************************************************************************
*/
void VoucherTransitLoad::VoucherRedemption::recvWorkDemand(
   int demander_rank,
   int message_length,
   const tbox::SAMRAI_MPI &mpi )
{
   d_mpi = mpi;
   d_demander_rank = demander_rank;

   std::vector<char> incoming_message(message_length);
   tbox::SAMRAI_MPI::Status status;
   d_mpi.Recv(
      static_cast<void*>(&incoming_message[0]),
      message_length,
      MPI_CHAR,
      d_demander_rank,
      VoucherTransitLoad_DEMANDTAG,
      &status );

   d_msg = boost::make_shared<tbox::MessageStream>(
      message_length, tbox::MessageStream::Read,
      static_cast<void*>(&incoming_message[0]), false);
   (*d_msg) >> d_id_gen >> d_voucher;
   TBOX_ASSERT( d_msg->endOfData() );
   TBOX_ASSERT( d_voucher.d_issuer_rank == mpi.getRank() );

   d_msg.reset();
}



/*
*************************************************************************
* Fulfill the voucher redemption by sending a supply of work to the
* demander.  The work is appropriated from a reserve.  Save mapping
* edges incident from local boxes.
*************************************************************************
*/
void VoucherTransitLoad::VoucherRedemption::sendWorkSupply(
   BoxTransitSet &reserve,
   hier::MappingConnector &unbalanced_to_balanced,
   hier::MappingConnector &balanced_to_unbalanced,
   const PartitioningParams &pparams )
{
   BoxTransitSet box_shipment(pparams);
   box_shipment.adjustLoad( reserve,
                            d_voucher.d_load,
                            d_voucher.d_load,
                            d_voucher.d_load );
   box_shipment.reassignOwnership( d_id_gen, d_demander_rank );

   d_msg = boost::make_shared<tbox::MessageStream>();
   box_shipment.putToMessageStream(*d_msg);

   d_mpi.Isend(
      (void*)(d_msg->getBufferStart()),
      static_cast<int>(d_msg->getCurrentSize()),
      MPI_CHAR,
      d_demander_rank,
      VoucherTransitLoad_SUPPLYTAG,
      &d_mpi_request);

   box_shipment.generateLocalBasedMapEdges(
      unbalanced_to_balanced,
      balanced_to_unbalanced);

   return;
}



/*
*************************************************************************
* Receive work supply from the issuer of the voucher.  Save mapping
* edges incident from local boxes.
*************************************************************************
*/
void VoucherTransitLoad::VoucherRedemption::recvWorkSupply(
   int message_length,
   hier::BoxLevel &balanced_box_level,
   hier::MappingConnector &unbalanced_to_balanced,
   hier::MappingConnector &balanced_to_unbalanced,
   const PartitioningParams &pparams )
{
   std::vector<char> incoming_message(message_length);
   tbox::SAMRAI_MPI::Status status;
   d_mpi.Recv( static_cast<void*>(&incoming_message[0]),
               message_length,
               MPI_CHAR,
               d_voucher.d_issuer_rank,
               VoucherTransitLoad_SUPPLYTAG,
               &status );

   d_msg = boost::make_shared<tbox::MessageStream>(
      message_length, tbox::MessageStream::Read,
      static_cast<void*>(&incoming_message[0]), false);
   BoxTransitSet box_shipment(pparams);
   box_shipment.getFromMessageStream(*d_msg);
   d_msg.reset();

   box_shipment.putInBoxLevel(balanced_box_level);
   box_shipment.generateLocalBasedMapEdges(
      unbalanced_to_balanced,
      balanced_to_unbalanced);
}



/*
*************************************************************************
*************************************************************************
*/
void VoucherTransitLoad::VoucherRedemption::finishSendRequest()
{
   if ( d_mpi_request != MPI_REQUEST_NULL ) {
      tbox::SAMRAI_MPI::Wait( &d_mpi_request, MPI_STATUS_IGNORE );
   }
   TBOX_ASSERT( d_mpi_request == MPI_REQUEST_NULL );
}




/*
*************************************************************************
* Adjust the VoucherTransitLoad by moving work between it (main_bin)
* and a hold_bin.  Try to bring the load to the specified ideal.
*
* The high_load and low_load define an acceptable range around the
* ideal_load.
*
* This method makes a best effort and returns the amount of load
* moved.  It can move Vouchers between given bins and, if needed,
* break some Vouchers up to move part of the work.
*
* This method is purely local--it reassigns the load but does not
* communicate the change to any remote process.
*
* Return amount of load moved to this object from hold_bin.  Negative
* amount means load moved from this object to hold_bin.
*************************************************************************
*/
VoucherTransitLoad::LoadType
VoucherTransitLoad::adjustLoad(
   TransitLoad& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load )
{
   VoucherTransitLoad& main_bin(*this);

   if (d_print_steps) {
      tbox::plog << "  adjustLoad attempting to bring main load from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " or within [" << low_load << ", " << high_load << "]."
                 << std::endl;
   }
   TBOX_ASSERT( low_load <= ideal_load );
   TBOX_ASSERT( high_load >= ideal_load );


   LoadType actual_transfer = 0;

   if ((main_bin.empty() && ideal_load <= 0 ) ||
       (hold_bin.empty() && main_bin.getSumLoad() < ideal_load )) {
      return actual_transfer;
   }

   d_object_timers->t_adjust_load->start();

   const LoadType change = ideal_load - main_bin.getSumLoad();

   if ( change > 0 ) {
      // Move load to main_bin.
      actual_transfer =
         raiseDstLoad( recastTransitLoad(hold_bin),
                       main_bin,
                       ideal_load );
   }
   else if ( change < 0 ) {
      // Move load to hold_bin.
      actual_transfer =
         -raiseDstLoad( main_bin,
                        recastTransitLoad(hold_bin),
                        hold_bin.getSumLoad() - change );
   }

   if ( d_print_steps ) {
      const LoadType point_miss = main_bin.getSumLoad() - ideal_load;
      const LoadType range_miss =
         main_bin.getSumLoad() > high_load ? main_bin.getSumLoad() - high_load :
         main_bin.getSumLoad() < low_load ? low_load - main_bin.getSumLoad() : 0;
      tbox::plog << "  adjustLoad point_miss=" << point_miss
                 << "  range_miss="
                 << (range_miss > 0 ? " ":"") // Add space if missed range
                 << (range_miss > 0.5*d_pparams->getMinBoxSize().getProduct() ? " ":"") // Add space if missed range by a lot
                 << range_miss
                 << "  " << main_bin.getSumLoad() << '/'
                 << ideal_load << " [" << low_load << ',' << high_load << ']'
                 << std::endl;
   }

   d_object_timers->t_adjust_load->stop();

   return actual_transfer;
}




/*
*************************************************************************
* Raise load of dst container by shifing load from src.
*************************************************************************
*/
VoucherTransitLoad::LoadType
VoucherTransitLoad::raiseDstLoad(
   VoucherTransitLoad& src,
   VoucherTransitLoad& dst,
   LoadType ideal_dst_load )
{
   TBOX_ASSERT( ideal_dst_load >= dst.getSumLoad() );

   if ( src.empty() ) {
      return 0;
      // No-op empty-container cases is not handled below.
   }

src.recursivePrint(tbox::plog, "src: ", 1);
dst.recursivePrint(tbox::plog, "dst: ", 1);

   /*
    * Decide whether to take work from the beginning or the end of
    * src, whichever is closer to the dst.  This is a minor
    * optimization to better preserve locality when the issuer rank
    * ranges of src and dst don't overlap much.
    */
   bool take_from_src_end = true;
   if ( !dst.empty() ) {

      const LoadType gap_at_src_end =
         tbox::MathUtilities<double>::Abs(
            src.rbegin()->d_issuer_rank - dst.begin()->d_issuer_rank);

      const LoadType gap_at_src_begin =
         tbox::MathUtilities<double>::Abs(
            dst.rbegin()->d_issuer_rank - src.begin()->d_issuer_rank);

      take_from_src_end = gap_at_src_end < gap_at_src_begin;
   }

   LoadType old_dst_load = dst.getSumLoad();

   while ( dst.getSumLoad() < ideal_dst_load && !src.empty() ) {

      iterator src_itr;
      if ( take_from_src_end ) {
         src_itr = src.end();
         --src_itr;
      }
      else {
         src_itr = src.begin();
      }

      Voucher free_voucher = *src_itr;
      src.erase(src_itr);

      if ( free_voucher.d_load < (ideal_dst_load - dst.getSumLoad()) ) {
         dst.insert(free_voucher);
      }
      else {
         Voucher partial_voucher((ideal_dst_load - dst.getSumLoad()), free_voucher);
         dst.insert(partial_voucher);
         src.insert(free_voucher);
      }

   }

   return (dst.getSumLoad() - old_dst_load);
}



/*
 ***********************************************************************
 ***********************************************************************
 */
bool
VoucherTransitLoad::eraseIssuer( int issuer_rank )
{
   Voucher tmp_voucher( 0.0, issuer_rank );
   const_iterator vi = d_voucher_set.lower_bound(tmp_voucher);
   if ( vi != d_voucher_set.end() && vi->d_issuer_rank == issuer_rank ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      const_iterator vi1 = vi;
      ++vi1;
      TBOX_ASSERT( vi1 == d_voucher_set.end() || vi1->d_issuer_rank != issuer_rank);
#endif
      d_sumload -= vi->d_load;
      d_voucher_set.erase(vi);
      return true;
   }
   return false;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
VoucherTransitLoad::Voucher
VoucherTransitLoad::yankVoucher( int issuer_rank )
{
   Voucher tmp_voucher( 0.0, issuer_rank );
   const_iterator vi = d_voucher_set.lower_bound(tmp_voucher);
   if ( vi != d_voucher_set.end() && vi->d_issuer_rank == issuer_rank ) {
      tmp_voucher.d_load = vi->d_load;
      erase(vi);
   }
   return tmp_voucher;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix));
   if (ti == s_static_timers.end()) {
      d_object_timers = &s_static_timers[timer_prefix];
      getAllTimers(timer_prefix, *d_object_timers);
   } else {
      d_object_timers = &(ti->second);
   }
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_adjust_load = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::adjustLoad()");

   timers.t_assign_content_to_local_process_and_generate_map = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::assignContentToLocalProcessAndPopulateMaps()");
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::recursivePrint(
   std::ostream &co,
   const std::string &border,
   int detail_depth ) const
{
   co << getSumLoad() << " units in " << size() << " vouchers";
   if ( detail_depth > 0 ) {
      size_t count = 0;
      co << ":\n";
      for ( VoucherTransitLoad::const_iterator vi=begin();
            vi!=end() && count < 50; ++vi, ++count ) {
         co << border << "    " << *vi << '\n';
      }
   }
   else {
      co << ".\n";
   }
   co.flush();
}


}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
