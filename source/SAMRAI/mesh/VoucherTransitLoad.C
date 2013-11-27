/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
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

const int VoucherTransitLoad::VoucherTransitLoad_EDGETAG0;
const int VoucherTransitLoad::VoucherTransitLoad_EDGETAG1;

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
void VoucherTransitLoad::insertAll( const hier::BoxContainer &other )
{
   for ( hier::BoxContainer::const_iterator bi=other.begin(); bi!=other.end(); ++bi ) {
      d_voucher_set.insert( Voucher( LoadType(bi->size()), bi->getOwnerRank() ) );
      d_sumload += bi->size();
   }
}


/*
*************************************************************************
*************************************************************************
*/
void VoucherTransitLoad::insertAll( const TransitLoad &other_transit_load )
{
   const VoucherTransitLoad &other = recastTransitLoad(other_transit_load);
   for ( std::set<Voucher,VoucherRankCompare>::iterator si=d_voucher_set.begin();
         si!=d_voucher_set.end(); ++si ) {
      std::set<Voucher,VoucherRankCompare>::iterator sj =
         d_voucher_set.lower_bound(*si);
      if ( sj->d_issuer_rank == si->d_issuer_rank ) {
         Voucher combined_voucher( *si, *sj );
         d_voucher_set.erase(sj++);
         d_voucher_set.insert(combined_voucher);
      }
      else {
         d_voucher_set.insert(*si);
      }
   }
   d_sumload += other.d_sumload;
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
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level
 * and put edges in balanced<==>unbalanced Connector).
 *
 * This method does some P2P communication to convert the vouchers to
 * boxes then delegates to BoxTransitSet to do the rest.
 *
 * This method does two things:
 * - Voucher redeemers request and receive work for their vouchers.
 * - Voucher fulfillers receive and fulfill redemption requests.
 * The code is writen to let each process be both redeemers and
 * fulfillers.  Logic should drop through correctly when the process
 * plays just one role.
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


   // 1. Request work for vouchers to be redeemed.

   std::map<int,VoucherRedemption> redemptions_to_request;

   const hier::LocalId last_local_id = unbalanced_box_level.getLastLocalId();
   int count = 0;
   for ( const_iterator si=d_voucher_set.begin();
         si!=d_voucher_set.end(); ++si, ++count ) {

      hier::SequentialLocalIdGenerator id_gen(
         last_local_id + count,
         hier::LocalId(static_cast<int>(d_voucher_set.size())) );

      redemptions_to_request[si->d_issuer_rank].sendWorkDemand( *si, id_gen, mpi );
   }


   /*
    * If there is unaccounted work, then some process must have our
    * vouchers.  Fulfill these redemption requests until we accounted
    * for everything.
    *
    * For now, work is same as cell count.
    */

   LoadType unaccounted_work = LoadType(unbalanced_box_level.getLocalNumberOfCells())
      - findIssuerValue(mpi.getRank());

   if ( d_print_edge_steps ) {
      tbox::plog << "unaccounted work is " << unaccounted_work << std::endl;
   }

   // Set up the reserve for fulfilling incoming redemption requests.
   BoxTransitSet reserve(*d_pparams);
   reserve.insertAll( unbalanced_box_level.getBoxes() );


   // 2. Receive redemption requests for us to fulfill.

   std::map<int,VoucherRedemption> redemptions_to_fulfill;

   while ( unaccounted_work > d_pparams->getLoadComparisonTol() ) {

      tbox::SAMRAI_MPI::Status status;
      mpi.Probe( MPI_ANY_SOURCE, VoucherTransitLoad_EDGETAG0, &status );

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count( &status, MPI_CHAR, &count );

      VoucherRedemption &vr = redemptions_to_fulfill[source];
      vr.recvWorkDemand( source, count, mpi );

      unaccounted_work -= vr.d_voucher.d_load;
      TBOX_ASSERT( unaccounted_work > -d_pparams->getLoadComparisonTol() );

   }


   // 3. Fulfill redemption requests.

   for ( std::map<int,VoucherRedemption>::iterator mi=redemptions_to_fulfill.begin();
         mi!=redemptions_to_fulfill.end(); ++mi ) {

      VoucherRedemption &vr = mi->second;
      vr.sendWorkSupply( reserve,
                         unbalanced_to_balanced,
                         balanced_to_unbalanced,
                         *d_pparams );
   }


   // 4. Receive work for vouchers to be redeemed.

   while ( !redemptions_to_request.empty() ) {

      tbox::SAMRAI_MPI::Status status;
      mpi.Probe( MPI_ANY_SOURCE, VoucherTransitLoad_EDGETAG0, &status );

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count( &status, MPI_CHAR, &count );

      VoucherRedemption &vr = redemptions_to_request[source];
      vr.recvWorkSupply( count,
                         balanced_box_level,
                         unbalanced_to_balanced,
                         balanced_to_unbalanced,
                         *d_pparams );

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
   (*d_msg) << d_voucher.d_load << d_id_gen;
   d_mpi.Isend(
      (void*)(d_msg->getBufferStart()),
      static_cast<int>(d_msg->getCurrentSize()),
      MPI_CHAR,
      d_voucher.d_issuer_rank,
      VoucherTransitLoad_EDGETAG0,
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
      VoucherTransitLoad_EDGETAG0,
      &status );

   d_msg = boost::make_shared<tbox::MessageStream>(
      message_length, tbox::MessageStream::Read,
      static_cast<void*>(&incoming_message[0]), false);
   (*d_msg) >> d_voucher.d_load >> d_id_gen;
   TBOX_ASSERT( d_msg->endOfData() );
   TBOX_ASSERT( d_voucher.d_issuer_rank == mpi.getRank() );
   d_msg.reset();
}



/*
*************************************************************************
* Fulfill the voucher redemption by sending a supply of work to the
* demander.  The work is appropriated from a reserve.  Mapping edges
* based on the supplied work are saved in MappingConnectors.
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
   box_shipment.reassignOwnership( d_id_gen, d_mpi.getRank() );

   box_shipment.putToMessageStream(*d_msg);
   d_mpi.Isend(
      (void*)(d_msg->getBufferStart()),
      static_cast<int>(d_msg->getCurrentSize()),
      MPI_CHAR,
      d_demander_rank,
      VoucherTransitLoad_EDGETAG0,
      &d_mpi_request);

   box_shipment.generateLocalBasedMapEdges(
      unbalanced_to_balanced,
      balanced_to_unbalanced);

   return;
}



/*
*************************************************************************
* Receive work supply from the issuer of the voucher.  Mapping edges
* based on the work shipped are saved in MappingConnectors.
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
               VoucherTransitLoad_EDGETAG1,
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
 *
 * This method adjusts the load in this VoucherTransitLoad by
 * moving work between it (main_bin) and a holding_bin.  It tries to bring
 * main_bin's load to the specified ideal_load.
 *
 * The high_load and low_load define an acceptable range around the
 * ideal_load.  As soon as the main load falls in this range, no
 * further change is tried, even if it may bring the load closer to
 * the ideal.
 *
 * This method makes a best effort and returns the amount of load
 * moved.  It can move Vouchers between given sets and, if needed,
 * break some Vouchers up to move part of the work.
 *
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 * Return amount of load moved to main_bin from hold_bin.  Negative
 * amount means load moved from main_bin to hold_bin.
 *
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
   int issuer_rank;
   LoadType load;
   for (size_t i = 0; i < num_vouchers; ++i) {
      msg >> issuer_rank >> load;
      const Voucher voucher(load, issuer_rank);
      std::set<Voucher,VoucherRankCompare>::iterator sj =
         d_voucher_set.lower_bound(voucher);
      if ( sj->d_issuer_rank == voucher.d_issuer_rank ) {
         const Voucher combined_voucher( *sj, voucher );
         d_voucher_set.erase(sj++);
         d_voucher_set.insert(combined_voucher);
      }
      else {
         d_voucher_set.insert(voucher);
      }
   }
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
   co << getSumLoad() << " units in " << size() << " boxes.";
   if ( detail_depth > 0 ) {
      size_t count = 0;
      co << ":\n";
      for ( VoucherTransitLoad::const_iterator vi=begin();
            vi!=end() && count < 50; ++vi, ++count ) {
         tbox::plog << border << "    " << vi->d_issuer_rank << ':' << vi->d_load << '\n';
      }
   }
   else {
      co << ".\n";
   }
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
