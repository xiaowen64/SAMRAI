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
   d_vouchers(),
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
      d_vouchers[bi->getOwnerRank()] += LoadType(bi->size());
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
   for ( std::map<int,LoadType>::const_iterator mi=other.d_vouchers.begin();
         mi!=other.d_vouchers.end(); ++mi ) {
      d_vouchers[mi->first] += mi->second;
   }
   d_sumload += other.d_sumload;
}


/*
*************************************************************************
*************************************************************************
*/
size_t VoucherTransitLoad::getNumberOfItems() const
{
   return d_vouchers.size();
}


/*
*************************************************************************
*************************************************************************
*/
size_t VoucherTransitLoad::getNumberOfOriginatingProcesses() const
{
   return d_vouchers.size();
}


/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level
 * and put edges in balanced<==>unbalanced Connector).
 *
 * We can generate balanced--->unbalanced edges for all boxes because
 * we have their origin info.  If the box originated locally, we can
 * generate the unbalanced--->balanced edge for them as well.
 * However, we can't generate these edges for boxes originating
 * remotely.  They are generated in
 * constructSemilocalUnbalancedToBalanced, which uses communication.
 */
void
VoucherTransitLoad::assignContentToLocalProcessAndGenerateMap(
   hier::BoxLevel& balanced_box_level,
   hier::MappingConnector &balanced_to_unbalanced,
   hier::MappingConnector &unbalanced_to_balanced ) const
{
   d_object_timers->t_assign_content_to_local_process_and_generate_map->start();

   if ( d_print_steps ) {
      tbox::plog << "VoucherTransitLoad::assignUnassignedToLocalProcessAndGenerateMap: entered." << std::endl;
   }

   BoxTransitSet transit_boxes(*d_pparams);

   TBOX_ERROR("Missing code to convert vouchers to fill transit_boxes.");

   transit_boxes.assignContentToLocalProcessAndGenerateMap(
      balanced_box_level, balanced_to_unbalanced, unbalanced_to_balanced );

   if ( d_print_steps ) {
      tbox::plog << "VoucherTransitLoad::assignUnassignedToLocalProcessAndGenerateMap: exiting." << std::endl;
   }

   d_object_timers->t_assign_content_to_local_process_and_generate_map->stop();
}





/*
 *************************************************************************
 * Construct semilocal relationships in unbalanced--->balanced
 * Connector.
 *
 * Determine relationships in unbalanced_to_balanced by sending
 * balanced boxes back to the owners of the unbalanced Boxes that
 * originated them.  We don't know what ranks will sending to us, so
 * we keep receiving messages from any rank until we have accounted
 * for all the cells in the unbalanced BoxLevel.
 *************************************************************************
 */
void
VoucherTransitLoad::constructSemilocalUnbalancedToBalanced(
   hier::MappingConnector &unbalanced_to_balanced,
   const VoucherTransitLoad &kept_imports ) const
{
   d_object_timers->t_construct_semilocal->start();

   if ( d_print_steps ) {
      tbox::plog << "VoucherTransitLoad::constructSemilocalUnbalancedToBalanced: entered."
                 << std::endl;
   }

   const hier::BoxLevel &balanced_box_level = unbalanced_to_balanced.getHead();
   const hier::BoxLevel &unbalanced_box_level = unbalanced_to_balanced.getBase();
   const tbox::SAMRAI_MPI &mpi = unbalanced_box_level.getMPI();

#if 1
   TBOX_ERROR("Unfinished code.");
#else
   // Stuff the imported boxes into buffers by their original owners.
   d_object_timers->t_pack_edge->start();
   std::map<int,boost::shared_ptr<tbox::MessageStream> > outgoing_messages;
   for ( const_iterator bi=kept_imports.begin(); bi!=kept_imports.end(); ++bi ) {
      const Voucher &bit = *bi;
      boost::shared_ptr<tbox::MessageStream> &mstream =
         outgoing_messages[bit.d_orig_box.getOwnerRank()];
      if ( !mstream ) {
         mstream.reset(new tbox::MessageStream);
      }
      bit.putToMessageStream(*mstream);
   }
   d_object_timers->t_pack_edge->stop();


   /*
    * Send outgoing_messages.  Optimization for mitigating contention:
    * Start by sending to the first recipient with a rank higher than
    * the local rank.
    */

   std::map<int,boost::shared_ptr<tbox::MessageStream> >::iterator recip_itr =
      outgoing_messages.upper_bound(mpi.getRank());
   if ( recip_itr == outgoing_messages.end() ) {
      recip_itr = outgoing_messages.begin();
   }

   int outgoing_messages_size = static_cast<int>(outgoing_messages.size());
   std::vector<tbox::SAMRAI_MPI::Request>
      send_requests( outgoing_messages_size, MPI_REQUEST_NULL );

   d_object_timers->t_construct_semilocal_send_edges->start();
   for ( int send_number = 0; send_number < outgoing_messages_size; ++send_number ) {

      int recipient = recip_itr->first;
      tbox::MessageStream &mstream = *recip_itr->second;

      if ( d_print_edge_steps ) {
         tbox::plog << "Accounting for cells on proc " << recipient << '\n';
      }

      mpi.Isend(
         (void*)(mstream.getBufferStart()),
         static_cast<int>(mstream.getCurrentSize()),
         MPI_CHAR,
         recipient,
         VoucherTransitLoad_EDGETAG0,
         &send_requests[send_number]);

      ++recip_itr;
      if ( recip_itr == outgoing_messages.end() ) {
         recip_itr = outgoing_messages.begin();
      }

   }
   d_object_timers->t_construct_semilocal_send_edges->stop();


   int num_cells_imported = 0;
   for ( const_iterator si=kept_imports.begin(); si!=kept_imports.end(); ++si ) {
      num_cells_imported += si->d_box.size();
   }

   int num_unaccounted_cells = static_cast<int>(
      unbalanced_box_level.getLocalNumberOfCells() + num_cells_imported
      - balanced_box_level.getLocalNumberOfCells() );

   if ( d_print_edge_steps ) {
      tbox::plog << num_unaccounted_cells << " unaccounted cells." << std::endl;
   }


   /*
    * Receive info about exported cells from processes that now own
    * those cells.  Receive until all cells are accounted for.
    */

   std::vector<char> incoming_message;
   Voucher balanced_box_in_transit(unbalanced_box_level.getDim());

   while ( num_unaccounted_cells > 0 ) {

      d_object_timers->t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Status status;
      mpi.Probe( MPI_ANY_SOURCE, VoucherTransitLoad_EDGETAG0, &status );

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count( &status, MPI_CHAR, &count );
      incoming_message.resize( count, -1 );

      mpi.Recv(
         static_cast<void*>(&incoming_message[0]),
         count,
         MPI_CHAR,
         source,
         VoucherTransitLoad_EDGETAG0,
         &status );
      d_object_timers->t_construct_semilocal_comm_wait->stop();

      tbox::MessageStream msg( incoming_message.size(),
                               tbox::MessageStream::Read,
                               static_cast<void*>(&incoming_message[0]),
                               false );
      const int old_count = num_unaccounted_cells;
      d_object_timers->t_unpack_edge->start();
      while ( !msg.endOfData() ) {

         balanced_box_in_transit.getFromMessageStream(msg);
         unbalanced_to_balanced.insertLocalNeighbor(
            balanced_box_in_transit.d_box,
            balanced_box_in_transit.d_orig_box.getBoxId() );
         num_unaccounted_cells -= balanced_box_in_transit.d_box.size();

      }
      d_object_timers->t_unpack_edge->stop();

      if ( d_print_edge_steps ) {
         tbox::plog << "Process " << source << " accounted for "
                    << (old_count-num_unaccounted_cells) << " cells, leaving "
                    << num_unaccounted_cells << " unaccounted.\n";
      }

      incoming_message.clear();
   }
   TBOX_ASSERT( num_unaccounted_cells == 0 );


   // Wait for the sends to complete before clearing outgoing_messages.
   if (send_requests.size() > 0) {
      std::vector<tbox::SAMRAI_MPI::Status> status(send_requests.size());
      d_object_timers->t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Waitall(
         static_cast<int>(send_requests.size()),
         &send_requests[0],
         &status[0]);
      d_object_timers->t_construct_semilocal_comm_wait->stop();
      outgoing_messages.clear();
   }
#endif

   if ( d_print_steps ) {
      tbox::plog << "VoucherTransitLoad::constructSemilocalUnbalancedToBalanced: exiting."
                 << std::endl;
   }

   d_object_timers->t_construct_semilocal->stop();

   return;
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
      getTimer(timer_prefix + "::assignContentToLocalProcessAndGenerateMap()");
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::putToMessageStream( tbox::MessageStream &msg ) const
{
   msg << static_cast<int>(d_vouchers.size());
   for (const_iterator ni = d_vouchers.begin(); ni != d_vouchers.end(); ++ni) {
      msg << ni->first << ni->second;
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
   int num_vouchers = 0;
   msg >> num_vouchers;
   int orig_rank;
   LoadType load;
   for (int i = 0; i < num_vouchers; ++i) {
      msg >> orig_rank >> load;
      d_vouchers[orig_rank] += load;
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
      for ( VoucherTransitLoad::const_iterator bi=begin();
            bi!=end() && count < 50; ++bi, ++count ) {
         tbox::plog << border << "    " << bi->first << ':' << bi->second << '\n';
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
