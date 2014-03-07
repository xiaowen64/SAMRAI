/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner_C
#define included_mesh_CascadePartitioner_C

#include "SAMRAI/mesh/CascadePartitioner.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {



/*
 *************************************************************************
 * If one half has a positive surplus and the other half has a
 * negative surplus, the former supplies work to the latter.  Amount
 * supplied is ideally the minimum of the supplier's surplus and the
 * demander's deficit.  (Actual ammounts are limited by load cutting
 * restrictions.)
 *
 * This method records estimates of the weight changes to the groups
 * it knows about, It doesn't record the actual weight changes because
 * that happens remotely.  Each process on the supply half send a
 * message received by its contact on the demand half.  The messages
 * has the actual work to be transfered.
 *************************************************************************
 */
void
CascadePartitioner::balanceConstituentHalves()
{

   if ( d_groups[inner_cycle].ourSurplus() >  d_pparams.getLoadComparisonTol() &&
        d_groups[inner_cycle].farSurplus() < -d_pparams.getLoadComparisonTol() ) {

      double supplied_weight = d_groups[inner_cycle].supplyLoadFromOurHalf(
         -d_groups[inner_cycle].farSurplus(), d_contact );

      d_groups[inner_cycle].addWeightToFarHalf(supplied_weight);

      /*
       * If local indicated to contact that local has surplus,
       * send the shipment set aside by groups[0].
       */
      if ( d_groups[0].canSend() ) {
         send_comm.beginSend();
      }
   }

   else if ( d_groups[inner_cycle].farSurplus() >  d_pparams.getLoadComparisonTol() &&
             d_groups[inner_cycle].ourSurplus() < -d_pparams.getLoadComparisonTol() ) {

      if ( d_contact_has_surplus ) {
         recv_comm.beginRecv( d_contact );
      }

      double supplied_weight = d_groups[inner_cycle].supplyLoadFromFarHalf(
         -d_groups[inner_cycle].ourSurplus() );

      d_groups[inner_cycle].addWeightToOurHalf(supplied_weight);

      recv_comm.completeCurrentOperation();
      // Unpack into d_groups[0].d_local_load (local_load).
   }

   return;
}



/*
 *************************************************************************
 * Supply specified amount of load to another group with a work
 * demand.  Any load supplied by local process is to be sent to the
 * designated taker, a process in the demand group.  Give priority
 * to supplies closest to the taker in rank space.
 *
 * 1. If group is single-process, try to remove load from d_local_load
 *    and put it in d_shipment.
 * 2. Else:
 *    a: Supply load from the half closest to the taker.
 *    b: If step a didn't supply enough, remove some from the other half.
 *
 * This method is recursive by the sequence
 * supplyLoad-supplyLoadFromOurHalf-supplyLoad.
 *************************************************************************
 */
double
CascadePartitioner::supplyLoad( double amount, int taker )
{
   TBOX_ASSERT( amount > 0.0 );
   double removed = 0.0;

   if ( d_cycle_num == 0 ) {
      // Group is single-process, not two halves.
      d_shipment = d_local_load->clone();
      d_shipment->adjustLoad( *d_local_load, amount, amount, amount );
      removed = d_shipment->getSumLoad();
      d_lower_weight -= removed;
      sendMyShipment(taker);
   }

   else if ( ( d_our_position == Lower && taker <  d_first_lower_rank ) ||
             ( d_our_position == upper && taker >= d_first_upper_rank ) ) {
      // Group is two-halves, with our half closer to taker.
      removed = supplyLoadFromOurHalf( amount, priority );
      if ( removed < amount ) {
         removed += supplyLoadFromFarHalf( amount-removed );
      }
   }
   else {
      // Group is two-halves, with far half closer to taker.
      removed = supplyLoadFromFarHalf( amount );
      if ( removed < amount ) {
         removed += supplyLoadFromOurHalf( amount-removed, taker );
      }
   }
   return removed;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitioner::sendMyShipment( int taker )
{
   tbox::MessageStream msg;
   msg << *shipment;
   d_comm.setPeerRank(taker);
   d_comm.beginSend( msg.getBufferStart(), msg.getCurrentSize() );
   return;
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
