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
 * If one half has a surplus and the other has a negative surplus,
 * move weight from the former to the latter.
 *************************************************************************
 */
void
CascadePartitioner::balanceHalves()
{

   if ( d_groups[inner_cycle].ourSurplus() >  d_pparams.getLoadComparisonTol() &&
        d_groups[inner_cycle].farSurplus() < -d_pparams.getLoadComparisonTol() ) {

      double weight = d_groups[inner_cycle].removeWeightFromOurHalf(
         -d_groups[inner_cycle].farSurplus() );

      d_groups[inner_cycle].addWeightToFarHalf(weight);

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

      double weight = d_groups[inner_cycle].removeWeightFromOtherHalf(
         -d_groups[inner_cycle].ourSurplus() );

      d_groups[inner_cycle].addWeightToOurHalf(weight);

      recv_comm.completeCurrentOperation();
      // Unpack into d_groups[0].d_local_load (local_load).
   }

   return;
}



/*
 *************************************************************************
 * 1. If group is single-process, try to remove load from d_local_load
 *    and put it in d_shipment.
 * 2. Else:
 *    a: Remove weight from the half matching the priority.
 *    b: If step a didn't remove enough, remove some from the other half.
 *************************************************************************
 */
double
CascadePartitioner::giveLoad( double amount, int taker )
{
   TBOX_ASSERT( amount > 0.0 );
   double removed = 0.0;
   if ( d_cycle_num == 0 ) {
      d_shipment = d_local_load->clone();
      shipment->adjustLoad( *d_local_load, amount, amount, amount );
      tbox::MessageStream msg;
      msg << *shipment;
      removed = d_shipment.getSumLoad();
      d_lower_weight -= removed;
   }

   const Position priority = taker < d_first_lower_rank ? Lower : Upper;

   else if ( priority == d_our_position ) {
      removed = giveLoadFromOurHalf( amount, priority );
      if ( removed < amount ) {
         removed += giveLoadFromFarHalf( amount-removed );
      }
   }
   else {
      removed = giveLoadFromFarHalf( amount );
      if ( removed < amount ) {
         removed += giveLoadFromOurHalf( amount-removed, taker );
      }
   }
   return removed;
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
