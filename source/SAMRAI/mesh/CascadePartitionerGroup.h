/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner
#define included_mesh_CascadePartitioner

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace mesh {


/*!
 * @brief A grouping of processes in the CascadePartitioner algorithm.
 *
 * Smallest gruops are single-process groups.  With each
 * CascadePartitioner cycle, two adjacent groups are combined to make
 * a bigger group.  The bigger group can shift load from its
 * overloaded half to its underloaded half.
 *
 */
class CascadePartitionerGroup {

   //! @brief Where a group falls in the next larger group.
   enum Position { Lower=0, Upper=1 };

   CascadePartitionerGroup() :
      d_mpi(),
      d_cycle_num(-1),
      d_group_size(-1),
      d_first_lower_rank(-1),
      d_first_upper_rank(-1),
      d_end_rank(-1),
      d_contact(-1),
      d_our_half(0),
      d_our_position(Lower),
      d_lower_weight(0.0),
      d_upper_weight(0.0),
      d_our_weight(0.0),
      d_far_weight(0.0),
      d_may_supply_from_our_half(true),
      d_may_supply_from_far_half(true),
      d_local_load(0) {}
      d_shipment(0) {}

   ~CascadePartitionerGroup() {}

   /*!
    * @brief Make a cycle-zero, single-process group.
    */
   void makeSingleProcessGroup( tbox::SAMRAI_MPI &mpi, TransitLoad &local_load,
                                TransitLoad &shipment ) {
      d_mpi = mpi;
      d_cycle_num = 0;
      d_group_size = 1;
      d_first_lower_rank = d_mpi.getRank();
      d_first_upper_rank = d_first_lower_rank+1;
      d_end_rank = d_first_lower_rank+1;
      d_contact = -1;
      d_our_half = 0;
      d_our_position = Lower;
      d_lower_weight = local_load.getSumLoad();
      d_upper_weight = 0.0;
      d_our_weight = &d_lower_weight;
      d_far_weight = 0;
      d_local_load = &local_load;
      d_shipment = &shipment;
   }

   /*!
    * @brief Make a combo group consistng of the given half-group and
    * its corresponding partner, which this method will determine.
    *
    * Local process provides the half containing the local rank
    * (our_half).  Data for other half (the far half) is obtained by
    * communication.
    */
   void makeComboGroup( CascadePartitionerGroup &our_half ) {

      /*
       * Set up the groups
       */

      d_our_half = &our_half;
      d_cycle_num = our_half.d_cycle_num+1;
      d_mpi = our_half.d_mpi;
      d_local_load = 0; // For single-process groups only.

      d_group_size = 1 << d_cycle_num;
      int group_num = d_mpi.getRank()/group_size;

      d_first_lower_rank = group_size*group_num;
      d_first_upper_rank = tbox::MathUtilities<int>::Min(
         first_lower_rank + group_size/2, d_mpi.getSize());
      d_end_rank = tbox::MathUtilities<int>::Min(
         d_first_lower_rnak + group_size, d_mpi.getSize());

      int relative_rank = d_mpi.getRank() - d_first_lower_rank;

      d_our_position = relative_rank >= group_size/2 ? Upper : Lower;
      d_contact = d_position == Lower ?
         d_mpi.getRank() + group_size/2 :
         d_mpi.getRank() - group_size/2;
      if ( d_contact >= d_mpi.getSize() ) {
         d_contact = -1;
      }
      d_our_weight = d_our_position == Lower ? &d_lower_weight : &d_upper_weight;
      d_far_weight = d_our_position == Lower ? &d_upper_weight : &d_lower_weight;

      /*
       * Record weights of the two halves.
       */
      double far_weight = 0.0;
      double our_weight = d_our_half->getComboWeight();
      if ( d_contact >= 0 ) {
         d_mpi.Sendrecv( &our_weight, 1, MPI_DOUBLE, d_contact, 0,
                         &far_weight, 1, MPI_DOUBLE, d_contact, 1 );
      }
      d_lower_weight = d_our_position == Lower ? our_weight : far_weight;
      d_upper_weight = d_our_position == Upper ? our_weight : far_weight;
   }


private:

   /*!
    * @brief Improve balance of the two halves of this group by moving
    * load from overloaded half to underloadef half.
    */
   void balanceConstituentHalves();

   //! @brief Position of the half containing local process.
   Position ourPosition() const { return d_our_position; }
   //! @brief Position of the half not containing local process.
   Position farPosition() const { return !d_our_position; }
   //! @brief Contact rank.
   int contact() const { return d_contact; }
   //! @brief Weight of the group (combined weight of its two halves).
   double getComboWeight() const { return d_lower_weight + d_upper_weight; }
   //! @brief Surplus of lower half.
   double lowerSurplus() const {
      return d_lower_weight - d_global_load_avg*(d_first_upper_rank-d_first_lower_rank);
   }
   //! @brief Surplus of upper half.
   double upperSurplus() const {
      return d_upper_weight - d_global_load_avg*(d_end_rank-d_first_upper_rank);
   }
   //! @brief Surplus of our half.
   double ourSurplus() const {
      return d_our_position == Lower ? getLowerSurplus() : getUpperSurplus();
   }
   //! @brief Surplus of far half.
   double farSurplus() const {
      return d_our_position == Lower ? getUpperSurplus() : getLowerSurplus();
   }

   /*!
    * @brief Try to supply an amount of weight by removing it from the
    * group, and return the (estimated) amount supplied.
    *
    * Removing weight from groups containing more than the local process
    * reports an estimate of the amount removed, based on what should have
    * been removed.  Due to load cutting restrictions, the actual amount
    * removed may differ.
    *
    * @param taker Representative of the group demanding this work.
    * Local process is to send to this representative any work it
    * personally supplies.
    */
   double supplyLoad( double amount, int taker );

   /*!
    * @brief Try to remove an amount of weight from our half of the group.
    *
    * The return value is exact if the group includes only the
    * local rank, otherwise, it's an estimate based on what the
    * requested removal amount was.
    *
    * This method is recursive by the sequence
    * supplyLoadFromOurHalf-supplyLoad-supplyLoadFromOurHalf.
    *
    * @param taker Rank of process taking load from local process.
    *
    * @return Amount removed (or an estimate)
    */
   double supplyLoadFromOurHalf( double amount, int taker ) {
      TBOX_ASSERT( amount > 0.0 );
      double removed = 0.0;
      if ( d_may_supply_from_our_half &&
           ourSurplus() >= d_pparams->getLoadComparisonTol() ) {

         removed = d_our_half->supplyLoad( amount, taker_rank );
         if ( d_our_position == Lower ) {
            d_lower_weight -= removed;
         } else {
            d_upper_weight -= removed;
         }
      }
      return removed;
   }

   /*!
    * @brief Symbolically try to remove an amount of weight from the
    * half of the group not containing the local process.
    *
    * No real work is exchanged because the local process is not in
    * the far half.  This method just estimates what the far half
    * would give away.
    *
    * @return Estimate of amount removed (actual value not locally
    * available)
    */
   double supplyLoadFromFarHalf( double amount ) {
      TBOX_ASSERT( amount > 0.0 );
      double removed = 0.0;
      if ( d_may_supply_from_far_half &&
           farSurplus() >= d_pparams->getLoadComparisonTol() ) {
         removed = tbox::MathUtilities<double>::Min( amount, d_far_weight );
         d_far_weight -= removed;
      }
      return removed;
   }

   /*!
    * @brief Record work amount received (estimated) by our half and
    * forbid it from becoming a supplier.
    */
   void recordDemandReceivedByOurHalf( double amount ) {
      if ( d_our_position == Lower ) { d_lower_weight += amount; }
      else { d_upper_weight += amount; }
      d_may_supply_from_our_half = false;
   }
   /*!
    * @brief Record work amount received (estimated) by far half and
    * forbid it from becoming a supplier.
    */
   void recordDemandReceivedByFarHalf( double amount ) {
      if ( d_our_position == Upper ) { d_lower_weight += amount; }
      else { d_upper_weight += amount; }
      d_may_supply_from_far_half = false;
   }

   void sendMyShipment( int taker );
   void unpackSuppliedLoad();

   tbox::SAMRAI_MPI &d_mpi;

   //! @brief Cycle number.  Group has 2^d_cycle_num ranks.
   int d_cycle_num;

   //! @brief Number of processes in the group.
   int d_group_size;

   //! @brief First rank in lower half.
   int d_first_lower_rank;

   //! @brief First rank in upper half.
   int d_first_upper_rank;

   //! @brief One past the last rank.
   int d_end_rank;

   //! Rank of contact (if any) in the far half of the group.
   int d_contact;

   //! @brief The half containing me.
   CascadePartitionerGroup *d_our_half;

   //! @brief Position of d_our_half in this group.
   Position d_our_position;

   //! @brief Sum of load held by lower half (or approxmimation).
   double d_lower_weight;

   //! @brief Sum of load held by upper half (or approxmimation).
   double d_upper_weight;

   //! @brief Points to either d_lower_weight or d_upper_weight.
   double *d_our_weight;

   //! @brief Points to either d_lower_weight or d_upper_weight.
   double *d_far_weight;

   //! @brief Whether our half may supply load.
   bool d_may_supply_from_our_half;

   //! @brief Whether far half may supply load.
   bool d_may_supply_from_far_half;


   //! @brief Load of local process, for single-process group.
   TransitLoad *d_local_load;
   //! @brief Cache for to be shipped, for single-process group.
   bool::shared_ptr<TransitLoad> d_shipment;
   tbox::AsyncCommPeer<char> d_comm;
};

}
}

#endif
