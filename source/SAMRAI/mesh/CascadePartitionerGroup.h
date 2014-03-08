/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitionerGroup
#define included_mesh_CascadePartitionerGroup

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"
#include "boost/shared_ptr.hpp"

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

public:

   CascadePartitionerGroup() :
      d_mpi(tbox::SAMRAI_MPI::getSAMRAIWorld()),
      d_cycle_num(-1),
      d_first_lower_rank(-1),
      d_first_upper_rank(-1),
      d_end_rank(-1),
      d_contact(-1),
      d_our_half(0),
      d_our_position(Lower),
      d_lower_weight(0.0),
      d_upper_weight(0.0),
      d_our_weight(0),
      d_far_weight(0),
      d_our_half_may_supply(true),
      d_far_half_may_supply(true),
      d_local_load(0),
      d_global_load_avg(0),
      d_pparams(0),
      d_shipment(),
      d_comm() {}

   /*!
    * @brief Copy constructor doesn't copy anything.  It is not used but is required
    * for the class to be used in an stl::vector.
    */
   CascadePartitionerGroup( const CascadePartitionerGroup &other ) :
      d_mpi(MPI_COMM_NULL),
      d_cycle_num(-1),
      d_first_lower_rank(-1),
      d_first_upper_rank(-1),
      d_end_rank(-1),
      d_contact(-1),
      d_our_half(0),
      d_our_position(Lower),
      d_lower_weight(0.0),
      d_upper_weight(0.0),
      d_our_weight(0),
      d_far_weight(0),
      d_our_half_may_supply(true),
      d_far_half_may_supply(true),
      d_local_load(0),
      d_global_load_avg(0),
      d_pparams(0),
      d_shipment(),
      d_comm() {}

   ~CascadePartitionerGroup() {}

   /*!
    * @brief Make a cycle-zero, single-process group.
    */
   void makeSingleProcessGroup(
      tbox::SAMRAI_MPI &mpi,
      const PartitioningParams &pparams,
      TransitLoad &local_load,
      double global_load );

   /*!
    * @brief Make a combo group consistng of the given half-group and
    * its corresponding partner, which this method will determine.
    */
   void makeComboGroup( CascadePartitionerGroup &our_half );

   /*!
    * @brief Improve balance of the two halves of this group by moving
    * load from overloaded half to underloadef half.
    */
   void balanceConstituentHalves();


private:

   //! @brief Where a group falls in the next larger group.
   enum Position { Lower=0, Upper=1 };

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
      return d_our_position == Lower ? lowerSurplus() : upperSurplus();
   }
   //! @brief Surplus of far half.
   double farSurplus() const {
      return d_our_position == Lower ? upperSurplus() : lowerSurplus();
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
      if ( d_our_half_may_supply &&
           ourSurplus() >= d_pparams->getLoadComparisonTol() ) {

         removed = d_our_half->supplyLoad( amount, taker );
         *d_our_weight -= removed;
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
      if ( d_far_half_may_supply &&
           farSurplus() >= d_pparams->getLoadComparisonTol() ) {
         removed = tbox::MathUtilities<double>::Min( amount, *d_far_weight );
         *d_far_weight -= removed;
      }
      return removed;
   }

   /*!
    * @brief Record estimated work amount received by our half-group and
    * that the half-group may not become a supplier.
    */
   void recordDemandReceivedByOurHalf( double amount ) {
      if ( d_our_position == Lower ) { d_lower_weight += amount; }
      else { d_upper_weight += amount; }
      d_our_half_may_supply = false;
   }
   /*!
    * @brief Record estimated work amount received by far half-group and
    * that the half-group may not become a supplier.
    */
   void recordDemandReceivedByFarHalf( double amount ) {
      if ( d_our_position == Upper ) { d_lower_weight += amount; }
      else { d_upper_weight += amount; }
      d_far_half_may_supply = false;
   }

   void sendMyShipment( int taker );
   void unpackSuppliedLoad();

   tbox::SAMRAI_MPI d_mpi;

   //! @brief Cycle number.  Group has 2^d_cycle_num ranks.
   int d_cycle_num;

   //! @brief First rank in lower half.
   int d_first_lower_rank;

   //! @brief First rank in upper half.
   int d_first_upper_rank;

   //! @brief One past the last rank.
   int d_end_rank;

   //! Rank of contact (if any) in the far half of the group.
   int d_contact;

   //! Whether contact may supply load.
   bool d_contact_may_supply;

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
   bool d_our_half_may_supply;

   //! @brief Whether far half may supply load.
   bool d_far_half_may_supply;


   //! @brief Load of local process, for single-process group.
   TransitLoad *d_local_load;

   double d_global_load_avg;

   const PartitioningParams *d_pparams;

   //! @brief Cache for to be shipped, for single-process group.
   boost::shared_ptr<TransitLoad> d_shipment;
   tbox::AsyncCommPeer<char> d_comm;
};

}
}

#endif
