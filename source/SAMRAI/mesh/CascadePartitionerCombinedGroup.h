/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitionerCombinedGroup
#define included_mesh_CascadePartitionerCombinedGroup

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"
#include "boost/shared_ptr.hpp"

namespace SAMRAI {
namespace mesh {

class CascadePartitioner;

/*!
 * @brief A grouping of processes in the CascadePartitioner algorithm,
 * which is either a single-process group or a combination of two
 * groups.
 *
 * @b Terminology: The first groups in the CascadePartitioner are
 * single-process groups.  With each CascadePartitioner cycle, two
 * adjacent groups are combined to make a bigger group.  Cycle c has
 * $2^c$ ranks in each group.  The two constituent groups in combined
 * group are sometimes called "halves", because each makes up half of
 * the combined group.  The lower half has smaller ranks than the
 * upper half.  In addition, the half containing the local process is
 * called "our half", while the one not containing the local process
 * is called the "far half".
 *
 * A group can shift load from its overloaded half to its underloaded
 * half, and this is how the CascadePartitioner balance loads.
 *
 */
class CascadePartitionerCombinedGroup {

public:

   CascadePartitionerCombinedGroup() :
      d_common(0),
      d_cycle_num(-1),
      d_lower_begin(-1),
      d_upper_begin(-1),
      d_upper_end(-1),
      d_contact(-1),
      d_our_half(0),
      d_our_position(Lower),
      d_lower_work(0.0),
      d_upper_work(0.0),
      d_lower_capacity(0.0),
      d_upper_capacity(0.0),
      d_our_work(0),
      d_far_work(0),
      d_local_may_supply(true),
      d_contact_may_supply(true),
      d_our_half_may_supply(true),
      d_far_half_may_supply(true) {}

   /*!
    * @brief Copy constructor doesn't copy anything.  It is not used but is required
    * for the class to be used in an stl::vector.
    */
   CascadePartitionerCombinedGroup( const CascadePartitionerCombinedGroup &other ) :
      d_common(0),
      d_cycle_num(-1),
      d_lower_begin(-1),
      d_upper_begin(-1),
      d_upper_end(-1),
      d_contact(-1),
      d_our_half(0),
      d_our_position(Lower),
      d_lower_work(0.0),
      d_upper_work(0.0),
      d_lower_capacity(0.0),
      d_upper_capacity(0.0),
      d_our_work(0),
      d_far_work(0),
      d_local_may_supply(true),
      d_contact_may_supply(true),
      d_our_half_may_supply(true),
      d_far_half_may_supply(true) {}

   ~CascadePartitionerCombinedGroup() {}

   /*!
    * @brief Make a cycle-zero, single-process group.
    */
   void makeSingleProcessGroup(
      const CascadePartitioner *common_data,
      TransitLoad &local_load );

   /*!
    * @brief Make a combined group consisting of the given half-group
    * and the other half-group, which this method will figure out.
    */
   void makeCombinedGroup( CascadePartitionerCombinedGroup &our_half );

   /*!
    * @brief Improve balance of the two halves of this group by
    * supplying load from overloaded half to underloadef half.
    *
    * Ideally, the work supplied is minimum of the overloaded half's
    * surplus and the underloaded half's deficit.  The ideal may not
    * be achieved due to load-cutting restrictions.
    */
   void balanceConstituentHalves();

   void printClassData( std::ostream &co, const std::string &border ) const;


private:

   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int CascadePartitionerCombinedGroup_TAG_InfoExchange = 1000;
   static const int CascadePartitionerCombinedGroup_TAG_LoadTransfer0 = 1001;
   static const int CascadePartitionerCombinedGroup_TAG_LoadTransfer1 = 1002;

   //! @brief Where a group falls in the next larger group.
   enum Position { Lower=0, Upper=1 };

   //! @brief Work of the group (combined work of its two halves).
   double getCombinedWork() const { return d_lower_work + d_upper_work; }

   //! @brief Surplus (estimated) of lower half.
   double lowerSurplus() const {
      return d_lower_work - d_lower_capacity;
   }
   //! @brief Surplus (estimated) of upper half.
   double upperSurplus() const {
      return d_upper_work - d_upper_capacity;
   }
   //! @brief Surplus (estimated) of our half.
   double ourSurplus() const {
      return d_our_position == Lower ? lowerSurplus() : upperSurplus();
   }
   //! @brief Surplus (estimated) of far half.
   double farSurplus() const {
      return d_our_position == Lower ? upperSurplus() : lowerSurplus();
   }

   /*!
    * @brief Try to supply the requested amount of work by removing
    * it from this group, and return the (estimated) amount supplied.
    *
    * Supplying work returns an estimate of the amount supplied, based
    * on the work available and assuming perfect load cutting.  Due to
    * restrictions such as in box cutting, the actual amount supplied
    * may differ.  Using correct values locally and estimates for far
    * groups will lead to discrpepancies in record-keeping.
    *
    * Single-process groups set aside any work it personally gives up
    * in d_common->d_shipment
    *
    * @param taker Representative of the group getting this work.
    *
    * @return Work supplied (or an estimate)
    */
   double supplyWork( double work_requested, int taker );

   /*!
    * @brief Try to supply the requested amount of work from our half
    * of the group.
    *
    * The return value is exact if the group includes only the local
    * rank, otherwise, it's an estimate based on the requested supply.
    *
    * @param taker Rank of process taking load from local process.
    *
    * @return Work supplied (or an estimate)
    */
   double supplyWorkFromOurHalf( double work_requested, int taker );

   /*!
    * @brief Symbolically to supply the requested amount of work from
    * the half of the group not containing the local process.
    *
    * No real work is exchanged because the local process is not in
    * the far half.  This method just estimates what the far half
    * could give away.
    *
    * @return Estimate of work supplied (actual value not locally
    * available)
    */
   double supplyWorkFromFarHalf( double work_requested );

   void sendMyShipment( int taker );
   void receiveAndUnpackSuppliedLoad();

   const CascadePartitioner *d_common;

   //! @brief Cycle number.  Group has 2^d_cycle_num ranks.
   int d_cycle_num;

   //! @brief First rank in lower half.
   int d_lower_begin;

   //! @brief First rank in upper half.
   int d_upper_begin;

   //! @brief One past the last rank.
   int d_upper_end;

   //! Rank of contact (if any) in the far half of the group.
   int d_contact;

   //! @brief The half containing me.
   CascadePartitionerCombinedGroup *d_our_half;

   //! @brief Position of d_our_half in this group.
   Position d_our_position;

   //! @brief Sum of load held by lower half (or approxmimation).
   double d_lower_work;

   //! @brief Sum of load held by upper half (or approxmimation).
   double d_upper_work;

   //! @brief Ideal load based on number of ranks in lower half.
   double d_lower_capacity;

   //! @brief Ideal load based on number of ranks in upper half.
   double d_upper_capacity;

   //! @brief Points to either d_lower_work or d_upper_work.
   double *d_our_work;

   //! @brief Points to either d_lower_work or d_upper_work.
   double *d_far_work;

   //@{
   //! @name For determining whether to send/receive, not for computing how much load to transfer.

   //! @brief Whether local process may supply load.
   bool d_local_may_supply;

   //! Whether contact may supply load.
   bool d_contact_may_supply;

   //! @brief Whether our half may supply load.
   bool d_our_half_may_supply;

   //! @brief Whether far half may supply load.
   bool d_far_half_may_supply;

   //@}

};

}
}

#endif
