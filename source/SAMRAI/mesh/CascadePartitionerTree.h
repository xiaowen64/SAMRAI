/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitionerTree
#define included_mesh_CascadePartitionerTree

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
 * @brief A binary-tree of process groups in the CascadePartitioner
 * algorithm,
 *
 * In this partitioner, the MPI ranks are recursively split into
 * groups, forming the nodes of a binary tree.  The root branch
 * contains all ranks.  The leaves are single-process groups.  Each
 * group is represented by a CascadePartitionerTree.
 * repre
 *
 * We balance loads one sibling pair at a time, starting with the
 * leaves, shifting load from overloaded groups to underloaded groups.
 * As we move toward the groups, the loads can be propagated into ever
 * bigger groups.
 *
 * @b Terminology: The root group includes all ranks.  It splits into
 * lower and upper groups (known as branches).  The lower branch has
 * the lower ranks.  If the group has an odd number of ranks, the
 * upper branch has one more rank in it.  The branch containing the
 * local process is also called "near branch", while the one not
 * containing the local process is the "far branch".
 *
 * A parent can shift load from its overloaded child to its
 * underloaded child, and this is how the CascadePartitioner balance
 * loads.
 *
 */
class CascadePartitionerTree {

public:

   /*!
    * @brief Construct the tree for the given CascadePartitioner by
    * constructing the root and recursively constructing its children.
    */
   CascadePartitionerTree( const CascadePartitioner &partitioner );

   ~CascadePartitionerTree();

   /*!
    * @brief Run the complete cascade partitioner algorithm.
    */
   void balanceAll();

   void printClassData( std::ostream &co, const std::string &border ) const;


   //! @brief Generation number (generation 0 contains all ranks).
   int generationNum() const { return d_gen_num; }

   //! @brief Cycle number (cycle i groups has 2^i processes).
   int cycleNum() const;

   //! @brief Size of group (number of processes in it).
   size_t size() const {
      return d_end-d_begin;
   }

   //! @brief Whether group contains a given rank.
   bool containsRank( int rank ) const {
      return d_begin <= rank && rank < d_end;
   }


private:

   //! @brief Where a group falls in the next larger group.
   enum Position { Lower=0, Upper=1 };

   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int CascadePartitionerTree_TAG_InfoExchange0 = 1000;
   static const int CascadePartitionerTree_TAG_InfoExchange1 = 1001;
   static const int CascadePartitionerTree_TAG_LoadTransfer0 = 1002;
   static const int CascadePartitionerTree_TAG_LoadTransfer1 = 1003;

   //! @brief Where a group falls in the next larger group.
   static const int s_lower = Lower;
   static const int s_upper = Upper;

   /*!
    * @brief Construct child node based on its position in the parent.
    */
   CascadePartitionerTree( CascadePartitionerTree &parent,
                           Position group_position );

   //! Allocate and set up the group's children.
   void makeChildren();

   /*!
    * @brief Combine near and far children data (using communication)
    * to compute work-related data for this group.
    */
   void combineChildren();

   /*!
    * @brief Improve balance of the two children of this group by
    * supplying load from overloaded child to underloaded child.
    *
    * Ideally, the work supplied is minimum of the overloaded child's
    * surplus and the underloaded child's deficit.  The ideal may not
    * be achieved due to load-cutting restrictions.
    */
   void balanceChildren();

   //! @brief Return sibling group.
   CascadePartitionerTree *sibling() {
      return d_parent->d_children[!d_position];
   }

   //! @brief Estimated surplus of the group.
   double estimatedSurplus() const {
      return d_work - d_obligation;
   }

   /*!
    * @brief Try to supply the requested amount of work by removing
    * it from this group, and return the (estimated) amount supplied.
    *
    * Supplying work returns an estimate of the amount supplied, based
    * on available surplus and assuming perfect load cutting.  Due to
    * restrictions such as in box cutting, the actual amount supplied
    * may differ.  Actual amount is available when the group contains
    * just the local process, but the estimate is always available.
    * We always use estimates to discrpepancies in record-keeping.
    *
    * Single-process groups set aside any work it personally gives up
    * in d_common->d_shipment
    *
    * @param taker Representative of the group getting this work.
    *
    * @return Work supplied estimate based on perfect load cutting
    */
   double supplyWork( double work_requested, int taker );

   void sendShipment( int taker );
   void receiveAndUnpackSuppliedLoad();

   //! @brief Recompute data for leaf groups.
   void recomputeLeafData();

   /*!
    * @brief Reset obligation recursively for all descendents.
    *
    * @param avg_load Average per-process load in the group.
    */
   void resetObligation( double avg_load );


   //! @brief Data the main CascadePartitioner shares with all parts of the tree.
   const CascadePartitioner *d_common;

   //@{
   //! @brief Group specification

   //! @brief Generation number.  Generation 0 contains all ranks.
   int d_gen_num;

   //! @brief First rank in group.
   int d_begin;

   //! @brief One past last rank in group.
   int d_end;

   //! @brief Position of this group in the parent.
   int d_position;

   /*!
    * @brief Rank of contacts in sibling branch.
    *
    * Communication between sibling groups is done between a process
    * in one group and its contact in the sibling group.  Contacting
    * pairs are assigned based on the process's relative rank in its
    * group.
    *
    * Most processes have just one contact in the sibling group.  The
    * only processes to have 2 contacts are in the lower sibling of a
    * parent that has an odd number of ranks.  In these cases, the
    * upper sibling has one more rank than the lower.  The last rank
    * of the lower sibling contacts the last 2 ranks in the upper.
    */
   int d_contact[2];

   //! @brief Parent group.
   CascadePartitionerTree *d_parent;

   //! @brief Lower and upper children branches.
   CascadePartitionerTree *d_children[2];

   //! @brief Near child branch (branch containing local process).
   CascadePartitionerTree *d_near;

   //! @brief Far child branch (branch not containing local process).
   CascadePartitionerTree *d_far;

   //! @brief Group containing just the local process.
   CascadePartitionerTree *d_leaf;

   //@}



   //@{
   //! @name Work measures

   //! @brief Estimated load of this branch.
   double d_work;

   //! @brief Amount of work the group obligated to have.
   double d_obligation;

   //@}



   //@{
   //! @name For determining participation in certain group activities.

   /*!
    * @brief Whether this group may supply work to its sibling.
    *
    * A group that received work may not later become a supplier.
    * (But its parent may still be if its sibling is.)
    * Note: This should permission should be propagated to ancestor
    * groups but not descendent groups, because work received should
    * be further distributed to descendents.  Work received should
    * not be supplied to sibling, because that allows load to move
    * back and forth for trivial imbalances.
    */
   bool d_group_may_supply;

   /*!
    * @brief Whether local process (if this is a near group) or
    * contact processes (if this is a far group) may supply load.
    *
    * If this is a near group, first value specifies whether the local
    * process may supply work and second value is unused.  If this is
    * a far group, the two values correspond to whether the two
    * contacts (d_contact) may supply work.
    *
    * A process that received work or is a member of a group that
    * received work may not later become a supplier.
    */
   bool d_process_may_supply[2];

   //@}

};

}
}

#endif
