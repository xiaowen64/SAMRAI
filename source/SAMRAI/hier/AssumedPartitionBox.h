/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Fast assumed partition for a single box.
 *
 ************************************************************************/
#ifndef included_hier_AssumedPartitionBox
#define included_hier_AssumedPartitionBox

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/IntVector.h"


namespace SAMRAI {
namespace hier {

/*!
 * @brief Compute an assumed partition of a box.  The assumed
 * partition should very fast to create and query and requires minimal
 * storage.
 *
 * An assumed partition should avoid extreme imalances, but it is not
 * meant for partitioning a PatchLevel.
 */
class AssumedPartitionBox
{

public:

   /*!
    * @brief Construct AssumedPartition from a box.
    *
    * @param[in] in_box Incoming box
    *
    * @param[in] rank_begin First rank
    *
    * @param[in] rank_end One past last rank
    *
    * @param[in] index_begin
    */
   AssumedPartitionBox(
      const Box& box,
      int rank_begin,
      int rank_end,
      int index_begin = 0 );

   /*!
    * @brief Destructor.
    */
   ~AssumedPartitionBox() {}

   //! @brief Number of box partitions.
   size_t getNumberOfParts() const {
      return d_partition_grid_size.getProduct();
   }

   //! @brief Return the owner for a box.
   int getOwner(int box_index) const;

   //! @brief Return box for given index.
   Box getBox(int box_index) const;

   //! @brief Return box for given partition's position in the partition grid.
   Box getBox(const IntVector &position) const;

   //! @brief Return index of first box.
   int begin() const {
      return d_first_index_with_2;
   }

   //! @brief Return one past index of last box.
   int end() const {
      return d_first_index_with_0;
   }

   //! @brief Return index of first box assigned to given rank.
   int beginOfRank(int rank) const;

   //! @brief Return one past index of last box assigned to given rank.
   int endOfRank(int rank) const;

   /*!
    * @brief Check the assumed partition for errors and
    * inconsistencies.  Write error diagnostics to plog.
    *
    * Failure indicates a bug in this class.
    */
   size_t selfCheck() const;

   /*!
    * @brief Find partitions overlapping the given box.
    *
    * @param[o] overlapping_boxes
    * @param[i] box
    *
    * @return Whether any partitions are found.
    */
   bool findOverlaps(
      BoxContainer &overlapping_boxes,
      const Box &box ) const;

   /*!
    * @brief Print info from this object
    *
    * @param[in,out] os The output stream
    * @param[in] border
    * @param[in] detail_depth
    */
   void
   recursivePrint(
      std::ostream& os,
      const std::string& border,
      int detail_depth = 2) const;


private:

   //! @brief Compute the partition lay-out.
   void computeLayout();

   //! @brief Box partitioned.
   Box d_box;
   //! @brief First rank.
   int d_rank_begin;
   //! @brief One past last rank.
   int d_rank_end;
   //! @brief Index for first box.
   int d_index_begin;

   //! @brief Size of each uniform partition.
   IntVector d_uniform_partition_size;
   //! @brief Number of partitions in each direction (size of partition grid).
   IntVector d_partition_grid_size;

   //! @brief Directions sorted from small to big, in d_partition_grid_size.
   IntVector d_major;
   //! @brief Box index stride in each direction.
   IntVector d_index_stride;

   /*
    * Each rank has 0, 1 or 2 partitions.  Lower ranks have more than higher ranks.
    * Ranks in [d_rank_begin,d_first_rank_with_1) have 2 partitions each.
    * Ranks in [d_first_rank_with_1,d_first_rank_with_0) have 1 partition each.
    * Ranks in [d_first_rank_with_0,d_rank_end) have 0 partition.
    */
   //! @brief First rank with 1 partition each.  Ranks < d_first_rank_with_1 has 2 partitions each.
   int d_first_rank_with_1;

   //! @brief First rank with 0 partition.  Ranks >= d_first_rank_with_0 has 0 partition.
   int d_first_rank_with_0;

   //! @brief First index for ranks that own 2 partitions each
   int d_first_index_with_2;

   //! @brief First index for ranks that own 1 partition each
   int d_first_index_with_1;

   //! @brief First index for ranks that own 0 partition each
   int d_first_index_with_0;

};

}
}

#endif  // included_hier_AssumedPartitionBox
