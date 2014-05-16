/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Fast assumed partition for a set of boxes.
 *
 ************************************************************************/
#ifndef included_hier_AssumedPartition
#define included_hier_AssumedPartition

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/AssumedPartitionBox.h"
#include "SAMRAI/hier/BoxContainer.h"


namespace SAMRAI {
namespace hier {

/*!
 * @brief Compute an assumed partition of a set of boxes.  The assumed
 * partition should very fast to create and query and requires minimal
 * storage.
 *
 * An assumed partition should avoid extreme imbalances, but its
 * purpose is not fine load balancing.
 *
 * See also AssumedPartitionBox.
 */
class AssumedPartition
{

public:

   /*!
    * @brief Construct AssumedPartition of a set of boxes.
    *
    * @param[in] boxes Incoming boxes
    *
    * @param[in] rank_begin First rank
    *
    * @param[in] rank_end One past last rank
    *
    * @param[in] index_begin
    *
    * @param[in] parts_per_rank See partition()
    *
    * @param[in] interleave See partition()
    */
   AssumedPartition(
      const BoxContainer& boxes,
      int rank_begin,
      int rank_end,
      int index_begin = 0,
      double parts_per_rank = 1.0,
      bool interleave = false );

   /*!
    * @brief Construct an empty AssumedPartition.
    */
   AssumedPartition();

   /*!
    * @brief Destructor.
    */
   ~AssumedPartition() {}

   /*!
    * @brief Partition a set of boxes.
    *
    * @param[in] boxes Incoming boxes
    *
    * @param[in] rank_begin First rank
    *
    * @param[in] rank_end One past last rank
    *
    * @param[in] parts_per_rank Algorithm normally tries to get one partition
    *   per rank.  This parameter is a request to change that.
    *
    * @param[in] interleave Algorithm normally assign consecutive box indices
    *   to a process.  This flag causes it to interleave (round-robin) the
    *   box assignments.
    */
   void partition(
      const BoxContainer& boxes,
      int rank_begin,
      int rank_end,
      int index_begin = 0,
      double parts_per_rank = 1.0,
      bool interleave = false );

   //! @brief Number of box partitions.
   size_t getNumberOfParts() const {
      return d_index_end - d_index_begin;
   }

   //! @brief Return the owner for a box.
   int getOwner(int box_index) const;

   //! @brief Return box for given index.
   Box getBox(int box_index) const;

   //! @brief Get all boxes.
   void getAllBoxes( BoxContainer &all_boxes ) const;

   //! @brief Get all boxes for a given rank.
   void getAllBoxes( BoxContainer &all_boxes, int rank ) const;

   //! @brief Return index of first box.
   int begin() const {
      return d_index_begin;
   }

   //! @brief Return one past index of last box.
   int end() const {
      return d_index_end;
   }

   //! @brief Return index of first box assigned to given rank.
   int beginOfRank(int rank) const;

   //! @brief Return one past index of last box assigned to given rank.
   int endOfRank(int rank) const;

   /*!
    * @brief Find partitions overlapping the given box, when there is
    * at most one original box, assuming all boxes are in the same
    * block.
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
    * @brief Find partitions overlapping the given box.
    *
    * @param[o] overlapping_boxes
    * @param[i] box
    * @param[i] grid_geometry
    * @param[i] refinement_ratio
    *
    * @return Whether any partitions are found.
    */
   bool findOverlaps(
      BoxContainer &overlapping_boxes,
      const Box &box,
      const BaseGridGeometry &grid_geometry,
      const IntVector &refinement_ratio ) const;

   /*!
    * @brief Check the assumed partition for errors and
    * inconsistencies.  Write error diagnostics to plog.
    *
    * Failure indicates a bug in this class.
    */
   size_t selfCheck() const;

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

   //! @brief Shorthand.
   typedef std::vector<AssumedPartitionBox> PartedBoxes;

   PartedBoxes d_parted_boxes;

   int d_rank_begin;
   int d_rank_end;
   int d_index_begin;
   int d_index_end;

};

}
}

#endif  // included_hier_AssumedPartition
