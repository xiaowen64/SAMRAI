/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Fast assumed partition for a set of boxes.
 *
 ************************************************************************/
#include "SAMRAI/hier/AssumedPartition.h"
#include "SAMRAI/hier/BaseGridGeometry.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {


/*
***************************************************************************************
***************************************************************************************
*/
AssumedPartition::AssumedPartition() :
   d_parted_boxes(),
   d_rank_begin(0),
   d_rank_end(0),
   d_index_begin(0),
   d_index_end(0)
{
}


/*
***************************************************************************************
***************************************************************************************
*/
AssumedPartition::AssumedPartition(
   const BoxContainer& boxes,
   int rank_begin,
   int rank_end,
   int index_begin,
   double parts_per_rank,
   bool interleave) :
   d_parted_boxes(),
   d_rank_begin(rank_begin),
   d_rank_end(rank_end),
   d_index_begin(index_begin),
   d_index_end(index_begin)
{
   partition(boxes, rank_begin, rank_end, index_begin, parts_per_rank, interleave);
}


/*
***************************************************************************************
* Partition the given boxes.
***************************************************************************************
*/
void
AssumedPartition::partition(
   const BoxContainer& boxes,
   int rank_begin,
   int rank_end,
   int index_begin,
   double parts_per_rank,
   bool interleave)
{
   TBOX_ASSERT( rank_end > rank_begin );
   d_parted_boxes.clear();
   d_rank_begin = rank_begin;
   d_rank_end = rank_end;
   d_index_begin = index_begin;
   d_index_end = index_begin;

   size_t num_cells = 0;
   for ( BoxContainer::const_iterator bi=boxes.begin(); bi!=boxes.end(); ++bi ) {
      num_cells += bi->size();
   }

   d_parted_boxes.reserve(boxes.size());

   const int num_ranks = rank_end - rank_begin;
   double rank_space_cut_lo = 0.0;
   double rank_space_cut_hi = 0.0;
   int box_rank_begin = rank_begin;
   int box_rank_end = box_rank_begin;
   size_t cell_count = 0;
   for ( BoxContainer::const_iterator bi=boxes.begin(); bi!=boxes.end(); ++bi ) {
      cell_count += bi->size();

      rank_space_cut_lo = rank_space_cut_hi;
      rank_space_cut_hi = static_cast<double>(cell_count)/num_cells;

      int box_rank_begin = static_cast<int>(rank_space_cut_lo * num_ranks + 0.5);
      int box_rank_end   = static_cast<int>(rank_space_cut_hi * num_ranks + 0.5);
      if ( box_rank_end == box_rank_begin ) {
         const double midpoint = 0.5*(rank_space_cut_lo + rank_space_cut_hi);
         box_rank_begin = static_cast<int>(midpoint * num_ranks + 0.5);
         box_rank_begin = tbox::MathUtilities<int>::Min( box_rank_begin, d_rank_end-1 );
         box_rank_end = box_rank_begin + 1;
      }

      d_parted_boxes.push_back(
         AssumedPartitionBox(*bi, box_rank_begin, box_rank_end, d_index_end,
                             parts_per_rank, interleave) );
      d_index_end = d_parted_boxes.back().end();
   }

}



/*
***************************************************************************************
* Return index of first box assigned to given rank.
***************************************************************************************
*/
int
AssumedPartition::beginOfRank(int rank) const
{
   for ( size_t i=0; i<d_parted_boxes.size(); ++i ) {
      const int begin = d_parted_boxes[i].beginOfRank(rank);
      const int end = d_parted_boxes[i].endOfRank(rank);
      if ( end > begin ) {
         return begin;
      }
   }
   return d_index_end;
}



/*
***************************************************************************************
* Return one past index of last box assigned to given rank.
***************************************************************************************
*/
int
AssumedPartition::endOfRank(int rank) const
{
   for ( PartedBoxes::const_reverse_iterator pi=d_parted_boxes.rbegin(); pi!=d_parted_boxes.rend(); ++pi ) {
      const int begin = pi->beginOfRank(rank);
      const int end = pi->endOfRank(rank);
      if ( end > begin ) {
         return end;
      }
   }
   return d_index_end;
}



/*
***************************************************************************************
* Compute the box with the given index.
***************************************************************************************
*/
int
AssumedPartition::getOwner(int box_index) const
{
   TBOX_ASSERT( box_index >= d_index_begin );
   TBOX_ASSERT( box_index < d_index_end );

   for ( size_t i=0; i<d_parted_boxes.size(); ++i ) {
      if ( box_index >= d_parted_boxes[i].begin() &&
           box_index < d_parted_boxes[i].end() ) {
         return d_parted_boxes[i].getOwner(box_index);
      }
   }

   TBOX_ERROR("AssumedPartition::getBox(): Should never be here.");
   return tbox::SAMRAI_MPI::getInvalidRank();
}



/*
***************************************************************************************
* Compute the box with the given index.
***************************************************************************************
*/
Box
AssumedPartition::getBox(int box_index) const
{
   TBOX_ASSERT( box_index >= d_index_begin );
   TBOX_ASSERT( box_index < d_index_end );

   for ( size_t i=0; i<d_parted_boxes.size(); ++i ) {
      if ( box_index >= d_parted_boxes[i].begin() &&
           box_index < d_parted_boxes[i].end() ) {
         return d_parted_boxes[i].getBox(box_index);
      }
   }

   TBOX_ERROR("AssumedPartition::getBox(): Should never be here.");
   return Box(tbox::Dimension(1));
}



/*
***************************************************************************************
* Compute the box with the given index.
***************************************************************************************
*/
void
AssumedPartition::getAllBoxes(BoxContainer &all_boxes) const
{
   for ( int id=d_index_begin; id<d_index_end; ++id ) {
      d_parted_boxes[id].getAllBoxes(all_boxes);
   }
   return;
}



/*
***************************************************************************************
* Compute the box with the given index.
***************************************************************************************
*/
void
AssumedPartition::getAllBoxes(BoxContainer &all_boxes, int rank) const
{
   const int id_begin = beginOfRank(rank);
   const int id_end = endOfRank(rank);
   for ( int id=id_begin; id<id_end; ++id ) {
      all_boxes.push_back(getBox(id));
   }
   return;
}



/*
***************************************************************************************
* Find all boxes intersecting the given box.  Return whether any boxes overlap.
***************************************************************************************
*/
bool
AssumedPartition::findOverlaps(
   BoxContainer &overlapping_boxes,
   const Box &box ) const
{
   TBOX_ASSERT( d_parted_boxes.size() < 2 );
   if ( !d_parted_boxes.empty() ) {
      return d_parted_boxes[0].findOverlaps(overlapping_boxes, box);
   }
   return false;
}



/*
***************************************************************************************
* Find all boxes intersecting the given box.  Return whether any boxes overlap.
***************************************************************************************
*/
bool
AssumedPartition::findOverlaps(
   BoxContainer &overlapping_boxes,
   const Box &box,
   const BaseGridGeometry &grid_geometry,
   const IntVector &refinement_ratio ) const
{
   size_t old_count = overlapping_boxes.size();

   for ( PartedBoxes::const_iterator pi=d_parted_boxes.begin(); pi!=d_parted_boxes.end(); ++pi ) {
      Box transformed_box = box;
      grid_geometry.transformBox(transformed_box,
                                 refinement_ratio,
                                 pi->getUnpartitionedBox().getBlockId(),
                                 box.getBlockId());
      pi->findOverlaps(overlapping_boxes, transformed_box);
   }

   return overlapping_boxes.size() > old_count;
}



/*
***************************************************************************************
* Check the assumed partition for errors and inconsistencies.  Write
* error diagnostics to plog.
*
* Return number of errors found.  Errors indicate a bug in this class.
***************************************************************************************
*/
size_t
AssumedPartition::selfCheck() const
{
   size_t nerr = 0;

   for ( size_t i=0; i<d_parted_boxes.size(); ++i ) {
      nerr += d_parted_boxes[i].selfCheck();
      if ( i > 0 ) {
         if ( d_parted_boxes[i].begin() < d_parted_boxes[i-1].end() ) {
            tbox::plog << "AssumedPartition::selfCheck(): index from unpartitioned box "
                       << i << " overlaps index from unpartitioned box " << i-1 << std::endl;
            ++nerr;
         }
      }
   }

   return nerr;
}



/*
***************************************************************************************
***************************************************************************************
*/
void
AssumedPartition::recursivePrint(
   std::ostream &co,
   const std::string& border,
   int detail_depth) const
{
   const char *to = "..";
   co << border
      << "indices (" << d_index_end-d_index_begin << "): "
      << d_index_begin << to << d_index_end-1
      << "    ranks (" << d_rank_end-d_rank_begin << "): "
      << d_rank_begin << to << d_rank_end-1
      << '\n';
   co << border << d_parted_boxes.size() << " pre-partition boxes:\n";
   for ( size_t i=0; i<d_parted_boxes.size(); ++i ) {
      co << border << "  pre-partition box " << i << ":\n";
      d_parted_boxes[i].recursivePrint(co, border+"    ", detail_depth-1);
   }
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
