/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Fast assumed partition for a box.
 *
 ************************************************************************/
#include "SAMRAI/hier/AssumedPartitionBox.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {


AssumedPartitionBox::AssumedPartitionBox(
   const Box& box,
   int rank_begin,
   int rank_end,
   int index_begin):
   d_box(box),
   d_rank_begin(rank_begin),
   d_rank_end(rank_end),
   d_uniform_partition_size(box.getDim()),
   d_partition_grid_size(box.getDim()),
   d_major(box.getDim()),
   d_index_stride(box.getDim()),
   d_index_begin(index_begin)
{
   computeLayout();
   assignToRanks();
   return;
}



/*
***************************************************************************************
* Return index of first box assigned to given rank.
***************************************************************************************
*/
int
AssumedPartitionBox::beginOfRank(int rank) const
{
   if ( rank >= d_rank_begin && rank < d_rank_end ) {
      int index = rank < d_first_rank_with_1 ? d_first_index_with_2 + (rank-d_first_rank_with_1)*2 : 0
                + rank < d_first_rank_with_0 ? d_first_index_with_1 + (rank-d_first_rank_with_0)   : 0;
      return index;
   }
   return d_index_begin + d_partition_grid_size.getProduct();
}



/*
***************************************************************************************
* Return one past index of last box assigned to given rank.
***************************************************************************************
*/
int
AssumedPartitionBox::endOfRank(int rank) const
{
   return beginOfRank(rank) + rank < d_first_rank_with_1 ? 2 : rank < d_first_rank_with_0 ? 1 : 0;
}



/*
***************************************************************************************
* Compute the owner of the given index.
***************************************************************************************
*/
int
AssumedPartitionBox::getOwner(int box_index) const
{
   TBOX_ASSERT( box_index >= d_index_begin );
   TBOX_ASSERT( box_index < d_index_begin + getNumberOfParts() );

   int owner = tbox::SAMRAI_MPI::getInvalidRank();
   if ( box_index < d_index_begin || box_index >= d_first_index_with_0 ) {
      // Not an index in this object.  Return invalid owner.
   }
   if ( box_index < d_first_index_with_1 ) {
      owner = d_rank_begin + (box_index - d_first_index_with_2)/2;
   }
   else if ( box_index < d_first_index_with_0 ) {
      owner = d_first_rank_with_1 + (box_index - d_first_index_with_1);
   }
   return owner;
}



/*
***************************************************************************************
* Compute the box with the given index.
***************************************************************************************
*/
Box
AssumedPartitionBox::getBox(int box_index) const
{
   TBOX_ASSERT( box_index >= d_index_begin );
   TBOX_ASSERT( box_index < d_index_begin + d_partition_grid_size.getProduct() );

   Box part(d_box);
   int box_index_diff = box_index - d_index_begin;
   for ( int d=d_box.getDim().getValue()-1; d>=0; --d ) {
      int dir = d_major[d];
      part.lower()[dir] = box_index_diff/d_index_stride[dir];
      box_index_diff -= part.lower()[dir]*d_index_stride[dir];
   }
   part.lower() *= d_uniform_partition_size;
   part.lower() += d_box.lower();
   part.upper() = part.lower() + d_uniform_partition_size - IntVector::getOne(d_box.getDim());
   part *= d_box;

   const int owner = getOwner(box_index);
   part.initialize( part, LocalId(box_index), owner );

   return part;
}



/*
***************************************************************************************
* Compute the box with the given position in the grid of partitions.
***************************************************************************************
*/
Box
AssumedPartitionBox::getBox(const IntVector &position) const
{
   TBOX_ASSERT( position >= hier::IntVector::getZero(d_box.getDim()) );
   TBOX_ASSERT( position < d_partition_grid_size );

   int box_index = d_index_begin;
   for ( int d=0; d<d_box.getDim().getValue(); ++d ) {
      box_index += position[d]*d_index_stride[d];
   }
   const int owner = getOwner(box_index);
   Box box( Index(position),
            Index(position + d_uniform_partition_size - IntVector::getOne(d_box.getDim())),
            d_box.getBlockId(),
            LocalId(box_index),
            owner );
   box.coarsen(d_uniform_partition_size);
   box *= d_box;
   return box;
}



/*
***************************************************************************************
* Find all boxes intersecting the given box.  Return whether any boxes overlap.
***************************************************************************************
*/
bool
AssumedPartitionBox::findOverlaps(
   BoxContainer &overlapping_boxes,
   const Box &box ) const
{
   Box coarsened_box = box;
   coarsened_box.coarsen( d_uniform_partition_size );
   coarsened_box *= Box( Index(IntVector::getZero(d_box.getDim())),
                         Index(d_partition_grid_size - IntVector::getOne(d_box.getDim())),
                         d_box.getBlockId() );
   for ( Box::iterator ci=coarsened_box.begin(); ci!=coarsened_box.end(); ++ci ) {
      overlapping_boxes.insert( getBox(*ci) );
   }
   return !coarsened_box.empty();
}



/*
***************************************************************************************
* Check the assumed partition for errors and inconsistencies.  Write
* error diagnostics to plog.
*
* Failure indicates a bug in this class.
***************************************************************************************
*/
size_t
AssumedPartitionBox::selfCheck() const
{
   size_t nerr = 0;

   BoxContainer tmp_boxes;

   BoxContainer all_parts;
   for ( int box_index=begin(); box_index!=end(); ++box_index ) {
      all_parts.pushBack( getBox(box_index) );
   }
   all_parts.makeTree();


   // Make sure no part intersects others.
   tmp_boxes.clear();
   for ( BoxContainer::const_iterator bi=all_parts.begin(); bi!=all_parts.end(); ++bi ) {
      const Box &box = *bi;
      tmp_boxes.clear();
      all_parts.findOverlapBoxes( tmp_boxes, box );
      tmp_boxes.order();
      if ( !tmp_boxes.empty() ) {
         BoxContainer::iterator self = tmp_boxes.find(box);
         if ( self != tmp_boxes.end() ) {
            tmp_boxes.erase(self);
         }
      }
      if ( !tmp_boxes.empty() ) {
         nerr += tmp_boxes.size();
         tbox::plog << "AssumedPartitionerBox::selfCheck(): Box "
                    << box << " unexpectedly overlaps these:\n"
                    << tmp_boxes.format("\t") << std::endl;
      }
   }

   // Part should cover all of d_box.
   BoxContainer box_leftover(d_box);
   box_leftover.removeIntersections(all_parts);
   if ( !box_leftover.empty() ) {
      nerr += box_leftover.size();
      tbox::plog << "AssumedPartitionerBox::selfCheck(): Partitions fail to cover box "
                 << d_box << "  Uncovered parts:\n"
                 << box_leftover.format("\t") << std::endl;
   }

   // Parts shoud cover no more than d_box.
   BoxContainer parts_leftover = all_parts;
   parts_leftover.removeIntersections(d_box);
   if ( !parts_leftover.empty() ) {
      nerr += parts_leftover.size();
      tbox::plog << "AssumedPartitionerBox::selfCheck(): Partitions overflow  box "
                 << d_box << "  Overflow parts:\n"
                 << parts_leftover.format("\t") << std::endl;
   }

   return nerr;
}



/*
***************************************************************************************
* Compute the partition lay-out.  We use a grid of uniform sized
* partitions whose union covers d_box and as little else as possible.
*
* TODO: experiment with other layouts that minimize overflowing d_box
* and have minimum aspect ratio in the partition size.
***************************************************************************************
*/
void
AssumedPartitionBox::computeLayout()
{
   const IntVector box_size = d_box.numberCells();

   const int num_ranks = d_rank_end - d_rank_begin;

   /*
    * Compute uniform partition size and how many partitions in each
    * direction.  There isn't one correct lay-out, but we try to avoid
    * excessive aspect ratios.
    */
   d_uniform_partition_size = box_size;
   d_partition_grid_size = hier::IntVector::getOne(d_box.getDim());
   int parts_count = d_partition_grid_size.getProduct();
   IntVector num_parts_can_increase(d_box.getDim(), 1);
   IntVector sorter(d_box.getDim());
   while ( parts_count < num_ranks &&
           num_parts_can_increase != hier::IntVector::getZero(d_box.getDim()) ) {
      sorter.sortIntVector(d_uniform_partition_size);
      int inc_dir = 0;
      for ( inc_dir=d_box.getDim().getValue()-1; inc_dir>=0; --inc_dir ) {
         if ( num_parts_can_increase[sorter[inc_dir]] ) break;
      }
      inc_dir = sorter[inc_dir];

      if ( 2*parts_count > num_ranks ) {
         const int cross_section = parts_count/d_partition_grid_size[inc_dir];
         d_partition_grid_size[inc_dir] = (num_ranks+cross_section-1)/cross_section;
         parts_count = d_partition_grid_size.getProduct();
      } else {
         d_partition_grid_size[inc_dir] *= 2;
         parts_count *= 2;
      }

      d_uniform_partition_size = IntVector::ceilingDivide(box_size, d_partition_grid_size);
      num_parts_can_increase[inc_dir] = d_uniform_partition_size[inc_dir] > 1;
   }
   TBOX_ASSERT( parts_count == d_partition_grid_size.getProduct() );
   TBOX_ASSERT( d_uniform_partition_size.getProduct() > 0 );

   // There can be partitions completele outside d_box.  Remove them.
   d_partition_grid_size = IntVector::ceilingDivide( box_size, d_uniform_partition_size );

   return;
}



/*
***************************************************************************************
* Compute rank assignment for the partition lay-out.
***************************************************************************************
*/
void
AssumedPartitionBox::assignToRanks()
{
   // There should not be more than 2 partitions per rank.
   TBOX_ASSERT( d_partition_grid_size.getProduct() <= 2*(d_rank_end-d_rank_begin) );

   /*
    * If there are more ranks than parts, the first getNumberOfParts()
    * ranks have 1 part each and the rest have none.  If there are
    * more parts than ranks, lower ranks have 2 parts and higher ranks
    * have 1.  In second case, index vs rank looks like this:
    *
    * index ^
    *       |
    *    i0 |             ......
    *       |          .
    *       |       .
    *    i1 |    .
    *       |   .
    *       |  .
    *       | .
    *       |.
    *    i2 +-----------------------> rank
    *       r2   r1       r0
    *
    * (In the first case, the figure degenerates with r2=r1 and i2=i1.)
    */
   if ( d_partition_grid_size.getProduct() <= d_rank_end-d_rank_begin ) {
      d_first_rank_with_1 = d_rank_begin;
      d_first_rank_with_0 = d_rank_begin + d_partition_grid_size.getProduct();
      d_first_index_with_1 = d_first_index_with_2 = d_index_begin;
   }
   else {
      d_first_index_with_2 = d_index_begin;
      d_first_rank_with_1 = d_partition_grid_size.getProduct() % (d_rank_end-d_rank_begin);
      d_first_index_with_1 = 2*(d_first_rank_with_1 - d_rank_begin);
      d_first_rank_with_0 = d_rank_end;
   }
   d_first_index_with_0 = d_index_begin + getNumberOfParts();

   d_major.sortIntVector(d_partition_grid_size);

   for ( int d=0; d<d_box.getDim().getValue(); ++d ) {
      int dir = d_major[d];
      d_index_stride[dir] = 1;
      for ( int d1=d-1; d1>=0; --d1 ) {
         d_index_stride[dir] *= d_partition_grid_size[d_major[d1]];
      }
   }

   return;
}



/*
***************************************************************************************
***************************************************************************************
*/
void
AssumedPartitionBox::recursivePrint(
   std::ostream &co,
   const std::string& border,
   int detail_depth) const
{
   const char *to = "..";
   co << border << "d_box = " << d_box
      << '\n' << border
      << "d_uniform_partition_size = " << d_uniform_partition_size
      << "  d_partition_grid_size = " << d_partition_grid_size
      << '\n' << border
      << "indices (" << d_partition_grid_size.getProduct() << "): "
      << d_index_begin << to << d_index_begin+d_partition_grid_size.getProduct()-1
      << "    ranks (" << d_rank_end-d_rank_begin << "): "
      << d_rank_begin << to << d_rank_end-1
      << '\n' << border
      << "ranks with 2 parts (" << d_first_rank_with_1-d_rank_begin << "): "
      << d_rank_begin << to << d_first_rank_with_1-1
      << "    ranks with 1 parts (" << d_first_rank_with_0-d_first_rank_with_1 << "): "
      << d_first_rank_with_1 << to << d_first_rank_with_0-1
      << "    ranks with 0 parts (" << d_rank_end-d_first_rank_with_0 << "): "
      << d_first_rank_with_0 << to << d_rank_end-1
      << std::endl;
   if ( detail_depth > 0 ) {
      BoxContainer parts;
      for ( int box_index=begin(); box_index!=end(); ++box_index ) {
         parts.pushBack( getBox(box_index) );
      }
      co << border << "Parts:" << parts.format(border);
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
