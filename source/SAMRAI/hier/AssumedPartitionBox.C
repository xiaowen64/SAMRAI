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
   int first_index):
   d_box(box),
   d_rank_begin(rank_begin),
   d_rank_end(rank_end),
   d_uniform_partition_size(box.getDim()),
   d_partition_grid_size(box.getDim()),
   d_major(box.getDim()),
   d_index_stride(box.getDim()),
   d_first_index(first_index)
{
   computeLayout();

   // There should not be more than 2 partitions per rank.
   TBOX_ASSERT( d_partition_grid_size.getProduct() <= 2*(d_rank_end-d_rank_begin) );

   /*
    * If there are more ranks than parts, the first getNumberOfParts()
    * ranks have 1 part each and the rest have none.  If there are
    * more parts than ranks, lower ranks have 2 parts and higher ranks
    * have 1.  In second case, index vs rank looks like this:
    *
    *       ^
    * index |
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
      d_first_rank_with_0 = d_rank_end;
      d_first_index_with_1 = d_first_index_with_2 = d_first_index;
   }
   else {
      d_first_index_with_2 = d_first_index;
      d_first_rank_with_1 = d_partition_grid_size.getProduct() % (d_rank_end-d_rank_begin);
      d_first_index_with_1 = 2*(d_first_rank_with_1 - d_rank_begin);
      d_first_rank_with_0 = d_rank_end;
   }
   d_first_index_with_0 = d_first_index + getNumberOfParts();

   d_major.sortIntVector(d_partition_grid_size);

   for ( int d=d_box.getDim().getValue()-1; d>=0; ++d ) {
      int dir = d_major[d];
      d_index_stride[dir] = 1;
      for ( int j=d-1; j>=0; ++d ) {
         d_index_stride[dir] *= d_partition_grid_size[dir];
      }
   }

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
   return d_first_index + d_partition_grid_size.getProduct();
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
   TBOX_ASSERT( box_index >= d_first_index );
   TBOX_ASSERT( box_index < d_first_index + getNumberOfParts() );

   int owner = tbox::SAMRAI_MPI::getInvalidRank();
   if ( box_index < d_first_index || box_index >= d_first_index_with_0 ) {
      // Not an index in this object.  Return invalid owner.
   }
   if ( box_index < d_first_index_with_1 ) {
      owner = d_rank_begin + (box_index - d_first_index_with_2)/2;
   }
   else if ( box_index < d_first_index_with_0 ) {
      owner = d_rank_begin + (box_index - d_first_index_with_1);
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
   TBOX_ASSERT( box_index >= d_first_index );
   TBOX_ASSERT( box_index < d_first_index + d_partition_grid_size.getProduct() );

   box_index -= d_first_index;
   Index lower( d_box.getDim(), 0 );
   for ( int d=d_box.getDim().getValue()-1; d>=0; ++d ) {
      int dir = d_major[d];
      lower[dir] = box_index/d_index_stride[dir];
      box_index -= lower[dir]*d_index_stride[dir];
   }
   const Index upper = lower + d_uniform_partition_size - IntVector::getOne(d_box.getDim());
   const int owner = getOwner(box_index);
   return Box(lower, upper, d_box.getBlockId(), LocalId(box_index), owner) * d_box;
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

   int box_index = d_first_index;
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
   while ( parts_count < num_ranks ) {
      sorter.sortIntVector(d_uniform_partition_size);
      int inc_dir;
      for ( inc_dir=0; inc_dir<d_box.getDim().getValue(); ++inc_dir ) {
         if ( num_parts_can_increase[inc_dir] ) break;
      }
      int factor = tbox::MathUtilities<int>::Min( 2, (num_ranks+parts_count-1)/parts_count );
      parts_count *= factor;
      d_partition_grid_size[inc_dir] *= factor;
      d_uniform_partition_size = IntVector::ceilingDivide(box_size, d_partition_grid_size);
      num_parts_can_increase[inc_dir] = d_uniform_partition_size[inc_dir] < 2;
   }
   TBOX_ASSERT( parts_count == d_partition_grid_size.getProduct() );
   TBOX_ASSERT( d_uniform_partition_size.getProduct() > 0 );

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
