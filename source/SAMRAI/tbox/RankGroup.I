/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   A class to manage groups of processor ranks 
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 ***********************************************************************
 * Default constructor creates RankGroup that uses all available ranks.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
RankGroup::RankGroup():
   d_min(-1),
   d_max(-1),
   d_ranks(0),
   d_storage(USING_ALL),
   d_samrai_mpi(SAMRAI_MPI::getSAMRAIWorld())
{
}

/*
 ***********************************************************************
 * Default constructor creates RankGroup that uses all available ranks in
 * a given communicator
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
RankGroup::RankGroup(
   const SAMRAI_MPI& samrai_mpi):
   d_min(-1),
   d_max(-1),
   d_ranks(0),
   d_storage(USING_ALL),
   d_samrai_mpi(samrai_mpi)
{
}

/*
 ***********************************************************************
 * Virtual destructor
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
RankGroup::~RankGroup()
{
}

/*
 ***********************************************************************
 * Return true if all ranks are used.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
bool RankGroup::containsAllRanks() const
{
   return d_storage == USING_ALL;
}

/*
 ***********************************************************************
 * Set min and max rank values.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void RankGroup::setMinMax(
   const int min,
   const int max)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   int nodes = 1;
   d_samrai_mpi.Comm_size(&nodes);
#endif

   TBOX_ASSERT(min >= 0);
   TBOX_ASSERT(min <= max);

   d_storage = USING_MIN_MAX;
   d_ranks.resizeArray(0);

   d_min = min;
   d_max = max;
}

/*
 ***********************************************************************
 * Determine if given rank is in the group.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
bool RankGroup::isMember(
   const int rank) const
{
   bool is_member = false;

   switch (d_storage) {

      case USING_ALL:
         is_member = true;
         break;

      case USING_ARRAY:
      {
         int lo = 0;
         int hi = d_ranks.size() - 1;
         while (!is_member && hi >= lo) {
            int i = (lo + hi) / 2;
            if (rank < d_ranks[i]) {
               hi = i - 1;
            } else if (rank > d_ranks[i]) {
               lo = i + 1;
            } else {
               is_member = true;
            }
         }
      }
      break;

      case USING_MIN_MAX:
         if (rank >= d_min && rank <= d_max) {
            is_member = true;
         }
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
   }

   return is_member;
}

/*
 ***********************************************************************
 * Size of RankGroup.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
int RankGroup::size() const
{
   int size = -tbox::MathUtilities<int>::getMax();

   switch (d_storage) {

      case USING_ALL:
         d_samrai_mpi.Comm_size(&size);
         break;

      case USING_ARRAY:
         size = d_ranks.size();
         break;

      case USING_MIN_MAX:
         size = d_max - d_min + 1;
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
   }

   return size;
}

/*
 ***********************************************************************
 * Get a rank using 1-to-1 mapping.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
int RankGroup::getMappedRank(
   const int index) const
{
   TBOX_ASSERT(index >= 0 && index < size());

   int mapped_rank = -tbox::MathUtilities<int>::getMax();

   switch (d_storage) {

      case USING_ALL:
         mapped_rank = index;
         break;

      case USING_ARRAY:
         mapped_rank = d_ranks[index];
         break;

      case USING_MIN_MAX:
         mapped_rank = d_min + index;
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
   }

   return mapped_rank;
}

/*
 ***********************************************************************
 * Get an integer identifier from 1-to-1 mapping.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
int RankGroup::getMapIndex(
   const int rank) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   int nodes = 1;
   d_samrai_mpi.Comm_size(&nodes);
   TBOX_ASSERT(rank >= 0 && rank < nodes);
#endif

   TBOX_ASSERT(isMember(rank));

   int map_id = -tbox::MathUtilities<int>::getMax();

   switch (d_storage) {

      case USING_ALL:
         map_id = rank;
         break;

      case USING_ARRAY:
      {
         int lo = 0;
         int hi = d_ranks.size() - 1;
         while (hi >= lo) {
            int i = (lo + hi) / 2;
            if (rank < d_ranks[i]) {
               hi = i - 1;
            } else if (rank > d_ranks[i]) {
               lo = i + 1;
            } else {
               map_id = i;
               break;
            }
         }
      }
      break;

      case USING_MIN_MAX:
         map_id = rank - d_min;
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
   }

   return map_id;
}

}
}
