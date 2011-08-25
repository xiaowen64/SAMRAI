/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A class to manage groups of processor ranks
 *
 ************************************************************************/

#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/List.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/tbox/RankGroup.I"
#endif

namespace SAMRAI {
namespace tbox {

/*
 ***********************************************************************
 * Constructor that takes a min and max rank.
 ***********************************************************************
 */

RankGroup::RankGroup(
   const int min,
   const int max,
   const SAMRAI_MPI& samrai_mpi):
   d_min(min),
   d_max(max),
   d_ranks(0),
   d_storage(USING_MIN_MAX),
   d_samrai_mpi(samrai_mpi)
{
   int nodes = 1;
   samrai_mpi.Comm_size(&nodes);

   TBOX_ASSERT(min >= 0);
   TBOX_ASSERT(max < nodes);
   TBOX_ASSERT(min <= max);

   /*
    * If min and max cover the full set of ranks, then switch to
    * "using_all" mode.
    */
   if (min == 0 && max == nodes - 1) {
      d_min = -1;
      d_max = -1;
      d_storage = USING_ALL;
   }
}

/*
 ***********************************************************************
 * Constructor that takes an array of ranks.
 ***********************************************************************
 */

RankGroup::RankGroup(
   const tbox::Array<int>& rank_group,
   const SAMRAI_MPI& samrai_mpi):
   d_min(-1),
   d_max(-1),
   d_ranks(rank_group),
   d_storage(USING_ARRAY),
   d_samrai_mpi(samrai_mpi)
{
   TBOX_ASSERT(rank_group.size() > 0);

   int nodes = 1;
   samrai_mpi.Comm_size(&nodes);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(rank_group.size() <= nodes);

   /*
    * Check that each entry in the array has a unique value and is increasing
    * order
    */
   for (int i = 0; i < rank_group.size(); i++) {
      TBOX_ASSERT(rank_group[i] >= 0);
      TBOX_ASSERT(rank_group[i] < nodes);
      if (i > 0) {
         TBOX_ASSERT(rank_group[i] > rank_group[i - 1]);
      }
   }
#endif

   /*
    * If array is the full set of ranks, then switch to "using_all" mode.
    */
   if (rank_group.size() == nodes) {
      d_ranks.resizeArray(0);
      d_storage = USING_ALL;
   }
}

}
}
