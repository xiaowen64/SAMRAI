/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A class to manage groups of processor ranks
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 ***********************************************************************
 * Return true if all ranks are used.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
bool
RankGroup::containsAllRanks() const
{
   return d_storage == USING_ALL;
}

/*
 ***********************************************************************
 * Set min and max rank values.
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
RankGroup::setMinMax(
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

}
}
