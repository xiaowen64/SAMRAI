/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Utilitiy for breaking boxes during partitioning.
 *
 ************************************************************************/

#ifndef included_mesh_PartitioningParams_C
#define included_mesh_PartitioningParams_C

#include "SAMRAI/mesh/PartitioningParams.h"

namespace SAMRAI {
namespace mesh {



/*
 *************************************************************************
 *************************************************************************
 */

bool
PartitioningParams::compareLoads(
   int flags[],
   double cur_load,
   double new_load,
   double ideal_load,
   double low_load,
   double high_load ) const
{
   double cur_range_miss = cur_load >= high_load ? cur_load-high_load :
      cur_load <= low_load ? low_load-cur_load : 0.0;
   double new_range_miss = new_load >= high_load ? new_load-high_load :
      new_load <= low_load ? low_load-new_load : 0.0;
   flags[0] = new_range_miss < (cur_range_miss-d_load_comparison_tol) ?
      1 : new_range_miss > cur_range_miss ? -1 : 0;

   double cur_diff = tbox::MathUtilities<double>::Abs(cur_load-ideal_load);
   double new_diff = tbox::MathUtilities<double>::Abs(new_load-ideal_load);

   flags[1] = new_diff < (cur_diff-d_load_comparison_tol) ?
      1 : new_diff > cur_diff ? -1 : 0;

   flags[2] = flags[0] != 0 ? flags[0] : flags[1] != 0 ? flags[1] : 0;

   return flags[2] == 1;
}



}
}


#endif
