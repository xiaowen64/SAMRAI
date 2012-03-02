//
// File:        $URL$
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
// Revision:    \f$       \f$
// Modified:    $LastChangedDate$
// Description: Singleton manager class for statistic objects.
//

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
const std::string&
StatisticRestartDatabase::getObjectName() const
{
   return d_object_name;
}

/*
 *************************************************************************
 *
 * Utility functions for retrieving statistic data from statistician.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
Statistician::getNumberProcessorStats() const
{
   return d_num_proc_stats;
}

SAMRAI_INLINE_KEYWORD
int
Statistician::getNumberPatchStats() const
{
   return d_num_patch_stats;
}

SAMRAI_INLINE_KEYWORD
int
Statistician::getMaximumNumberOfStatistics() {
   return d_proc_statistics.getSize();
}

SAMRAI_INLINE_KEYWORD
void
Statistician::setMaximumNumberOfStatistics(
   const int size) {
   if (size > d_proc_statistics.getSize()) {
      d_proc_statistics.resizeArray(size);
      d_patch_statistics.resizeArray(size);
   }
}

}
}
