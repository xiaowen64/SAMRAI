//
// File:        $URL$
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Load balance routines for uniform and non-uniform workloads.
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <cassert>
#endif
#endif

namespace SAMRAI {
namespace mesh {

SAMRAI_INLINE_KEYWORD
void
ChopAndPackLoadBalancer::setIgnoreLevelDomainIsSingleBox(
   bool flag)
{
   d_ignore_level_box_union_is_single_box = flag;
}

SAMRAI_INLINE_KEYWORD
const std::string&
ChopAndPackLoadBalancer::getObjectName() const
{
   return d_object_name;
}

/*
 *************************************************************************
 *
 * Private utility functions to determine parameter values for level.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
ChopAndPackLoadBalancer::getWorkloadDataId(
   int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_number >= 0);
#endif
   int wrk_id = (level_number < d_workload_data_id.getSize() ?
                 d_workload_data_id[level_number] :
                 d_master_workload_data_id);

   return wrk_id;
}

SAMRAI_INLINE_KEYWORD
double
ChopAndPackLoadBalancer::getMaxWorkloadFactor(
   int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_number >= 0);
#endif
   double factor = (level_number < d_max_workload_factor.getSize() ?
                    d_max_workload_factor[level_number] :
                    d_master_max_workload_factor);

   return factor;
}

SAMRAI_INLINE_KEYWORD
double
ChopAndPackLoadBalancer::getWorkloadTolerance(
   int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_number >= 0);
#endif
   double tolerance = (level_number < d_workload_tolerance.getSize() ?
                       d_workload_tolerance[level_number] :
                       d_master_workload_tolerance);

   return tolerance;
}

SAMRAI_INLINE_KEYWORD
std::string
ChopAndPackLoadBalancer::getBinPackMethod(
   int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_number >= 0);
#endif
   std::string factor = (level_number < d_bin_pack_method.getSize() ?
                         d_bin_pack_method[level_number] :
                         d_master_bin_pack_method);

   return factor;
}

}
}
