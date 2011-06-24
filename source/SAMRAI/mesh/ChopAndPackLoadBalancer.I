//
// File:        $URL$
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2010 Lawrence Livermore National Security, LLC
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
const std::string& ChopAndPackLoadBalancer::getObjectName() const
{
   return d_object_name;
}

}
}
