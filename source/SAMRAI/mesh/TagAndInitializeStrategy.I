//
// File:        $URL$
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Strategy interface for params, tagging, init for gridding.
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
const std::string& TagAndInitializeStrategy::getObjectName() const
{
   return d_object_name;
}

}
}
