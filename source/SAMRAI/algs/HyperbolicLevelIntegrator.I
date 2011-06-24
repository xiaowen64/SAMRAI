//
// File:        $URL$
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Integration routines for single level in AMR hierarchy
//              (basic hyperbolic systems).
//

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
const std::string& HyperbolicLevelIntegrator::getObjectName() const
{
   return d_object_name;
}

}
}
