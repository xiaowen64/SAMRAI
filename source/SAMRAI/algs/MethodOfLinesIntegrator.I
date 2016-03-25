//
// File:        $URL$
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Basic method-of-lines time integration algorithm.
//

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
const std::string& MethodOfLinesIntegrator::getObjectName() const
{
   return d_object_name;
}

}
}
