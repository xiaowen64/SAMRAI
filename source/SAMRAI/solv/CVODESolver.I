/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Wrapper class for CVODE solver function calls and data.
 */

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
const std::string& CVODESolver::getObjectName() const {
   return d_object_name;
}

}
}
