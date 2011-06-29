/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Level solver for diffusion-like elliptic problems.
 */

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
const std::string& SimpleCellRobinBcCoefs::getObjectName() const {
   return d_object_name;
}

}
}
