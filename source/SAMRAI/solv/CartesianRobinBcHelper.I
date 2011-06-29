/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Robin boundary condition support on cartesian grids.
 */

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
const std::string& CartesianRobinBcHelper::getObjectName() const {
   return d_object_name;
}

}
}
