/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Wrapper for SNES solver for use in a SAMRAI-based application.
 */

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
const std::string& SNES_SAMRAIContext::getObjectName() const {
   return d_object_name;
}

}
}
