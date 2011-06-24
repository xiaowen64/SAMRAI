/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: KINSOL solver for use within a SAMRAI-based application.
 */

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
const std::string& KINSOL_SAMRAIContext::getObjectName() const {
   return d_object_name;
}

}
}
