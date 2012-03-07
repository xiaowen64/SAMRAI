/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: KINSOL solver for use within a SAMRAI-based application.
 */

namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
const std::string&
KINSOL_SAMRAIContext::getObjectName() const
{
   return d_object_name;
}

SAMRAI_INLINE_KEYWORD
KINSOLSolver*
KINSOL_SAMRAIContext::getKINSOLSolver()
{
   return d_KINSOL_solver;
}

}
}
