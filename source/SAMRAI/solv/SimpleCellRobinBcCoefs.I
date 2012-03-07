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
void
SimpleCellRobinBcCoefs::setDiffusionCoefId(
   int diffusion_coef_id)
{
   d_diffusion_coef_id = diffusion_coef_id;
   d_diffusion_coef_constant = 0.0;
}

SAMRAI_INLINE_KEYWORD
void
SimpleCellRobinBcCoefs::setDiffusionCoefConstant(
   double diffusion_coef_constant)
{
   d_diffusion_coef_constant = diffusion_coef_constant;
   d_diffusion_coef_id = -1;
}

SAMRAI_INLINE_KEYWORD
const std::string&
SimpleCellRobinBcCoefs::getObjectName() const
{
   return d_object_name;
}

}
}
