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
const std::string&
CartesianRobinBcHelper::getObjectName() const
{
   return d_object_name;
}

/*
 ************************************************************************
 * Set the coefficient strategy pointer that will be used to get
 * Robin bc coefficients.  It should be some external implementation.
 * This function implies that the simple mappings for
 * parallelpiped domains are not used and resets those arrays to null.
 * is a parallelpiped (not checked) and that the boundary condition
 * coefficients are functions only of the location index of the
 * boundary.
 ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
CartesianRobinBcHelper::setCoefImplementation(
   const RobinBcCoefStrategy* coef_strategy)
{
   if (!coef_strategy) {
      TBOX_ERROR(d_object_name << ": Invalid pointer value"
                               << std::endl);
   }
   d_coef_strategy = coef_strategy;
}

SAMRAI_INLINE_KEYWORD
void
CartesianRobinBcHelper::setTargetDataId(
   int target_data_id)
{
   d_target_data_id = target_data_id;
}

SAMRAI_INLINE_KEYWORD
void
CartesianRobinBcHelper::setHomogeneousBc(
   bool is_homogeneous)
{
   d_homogeneous_bc = is_homogeneous;
}

}
}
