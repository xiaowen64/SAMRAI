/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Robin boundary condition problem-dependent interfaces.
 */

namespace SAMRAI {
namespace solv {

/*
 ************************************************************************
 * Set the boundary value for a Dirichlet boundary condition.
 ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
LocationIndexRobinBcCoefs::setBoundaryValue(
   int location_index,
   double value)
{
   if (location_index < 0 || location_index >= 2 * d_dim.getValue()) {
      TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                      << "in [0," << 2 * d_dim.getValue() - 1 << "].\n");
   }
   d_a_map[location_index] = 1.0;
   d_b_map[location_index] = 0.0;
   d_g_map[location_index] = value;
}

/*
 ************************************************************************
 * Set the slpe for a Neumann boundary condition.
 ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
LocationIndexRobinBcCoefs::setBoundarySlope(
   int location_index,
   double slope)
{
   if (location_index >= 2 * d_dim.getValue()) {
      TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                      << "in [0," << 2 * d_dim.getValue() - 1 << "].\n");
   }
   d_a_map[location_index] = 0.0;
   d_b_map[location_index] = 1.0;
   d_g_map[location_index] = slope;
}

/*
 ************************************************************************
 * Set the raw bc coefficients.
 ************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
LocationIndexRobinBcCoefs::setRawCoefficients(
   int location_index,
   double a,
   double b,
   double g)
{
   if (location_index >= 2 * d_dim.getValue()) {
      TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                      << "in [0," << 2 * d_dim.getValue() - 1 << "].\n");
   }
   d_a_map[location_index] = a;
   d_b_map[location_index] = b;
   d_g_map[location_index] = g;
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
LocationIndexRobinBcCoefs::getCoefficients(
   int i,
   double& a,
   double& b,
   double& g) const
{
   a = d_a_map[i];
   b = d_b_map[i];
   g = d_g_map[i];
}

/*
 ************************************************************************
 ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
LocationIndexRobinBcCoefs::getObjectName() const
{
   return d_object_name;
}

}
}
