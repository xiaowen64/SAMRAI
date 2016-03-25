/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/CellPoissonHypreSolver.I $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2043 $
 * Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
 * Description: Level solver for diffusion-like elliptic problems.
 */


#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace solv {


template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setSolnIdDepth( const int depth ) {
   d_soln_depth = depth;
   return;
}

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setRhsIdDepth( const int depth ) {
   d_rhs_depth = depth;
   return;
}




template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setUseSMG( bool use_smg ) {
   d_use_smg = use_smg;
   return;
}




/*
********************************************************************
* Specify bc using the default internal bc coefficient object.     *
* Clear up data supporting external bc coefficient setter.         *
********************************************************************
*/

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setBoundaries(
   const std::string& boundary_type,
   const int fluxes,
   const int flags,
   int* bdry_types)
{
   d_physical_bc_simple_case.setBoundaries(boundary_type,
                                           fluxes,
                                           flags,
                                           bdry_types);
   d_physical_bc_coef_strategy = &d_physical_bc_simple_case;
   d_physical_bc_variable.setNull();
   return;
}



/*
********************************************************************
* Set the physical boundary condition object.                      *
********************************************************************
*/

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setPhysicalBcCoefObject(
   const RobinBcCoefStrategy<DIM> *physical_bc_coef_strategy,
   const tbox::Pointer< hier::Variable<DIM> > variable )
{
   d_physical_bc_coef_strategy = physical_bc_coef_strategy;
   d_physical_bc_variable = variable;
   return;
}




template<int DIM> inline
int CellPoissonHypreSolver<DIM>::getNumberOfIterations() const
{
   return(d_number_iterations);
}

template<int DIM> inline
double CellPoissonHypreSolver<DIM>::getRelativeResidualNorm() const
{
   return(d_relative_residual_norm);
}

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setNumPreRelaxSteps(const int steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
#endif
   d_num_pre_relax_steps = steps;
}

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setNumPostRelaxSteps(const int steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
#endif
   d_num_post_relax_steps = steps;
}

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setPrintSolverInfo(const bool print)
{
   d_print_solver_info = print;
}

template<int DIM> inline
void CellPoissonHypreSolver<DIM>::setStoppingCriteria(
   const int max_iterations,
   const double residual_tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_iterations >= 0);
   TBOX_ASSERT(residual_tol      >= 0.0);
#endif
   d_max_iterations = max_iterations;
   d_relative_residual_tol      = residual_tol;
}


}
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
