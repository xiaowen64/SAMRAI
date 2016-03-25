/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/CellPoissonFACSolver.I $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 2043 $
 * Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
 * Description: High-level solver (wrapper) for scalar poisson equation.
 */


#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace solv {


template<int DIM> inline
void CellPoissonFACSolver<DIM>::setDPatchDataId( int id ) {
   d_poisson_spec.setDPatchDataId(id);
   return;
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setDConstant( double scalar ) {
   d_poisson_spec.setDConstant(scalar);
   return;
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setCPatchDataId( int id ) {
   d_poisson_spec.setCPatchDataId(id);
   return;
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setCConstant( double scalar ) {
   if ( scalar == 0.0 ) {
      d_poisson_spec.setCZero();
   }
   else {
      d_poisson_spec.setCConstant(scalar);
   }
   return;
}





template<int DIM> inline
void CellPoissonFACSolver<DIM>::setProlongationMethod(
   const std::string &prolongation_method )
{
   d_fac_ops.setProlongationMethod(prolongation_method);
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setCoarsestLevelSolverChoice(
   const std::string &choice )
{
   d_fac_ops.setCoarsestLevelSolverChoice(choice);
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setCoarsestLevelSolverTolerance(
   double tol )
{
   d_fac_ops.setCoarsestLevelSolverTolerance(tol);
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setCoarsestLevelSolverMaxIterations(
   int max_iterations )
{
   d_fac_ops.setCoarsestLevelSolverMaxIterations(max_iterations);
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setCoarseFineDiscretization(
   const std::string &coarsefine_method )
{
   d_fac_ops.setCoarseFineDiscretization(coarsefine_method);
}

#ifdef HAVE_HYPRE
template<int DIM> inline
void CellPoissonFACSolver<DIM>::setUseSMG( bool use_smg )
{
   if ( d_solver_is_initialized ) {
      TBOX_ERROR(d_object_name << ": setUseSMG(bool) may NOT be called\n"
                 << "while the solver state is initialized, as that\n"
                 << "would lead to a corrupted solver state.\n");
   }
   d_fac_ops.setUseSMG(use_smg);
   return;
}
#endif







template<int DIM> inline
void CellPoissonFACSolver<DIM>::setPresmoothingSweeps( int num_pre_sweeps ) {
   d_fac_precond.setPresmoothingSweeps( num_pre_sweeps );
   return;
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setPostsmoothingSweeps( int num_post_sweeps ) {
   d_fac_precond.setPostsmoothingSweeps( num_post_sweeps );
   return;
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setMaxCycles( int max_cycles ) {
   d_fac_precond.setMaxCycles( max_cycles );
   return;
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::setResidualTolerance( double residual_tol ) {
   d_fac_precond.setResidualTolerance( residual_tol );
   return;
}





template<int DIM> inline
int CellPoissonFACSolver<DIM>::getNumberOfIterations() const
{
   return d_fac_precond.getNumberOfIterations();
}

template<int DIM> inline
double CellPoissonFACSolver<DIM>::getResidualNorm() const
{
   return d_fac_precond.getResidualNorm();
}

template<int DIM> inline
void CellPoissonFACSolver<DIM>::getConvergenceFactors(double &avg_factor,
                                                       double &final_factor)
   const
{
  d_fac_precond.getConvergenceFactors(avg_factor,final_factor);
  return;
}





}
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
