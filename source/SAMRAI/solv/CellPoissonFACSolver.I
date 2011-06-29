/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar poisson equation. 
 *
 ************************************************************************/
namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setDPatchDataId(
   int id) {
   d_poisson_spec.setDPatchDataId(id);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setDConstant(
   double scalar) {
   d_poisson_spec.setDConstant(scalar);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setCPatchDataId(
   int id) {
   d_poisson_spec.setCPatchDataId(id);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setCConstant(
   double scalar) {
// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
   if (scalar == 0.0) {
      d_poisson_spec.setCZero();
   } else {
      d_poisson_spec.setCConstant(scalar);
   }
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setProlongationMethod(
   const std::string& prolongation_method)
{
   d_fac_ops.setProlongationMethod(prolongation_method);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setCoarsestLevelSolverChoice(
   const std::string& choice)
{
   d_fac_ops.setCoarsestLevelSolverChoice(choice);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setCoarsestLevelSolverTolerance(
   double tol)
{
   d_fac_ops.setCoarsestLevelSolverTolerance(tol);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setCoarsestLevelSolverMaxIterations(
   int max_iterations)
{
   d_fac_ops.setCoarsestLevelSolverMaxIterations(max_iterations);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setCoarseFineDiscretization(
   const std::string& coarsefine_method)
{
   d_fac_ops.setCoarseFineDiscretization(coarsefine_method);
}

#ifdef HAVE_HYPRE
SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setUseSMG(
   bool use_smg)
{
   if (d_solver_is_initialized) {
      TBOX_ERROR(
         d_object_name << ": setUseSMG(bool) may NOT be called\n"
         <<
         "while the solver state is initialized, as that\n"
         << "would lead to a corrupted solver state.\n");
   }
   d_fac_ops.setUseSMG(use_smg);
}
#endif

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setPresmoothingSweeps(
   int num_pre_sweeps) {
   d_fac_precond.setPresmoothingSweeps(num_pre_sweeps);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setPostsmoothingSweeps(
   int num_post_sweeps) {
   d_fac_precond.setPostsmoothingSweeps(num_post_sweeps);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setMaxCycles(
   int max_cycles) {
   d_fac_precond.setMaxCycles(max_cycles);
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::setResidualTolerance(
   double residual_tol) {
   d_fac_precond.setResidualTolerance(residual_tol);
}

SAMRAI_INLINE_KEYWORD
int CellPoissonFACSolver::getNumberOfIterations() const
{
   return d_fac_precond.getNumberOfIterations();
}

SAMRAI_INLINE_KEYWORD
double CellPoissonFACSolver::getResidualNorm() const
{
   return d_fac_precond.getResidualNorm();
}

SAMRAI_INLINE_KEYWORD
void CellPoissonFACSolver::getConvergenceFactors(
   double& avg_factor,
   double& final_factor)
const
{
   d_fac_precond.getConvergenceFactors(avg_factor, final_factor);
}

SAMRAI_INLINE_KEYWORD
const std::string& CellPoissonFACSolver::getObjectName() const {
   return d_object_name;
}

}
}
