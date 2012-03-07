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

/*
 *************************************************************************
 *
 * Static member functions that provide linkage with PETSc/SNES package.
 * See header file for SNESAbstractFunctions for more information.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::SNESFuncEval(
   SNES snes,
   Vec x,
   Vec f,
   void* ctx)
{
   NULL_USE(snes);
   ((SNES_SAMRAIContext *)ctx)->getSNESFunctions()->
   evaluateNonlinearFunction(x, f);
   return 0;
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::SNESJacobianTimesVector(
   Mat M,
   Vec xin,
   Vec xout)
{
   void* ctx;

   int ierr = MatShellGetContext(M, &ctx);
   PETSC_SAMRAI_ERROR(ierr);
   return ((SNES_SAMRAIContext *)ctx)->
          getSNESFunctions()->jacobianTimesVector(xin, xout);
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::SNESsetupPreconditioner(
   void* ctx)
{
   Vec current_solution;
   int ierr = SNESGetSolution(((SNES_SAMRAIContext *)ctx)->getSNESSolver(),
         &current_solution);
   PETSC_SAMRAI_ERROR(ierr);
   return ((SNES_SAMRAIContext *)ctx)->
          getSNESFunctions()->setupPreconditioner(current_solution);
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::SNESapplyPreconditioner(
   void* ctx,
   Vec r,
   Vec z)
{
   return ((SNES_SAMRAIContext *)ctx)->
          getSNESFunctions()->applyPreconditioner(r, z);
}

/*
 *************************************************************************
 *
 * Access functions for PETSc objects, and user-supplied functions.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
SNES
SNES_SAMRAIContext::getSNESSolver() const
{
   return d_SNES_solver;
}

SAMRAI_INLINE_KEYWORD
SNESAbstractFunctions*
SNES_SAMRAIContext::getSNESFunctions() const
{
   return d_SNES_functions;
}

SAMRAI_INLINE_KEYWORD
KSP
SNES_SAMRAIContext::getKrylovSolver() const
{
   return d_krylov_solver;
}

SAMRAI_INLINE_KEYWORD
Mat
SNES_SAMRAIContext::getJacobianMatrix() const
{
   return d_jacobian;
}

/*
 *************************************************************************
 *
 *  Access functions for parameters that control solver behavior.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getAbsoluteTolerance() const
{
   return d_absolute_tolerance;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setAbsoluteTolerance(
   double abs_tol)
{
   d_absolute_tolerance = abs_tol;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getRelativeTolerance() const
{
   return d_relative_tolerance;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setRelativeTolerance(
   double rel_tol)
{
   d_relative_tolerance = rel_tol;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getStepTolerance() const
{
   return d_step_tolerance;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setStepTolerance(
   double step_tol)
{
   d_step_tolerance = step_tol;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::getMaxNonlinearIterations() const
{
   return d_maximum_nonlinear_iterations;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setMaxNonlinearIterations(
   int max_nli)
{
   d_maximum_nonlinear_iterations = max_nli;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::getMaxFunctionEvaluations() const
{
   return d_maximum_function_evals;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setMaxFunctionEvaluations(
   int max_feval)
{
   d_maximum_function_evals = max_feval;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
std::string
SNES_SAMRAIContext::getForcingTermStrategy() const
{
   return d_forcing_term_strategy;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setForcingTermStrategy(
   std::string& strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(strategy == "CONSTANT" ||
      strategy == "EWCHOICE1" ||
      strategy == "EWCHOICE2");
#endif
   d_forcing_term_strategy = strategy;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getConstantForcingTerm() const
{
   return d_constant_forcing_term;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setConstantForcingTerm(
   double fixed_eta)
{
   d_constant_forcing_term = fixed_eta;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getInitialForcingTerm() const
{
   return d_initial_forcing_term;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setInitialForcingTerm(
   double initial_eta)
{
   d_initial_forcing_term = initial_eta;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getMaximumForcingTerm() const
{
   return d_maximum_forcing_term;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setMaximumForcingTerm(
   double max_eta)
{
   d_maximum_forcing_term = max_eta;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getEWChoice2Exponent() const
{
   return d_EW_choice2_alpha;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setEWChoice2Exponent(
   double alpha)
{
   d_EW_choice2_alpha = alpha;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getEWChoice2SafeguardExponent() const
{
   return d_EW_safeguard_exponent;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setEWChoice2SafeguardExponent(
   double beta)
{
   d_EW_safeguard_exponent = beta;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getEWChoice2ScaleFactor() const
{
   return d_EW_choice2_gamma;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setEWChoice2ScaleFactor(
   double gamma)
{
   d_EW_choice2_gamma = gamma;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getEWSafeguardThreshold() const
{
   return d_EW_safeguard_disable_threshold;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setEWSafeguardThreshold(
   double threshold)
{
   d_EW_safeguard_disable_threshold = threshold;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
std::string
SNES_SAMRAIContext::getLinearSolverType() const
{
   return d_linear_solver_type;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setLinearSolverType(
   std::string& type)
{
   d_linear_solver_type = type;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
bool
SNES_SAMRAIContext::getUsesPreconditioner() const
{
   return d_uses_preconditioner;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setUsesPreconditioner(
   bool uses_preconditioner)
{
   d_uses_preconditioner = uses_preconditioner;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getLinearSolverAbsoluteTolerance() const
{
   return d_linear_solver_absolute_tolerance;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setLinearSolverAbsoluteTolerance(
   double abs_tol)
{
   d_linear_solver_absolute_tolerance = abs_tol;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getLinearSolverDivergenceTolerance() const
{
   return d_linear_solver_divergence_tolerance;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setLinearSolverDivergenceTolerance(
   double div_tol)
{
   d_linear_solver_divergence_tolerance = div_tol;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::getMaximumLinearIterations() const
{
   return d_maximum_linear_iterations;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setMaximumLinearIterations(
   int max_li)
{
   d_maximum_linear_iterations = max_li;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::getMaximumGMRESKrylovDimension() const
{
   return d_maximum_gmres_krylov_dimension;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setMaximumGMRESKrylovDimension(
   int d)
{
   d_maximum_gmres_krylov_dimension = d;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
std::string
SNES_SAMRAIContext::getGMRESOrthogonalizationMethod() const
{
   return d_gmres_orthogonalization_algorithm;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setGMRESOrthogonalizationMethod(
   std::string& method)
{
   d_gmres_orthogonalization_algorithm = method;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
bool
SNES_SAMRAIContext::getUsesExplicitJacobian() const
{
   return d_uses_explicit_jacobian;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setUsesExplicitJacobian(
   bool use_jac)
{
   d_uses_explicit_jacobian = use_jac;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
std::string
SNES_SAMRAIContext::getDifferencingParameterMethod() const
{
   return d_differencing_parameter_strategy;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setDifferencingParameterMethod(
   std::string& method)
{
   d_differencing_parameter_strategy = method;
   d_context_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
double
SNES_SAMRAIContext::getFunctionEvaluationError() const
{
   return d_function_evaluation_error;
}

SAMRAI_INLINE_KEYWORD
void
SNES_SAMRAIContext::setFunctionEvaluationError(
   double evaluation_error)
{
   d_function_evaluation_error = evaluation_error;
   d_context_needs_initialization = true;
}

/*
 *************************************************************************
 *
 * Get the number of nonlinear iterations used in last solve.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::getNumberOfNonlinearIterations() const
{
   return d_nonlinear_iterations;
}

/*
 *************************************************************************
 *
 * Get the total number of linear iterations used in last solve.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
SNES_SAMRAIContext::getTotalNumberOfLinearIterations() const
{
   int linear_itns;
   int ierr = SNESGetLinearSolveIterations(d_SNES_solver,
         &linear_itns);
   PETSC_SAMRAI_ERROR(ierr);
   return linear_itns;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
SNES_SAMRAIContext::getObjectName() const
{
   return d_object_name;
}

}
}
