/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Wrapper class for KINSOL solver function calls and data.
 */

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Static member functions that provide linkage with KINSOL package.
 * See header file for KINSOLAbstractFunctions for more information.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::KINSOLFuncEval(
   N_Vector soln,
   N_Vector fval,
   void* my_solver)
{
   int success = 0;
   // SGS why no error condition?
   ((KINSOLSolver *)my_solver)->getKINSOLFunctions()->
   evaluateNonlinearFunction(SABSVEC_CAST(soln), SABSVEC_CAST(fval));

   return success;
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::KINSOLPrecondSet(
   N_Vector uu,
   N_Vector uscale,
   N_Vector fval,
   N_Vector fscale,
   void* my_solver,
   N_Vector vtemp1,
   N_Vector vtemp2)
{
   ((KINSOLSolver *)my_solver)->initializeKINSOL();

   int num_feval = 0;
   int success = ((KINSOLSolver *)my_solver)->getKINSOLFunctions()->
      precondSetup(SABSVEC_CAST(uu),
         SABSVEC_CAST(uscale),
         SABSVEC_CAST(fval),
         SABSVEC_CAST(fscale),
         SABSVEC_CAST(vtemp1),
         SABSVEC_CAST(vtemp2),
         num_feval);
   return success;
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::KINSOLPrecondSolve(
   N_Vector uu,
   N_Vector uscale,
   N_Vector fval,
   N_Vector fscale,
   N_Vector vv,
   void* my_solver,
   N_Vector vtemp)

{
   int num_feval = 0;
   int success = ((KINSOLSolver *)my_solver)->getKINSOLFunctions()->
      precondSolve(SABSVEC_CAST(uu),
         SABSVEC_CAST(uscale),
         SABSVEC_CAST(fval),
         SABSVEC_CAST(fscale),
         SABSVEC_CAST(vv),
         SABSVEC_CAST(vtemp),
         num_feval);
   return success;
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::KINSOLJacobianTimesVector(
   N_Vector v,
   N_Vector Jv,
   N_Vector uu,
   int* new_uu,
   void* my_solver)
{
   bool soln_changed = true;
   if (*new_uu == 0) {
      soln_changed = false;
   }

   int success = ((KINSOLSolver *)my_solver)->
      getKINSOLFunctions()->jacobianTimesVector(SABSVEC_CAST(v),
         SABSVEC_CAST(Jv),
         soln_changed,
         SABSVEC_CAST(uu));

   return success;
}

/*
 *************************************************************************
 *
 * Accessory functions for setting user-defined function information.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setKINSOLFunctions(
   KINSOLAbstractFunctions* my_functions,
   const int uses_preconditioner,
   const int uses_jac_times_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(my_functions == (KINSOLAbstractFunctions *)NULL));
#endif

   d_KINSOL_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_uses_jac_times_vector = uses_jac_times_vector;

   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setPreconditioner(
   const int uses_preconditioner)
{
   d_uses_preconditioner = uses_preconditioner;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setJacobianTimesVector(
   const int uses_jac_times_vector)
{
   d_uses_jac_times_vector = uses_jac_times_vector;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
KINSOLAbstractFunctions*
KINSOLSolver::getKINSOLFunctions() const
{
   return d_KINSOL_functions;
}

/*
 *************************************************************************
 *
 * Accessory function for setting constraints for nonlinear system.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setConstraintVector(
   SundialsAbstractVector* constraints)
{
   d_constraints = constraints;
}

/*
 *************************************************************************
 *
 * Accessory function for setting nonlinear solver parameters.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setResidualStoppingTolerance(
   const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(tol >= 0.0);
#endif
   d_residual_tol = tol;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxIterations(
   const int maxits)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxits >= 0);
#endif
   d_max_iter = maxits;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxKrylovDimension(
   const int kdim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(kdim >= 0);
#endif
   d_krylov_dimension = kdim;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setGlobalStrategy(
   const int global)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(global == KIN_NONE || global == KIN_LINESEARCH);
#endif
   d_global_strategy = global;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxNewtonStep(
   const double maxstep)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxstep > 0.0);
#endif
   d_max_newton_step = maxstep;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setNonlinearStepTolerance(
   const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(tol >= 0.0);
#endif
   d_step_tol = tol;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setRelativeFunctionError(
   const double reserr)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(reserr > 0.0);
#endif
   d_relative_function_error = reserr;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setLinearSolverConvergenceTest(
   const int conv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(
      conv == KIN_ETACONSTANT || conv == KIN_ETACHOICE1 || conv ==
      KIN_ETACHOICE2);
#endif
   d_eta_choice = conv;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setEisenstatWalkerParameters(
   const double alpha,
   const double gamma)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(alpha >= 0.0);
   TBOX_ASSERT(gamma >= 0.0);
#endif
   // sgs
   d_eta_alpha = alpha;
   d_eta_gamma = gamma;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setLinearSolverConstantTolerance(
   const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(tol >= 0.0);
#endif
   //
   d_eta_constant = tol;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxSubSetupCalls(
   const int maxsub)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxsub >= 0);
#endif
   d_maxsub = maxsub;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setNoInitialSetup(
   const bool flag)
{
   if (flag) {
      d_no_initial_setup = 1;
   } else {
      d_no_initial_setup = 0;
   }
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setNoResidualMonitoring(
   const bool flag)
{
   if (flag) {
      d_no_residual_monitoring = 1;
   } else {
      d_no_residual_monitoring = 0;
   }
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setResidualMonitoringParams(
   const double omega_min,
   const double omega_max)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(omega_min >= 0);
   TBOX_ASSERT(omega_max >= 0);
   TBOX_ASSERT(omega_max >= omega_min);
#endif
   d_omega_min = omega_min;
   d_omega_max = omega_max;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setResidualMonitoringConstant(
   const double omega)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(omega >= 0);
#endif
   d_omega = omega;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setNoMinEps(
   const bool flag)
{
   if (flag) {
      d_no_min_eps = 1;
   } else {
      d_no_min_eps = 0;
   }
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxBetaFails(
   const int max_beta_fails)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_beta_fails >= 0);
#endif
   d_max_beta_fails = max_beta_fails;
   d_KINSOL_needs_initialization = true;
}

/*
 *************************************************************************
 *
 * Accessory function for setting preconditioner parameters.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxStepsWithNoPrecondSetup(
   const int maxsolv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxsolv > 0);
#endif
   d_max_solves_no_set = maxsolv;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
KINSOLSolver::setMaxLinearSolveRestarts(
   const int restarts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(restarts >= 0);
#endif
   d_max_restarts = restarts;
   d_KINSOL_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::getTotalNumberOfNonlinearIterations() const
{
   long int num;
   int ierr = KINGetNumNonlinSolvIters(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return static_cast<int>(num);
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::getTotalNumberOfFunctionCalls() const
{
   long int num;
   int ierr = KINGetNumFuncEvals(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return static_cast<int>(num);
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::getTotalNumberOfBetaConditionFailures() const
{
   long int num;
   int ierr = KINGetNumBetaCondFails(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return static_cast<int>(num);
}

SAMRAI_INLINE_KEYWORD
int
KINSOLSolver::getTotalNumberOfBacktracks() const
{
   long int num;
   int ierr = KINGetNumBacktrackOps(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return static_cast<int>(num);
}

SAMRAI_INLINE_KEYWORD
double
KINSOLSolver::getScaledResidualNorm() const
{
   realtype norm;
   int ierr = KINGetFuncNorm(d_kin_mem, &norm);
   KINSOL_SAMRAI_ERROR(ierr);
   return norm;
}

SAMRAI_INLINE_KEYWORD
double
KINSOLSolver::getNewtonStepLength() const
{
   realtype step_length;
   int ierr = KINGetStepLength(d_kin_mem, &step_length);
   KINSOL_SAMRAI_ERROR(ierr);
   return step_length;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
KINSOLSolver::getObjectName() const
{
   return d_object_name;
}

}
}
