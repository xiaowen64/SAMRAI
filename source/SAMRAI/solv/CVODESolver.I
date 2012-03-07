/*
 * File:        $URL$
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision$
 * Modified:    $LastChangedDate$
 * Description: Wrapper class for CVODE solver function calls and data.
 */

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Static member functions that provide linkage with CVODE package.
 * See header file for CVODEAbstractFunctions for more information.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
CVODESolver::CVODERHSFuncEval(
   realtype t,
   N_Vector y,
   N_Vector y_dot,
   void* my_solver)
{
   return ((CVODESolver *)my_solver)->getCVODEFunctions()->
          evaluateRHSFunction(t, SABSVEC_CAST(y), SABSVEC_CAST(y_dot));
}

/*
 *************************************************************************
 *
 * Static member functions that provide linkage with CVSpgmr package.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
CVODESolver::CVSpgmrPrecondSet(
   realtype t,
   N_Vector y,
   N_Vector fy,
   int jok,
   booleantype* jcurPtr,
   realtype gamma,
   void* my_solver,
   N_Vector vtemp1,
   N_Vector vtemp2,
   N_Vector vtemp3)
{
   int success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
      CVSpgmrPrecondSet(t,
         SABSVEC_CAST(y),
         SABSVEC_CAST(fy),
         jok,
         jcurPtr,
         gamma,
         SABSVEC_CAST(vtemp1),
         SABSVEC_CAST(vtemp2),
         SABSVEC_CAST(vtemp3));

   return success;
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::CVSpgmrPrecondSolve(
   realtype t,
   N_Vector y,
   N_Vector fy,
   N_Vector r,
   N_Vector z,
   realtype gamma,
   realtype delta,
   int lr,
   void* my_solver,
   N_Vector vtemp)
{
   int success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
      CVSpgmrPrecondSolve(t,
         SABSVEC_CAST(y),
         SABSVEC_CAST(fy),
         SABSVEC_CAST(r),
         SABSVEC_CAST(z),
         gamma,
         delta,
         lr,
         SABSVEC_CAST(vtemp));

   return success;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
CVODESolver::initialize(
   SundialsAbstractVector* solution)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(solution == (SundialsAbstractVector *)NULL));
   TBOX_ASSERT(d_solution_vector == (SundialsAbstractVector *)NULL);
#endif
   d_solution_vector = solution;
   d_CVODE_needs_initialization = true;
   initializeCVODE();
}

/*
 *************************************************************************
 *
 * Integrate system of ODEs to d_t_f.  If necessary, re-initialize
 * CVODE.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
CVODESolver::solve()
{
   initializeCVODE();

   /*
    * Check to make sure that user specified final value for t
    * is greater than initial value for t.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_user_t_f > d_t_0);
#endif

   /*
    * See cvode.h header file for definition of return types.
    */

   int retval = CVode(d_cvode_mem,
         d_user_t_f,
         d_solution_vector->getNVector(),
         &d_actual_t_f,
         d_stepping_method);

   return retval;

}

/*
 *************************************************************************
 *
 * Setting CVODE log file name and print flag for CVODE statistics.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setLogFileData(
   const std::string& log_fname)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!log_fname.empty());
#endif
   if (!(log_fname == d_cvode_log_file_name)) {
      d_cvode_log_file_name = log_fname;
      d_CVODE_needs_initialization = true;
   }
}

/*
 *************************************************************************
 *
 * Accessor functions for user-defined function and linear solver.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setCVODEFunctions(
   CVODEAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(my_functions == (CVODEAbstractFunctions *)NULL));
#endif

   d_cvode_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
CVODEAbstractFunctions*
CVODESolver::getCVODEFunctions() const
{
   return d_cvode_functions;
}
/*
 *************************************************************************
 *
 * Accessor functions for CVODE integration parameters.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
CVODESolver::setLinearMultistepMethod(
   int linear_multistep_method)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((linear_multistep_method == CV_ADAMS) ||
      (linear_multistep_method == CV_BDF));
#endif
   d_linear_multistep_method = linear_multistep_method;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setIterationType(
   int iteration_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((iteration_type == CV_FUNCTIONAL) ||
      (iteration_type == CV_NEWTON));
#endif
   d_iteration_type = iteration_type;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setToleranceType(
   int tolerance_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((tolerance_type == CV_SS) ||
      (tolerance_type == CV_SV));
#endif
   d_tolerance_type = tolerance_type;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setRelativeTolerance(
   double relative_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(relative_tolerance >= 0.0);
#endif

   d_relative_tolerance = relative_tolerance;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setAbsoluteTolerance(
   double absolute_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(absolute_tolerance >= 0.0);
#endif
   d_absolute_tolerance_scalar = absolute_tolerance;
   d_use_scalar_absolute_tolerance = true;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setAbsoluteTolerance(
   SundialsAbstractVector* absolute_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(absolute_tolerance == (SundialsAbstractVector *)NULL));
   TBOX_ASSERT(absolute_tolerance->vecMin() >= 0.0);
#endif
   d_absolute_tolerance_vector = absolute_tolerance;
   d_use_scalar_absolute_tolerance = false;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setSteppingMethod(
   int stepping_method)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((stepping_method == CV_NORMAL) ||
      (stepping_method == CV_ONE_STEP));
#endif
   d_stepping_method = stepping_method;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setInitialValueOfIndependentVariable(
   double t_0)
{
   d_t_0 = t_0;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setFinalValueOfIndependentVariable(
   double t_f,
   bool cvode_needs_initialization)
{
   d_user_t_f = t_f;
   d_CVODE_needs_initialization = cvode_needs_initialization;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setInitialConditionVector(
   SundialsAbstractVector* ic_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(ic_vector == (SundialsAbstractVector *)NULL));
#endif
   d_ic_vector = ic_vector;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setMaximumLinearMultistepMethodOrder(
   int max_order)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_order >= 0);
#endif
   d_max_order = max_order;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setMaximumNumberOfInternalSteps(
   int max_num_internal_steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_num_internal_steps >= 0);
#endif

   d_max_num_internal_steps = max_num_internal_steps;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setMaximumNumberOfNilStepWarnings(
   int max_num_warnings)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_num_warnings >= 0);
#endif

   d_max_num_warnings = max_num_warnings;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setInitialStepSize(
   double init_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(init_step_size >= 0.0);
#endif
   d_init_step_size = init_step_size;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setMaximumAbsoluteStepSize(
   double max_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_step_size >= 0.0);
#endif
   d_max_step_size = max_step_size;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setMinimumAbsoluteStepSize(
   double min_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(min_step_size >= 0.0);
#endif
   d_min_step_size = min_step_size;
   d_CVODE_needs_initialization = true;
}

/*
 *************************************************************************
 *
 * Accessor functions for CVSpgmr parameters.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
CVODESolver::setPreconditioningType(
   int precondition_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((precondition_type == PREC_NONE) ||
      (precondition_type == PREC_LEFT) ||
      (precondition_type == PREC_RIGHT) ||
      (precondition_type == PREC_BOTH));
#endif
   d_precondition_type = precondition_type;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setGramSchmidtType(
   int gs_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((gs_type == CLASSICAL_GS) ||
      (gs_type == MODIFIED_GS));
#endif
   d_gram_schmidt_type = gs_type;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setMaxKrylovDimension(
   int max_krylov_dim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_krylov_dim >= 0);
#endif
   d_max_krylov_dim = max_krylov_dim;
   d_CVODE_needs_initialization = true;
}

SAMRAI_INLINE_KEYWORD
void
CVODESolver::setCVSpgmrToleranceScaleFactor(
   double tol_scale_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(tol_scale_factor >= 0);
#endif
   d_tol_scale_factor = tol_scale_factor;
   d_CVODE_needs_initialization = true;
}

/*
 *************************************************************************
 *
 * Accessor functions for results of CVODE integration step.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
SundialsAbstractVector*
CVODESolver::getSolutionVector() const
{
   return d_solution_vector;
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getDkyVector(
   double t,
   int k,
   SundialsAbstractVector* dky) const
{
   int return_code = CVodeGetDky(d_cvode_mem, t, k, dky->getNVector());
   return return_code;
}

SAMRAI_INLINE_KEYWORD
double
CVODESolver::getActualFinalValueOfIndependentVariable() const
{
   return d_actual_t_f;
}

/*
 *************************************************************************
 *
 * Access methods for CVODE statistics.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
CVODESolver::printStatistics(
   std::ostream& os) const
{
   printCVODEStatistics(os);
   printCVSpgmrStatistics(os);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfInternalStepsTaken() const
{
   long int r;
   int ierr = CVodeGetNumSteps(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfRHSFunctionCalls() const
{
   long int r;
   int ierr = CVodeGetNumRhsEvals(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfLinearSolverSetupCalls() const
{
   long int r;
   int ierr = CVodeGetNumLinSolvSetups(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfNewtonIterations() const
{
   long int r;
   int ierr = CVodeGetNumNonlinSolvIters(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfNonlinearConvergenceFailures() const
{
   long int r;
   int ierr = CVodeGetNumNonlinSolvConvFails(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfLocalErrorTestFailures() const
{
   long int r;
   int ierr = CVodeGetNumErrTestFails(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getOrderUsedDuringLastInternalStep() const
{
   int r;
   int ierr = CVodeGetLastOrder(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getOrderToBeUsedDuringNextInternalStep() const
{
   int r;
   int ierr = CVodeGetCurrentOrder(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getCVODEMemoryUsageForDoubles() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r1);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getCVODEMemoryUsageForIntegers() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r2);
}

SAMRAI_INLINE_KEYWORD
double
CVODESolver::getStepSizeForLastInternalStep() const
{
   realtype r;
   int ierr = CVodeGetLastStep(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

SAMRAI_INLINE_KEYWORD
double
CVODESolver::getStepSizeForNextInternalStep() const
{
   realtype r;
   int ierr = CVodeGetCurrentStep(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

SAMRAI_INLINE_KEYWORD
double
CVODESolver::getCurrentInternalValueOfIndependentVariable() const
{
   realtype r;
   int ierr = CVodeGetCurrentStep(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

SAMRAI_INLINE_KEYWORD
double
CVODESolver::getCVODESuggestedToleranceScalingFactor() const
{
   realtype r;
   int ierr = CVodeGetTolScaleFactor(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

/*
 *************************************************************************
 *
 * Access methods for CVSpgmr statistics.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfPreconditionerEvaluations() const
{
   long int r;
   int ierr = CVSpilsGetNumPrecEvals(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfLinearIterations() const
{
   long int r;
   int ierr = CVSpilsGetNumLinIters(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfPrecondSolveCalls() const
{
   long int r;
   int ierr = CVSpilsGetNumPrecSolves(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getNumberOfLinearConvergenceFailures() const
{
   long int r;
   int ierr = CVSpilsGetNumConvFails(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getCVSpgmrMemoryUsageForDoubles() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r1);
}

SAMRAI_INLINE_KEYWORD
int
CVODESolver::getCVSpgmrMemoryUsageForIntegers() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r2);
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string&
CVODESolver::getObjectName() const
{
   return d_object_name;
}

}
}
