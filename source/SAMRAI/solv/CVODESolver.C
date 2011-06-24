/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   C++ Wrapper class for CVODE solver package 
 *
 ************************************************************************/

#include "SAMRAI/solv/CVODESolver.h"

#ifdef HAVE_SUNDIALS

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

// CVODE includes
extern "C" {
#include "cvode/cvode.h"
}

extern "C" {
#include "cvode/cvode_spgmr.h"
}

#ifndef SAMRAI_INLINE
#include "SAMRAI/solv/CVODESolver.I"
#endif

namespace SAMRAI {
namespace solv {

const int CVODESolver::STAT_OUTPUT_BUFFER_SIZE = 256;

#define SABSVEC_CAST(v) \
   (static_cast<SAMRAI::solv::SundialsAbstractVector *>(v \
                                                        -> \
                                                        content))

/*
 *************************************************************************
 *                                                                       *
 * Static member functions that provide linkage with CVODE package.      *
 * See header file for CVODEAbstractFunctions for more information. *
 *                                                                       *
 *************************************************************************
 */

int CVODESolver::CVODERHSFuncEval(
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
 *                                                                       *
 * Static member functions that provide linkage with CVSpgmr package.	*
 *                                                                       *
 *************************************************************************
 */

int CVODESolver::CVSpgmrPrecondSet(
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

   int success;
   success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
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

int CVODESolver::CVSpgmrPrecondSolve(
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

   int success;
   success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
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
 *                                                                       *
 * CVODESolver constructor and destructor.                         *
 *                                                                       *
 *************************************************************************
 */
CVODESolver::CVODESolver(
   const std::string& object_name,
   CVODEAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!(my_functions == (CVODEAbstractFunctions *)NULL));
#endif

   d_object_name = object_name;
   d_cvode_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;

   d_solution_vector = (SundialsAbstractVector *)NULL;

   /*
    * Set default parameters to safe values or to CVODE/CVSpgmr defaults.
    */

   /*
    * CVODE memory record and log file.
    */
   d_cvode_mem = NULL;
   d_cvode_log_file = NULL;
   d_cvode_log_file_name = "cvode.log";

   /*
    * ODE parameters.
    */
   d_t_0 = 0.0;
   d_user_t_f = 0.0;
   d_actual_t_f = 0.0;
   d_ic_vector = ((SundialsAbstractVector *)NULL);

   /*
    * ODE integration parameters.
    */

   setLinearMultistepMethod(CV_BDF);
   setIterationType(CV_FUNCTIONAL);
   setToleranceType(CV_SS);
   setRelativeTolerance(0.0);
   setAbsoluteTolerance(0.0);
   d_absolute_tolerance_vector = (SundialsAbstractVector *)NULL;
   setSteppingMethod(CV_NORMAL);

   d_max_order = -1;
   d_max_num_internal_steps = -1;
   d_max_num_warnings = -1;
   d_init_step_size = -1;
   d_max_step_size = -1;
   d_min_step_size = -1;

   /*
    * CVSpgmr parameters.
    *
    * Note that when the maximum krylov dimension and CVSpgmr
    * tolerance scale factor are set to 0, CVSpgmr uses its
    * internal default values.  These are described in the header for
    * this class.
    */
   setPreconditioningType(PREC_NONE);
   setGramSchmidtType(MODIFIED_GS);
   setMaxKrylovDimension(0);
   setCVSpgmrToleranceScaleFactor(0);

   d_CVODE_needs_initialization = true;
}

CVODESolver::~CVODESolver()
{
   if (d_cvode_log_file) fclose(d_cvode_log_file);
   if (d_cvode_mem) CVodeFree(&d_cvode_mem);
}

/*
 *************************************************************************
 *                                                                       *
 * Functions to initialize linear solver and reset CVODE structure.      *
 *                                                                       *
 *************************************************************************
 */

void CVODESolver::initialize(
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

void CVODESolver::initializeCVODE()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(d_solution_vector == (SundialsAbstractVector *)NULL));
#endif

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   if (d_CVODE_needs_initialization) {

      /*
       * Set CVODE log file.
       */
      if (d_cvode_log_file) {
         fclose(d_cvode_log_file);
      }
      d_cvode_log_file = fopen(d_cvode_log_file_name.c_str(), "w");

      /*
       * Make sure that either the relative tolerance or the
       * absolute tolerance has been set to a nonzero value.
       */
      bool tolerance_error = false;
      if (d_use_scalar_absolute_tolerance) {
         if ((d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_scalar == 0.0)) {
            tolerance_error = true;
         }
      } else {
         if ((d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_vector->maxNorm() == 0.0)) {
            tolerance_error = true;
         }
      }

      if (tolerance_error && d_cvode_log_file) {
         fprintf(d_cvode_log_file,
            "%s: Both relative and absolute tolerance have value 0.0",
            d_object_name.c_str());
      }

      /*
       * CVODE function pointer.
       */
      CVRhsFn RHSFunc = CVODESolver::CVODERHSFuncEval;

      /*
       * Free previously allocated CVode memory.  Note that the
       * CVReInit() function is not used since the d_neq variable
       * might have been changed from the previous initializeCVODE()
       * call.
       */
      if (d_cvode_mem) CVodeFree(&d_cvode_mem);

      /*
       * Allocate main memory for CVODE package.
       */

      d_cvode_mem = CVodeCreate(d_linear_multistep_method, d_iteration_type);

      /*
       * Set tolerance parameter based on type
       */
      void* abstol;
      if (d_tolerance_type == CV_SV) {
         abstol = d_absolute_tolerance_vector->getNVector();
      } else {
         abstol = &d_absolute_tolerance_scalar;
      }

      int ierr = CVodeMalloc(d_cvode_mem,
            RHSFunc,
            d_t_0,
            d_ic_vector->getNVector(),
            d_tolerance_type,
            d_relative_tolerance,
            abstol);
      CVODE_SAMRAI_ERROR(ierr);

      ierr = CVodeSetFdata(d_cvode_mem, this);
      CVODE_SAMRAI_ERROR(ierr);

      /*
       * If the iteration type is set to NEWTON, then initialize
       * the CVSpgmr linear solver.
       */
      if (d_iteration_type == CV_NEWTON) {

         ierr = CVSpgmr(d_cvode_mem,
               d_precondition_type,
               d_max_krylov_dim);
         CVODE_SAMRAI_ERROR(ierr);

         if (!(d_max_order < 1)) {
            ierr = CVodeSetMaxOrd(d_cvode_mem, d_max_order);
            CVODE_SAMRAI_ERROR(ierr);
         }

         /*
          * Setup CVSpgmr function pointers.
          */
         CVSpilsPrecSetupFn precond_set = NULL;
         CVSpilsPrecSolveFn precond_solve = NULL;

         if (d_uses_preconditioner) {
            precond_set = CVODESolver::CVSpgmrPrecondSet;
            precond_solve = CVODESolver::CVSpgmrPrecondSolve;
            CVSpilsSetPreconditioner(d_cvode_mem, precond_set,
               precond_solve, this);
         }

         if (!(d_max_num_internal_steps < 0)) {
            ierr = CVodeSetMaxNumSteps(d_cvode_mem, d_max_num_internal_steps);
            CVODE_SAMRAI_ERROR(ierr);
         }

         if (!(d_max_num_warnings < 0)) {
            ierr = CVodeSetMaxHnilWarns(d_cvode_mem, d_max_num_warnings);
            CVODE_SAMRAI_ERROR(ierr);
         }

         if (!(d_init_step_size < 0)) {
            ierr = CVodeSetInitStep(d_cvode_mem, d_init_step_size);
            CVODE_SAMRAI_ERROR(ierr);
         }

         if (!(d_max_step_size < 0)) {
            ierr = CVodeSetMaxStep(d_cvode_mem, d_max_step_size);
            CVODE_SAMRAI_ERROR(ierr);
         }

         if (!(d_min_step_size < 0)) {
            ierr = CVodeSetMinStep(d_cvode_mem, d_min_step_size);
            CVODE_SAMRAI_ERROR(ierr);
         }
      }

   } // if no need to initialize CVODE, function does nothing

   d_CVODE_needs_initialization = false;
}

/*
 *************************************************************************
 *                                                                       *
 * Integrate system of ODEs to d_t_f.  If necessary, re-initialize       *
 * CVODE.                                                                *
 *                                                                       *
 *************************************************************************
 */
int CVODESolver::solve()
{

   int retval = CV_SUCCESS;

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

   retval = CVode(d_cvode_mem,
         d_user_t_f,
         d_solution_vector->getNVector(),
         &d_actual_t_f,
         d_stepping_method);

   return retval;

}

/*
 *************************************************************************
 *                                                                       *
 * Setting CVODE log file name and print flag for CVODE statistics.    *
 *                                                                       *
 *************************************************************************
 */

void CVODESolver::setLogFileData(
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
 *                                                                       *
 * Accessor functions for user-defined function and linear solver.       *
 *                                                                       *
 *************************************************************************
 */

void CVODESolver::setCVODEFunctions(
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

CVODEAbstractFunctions *CVODESolver::getCVODEFunctions() const
{
   return d_cvode_functions;
}
/*
 *************************************************************************
 *                                                                       *
 * Accessor functions for CVODE integration parameters.                  *
 *                                                                       *
 *************************************************************************
 */
void CVODESolver::setLinearMultistepMethod(
   int linear_multistep_method)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((linear_multistep_method == CV_ADAMS) ||
      (linear_multistep_method == CV_BDF));
#endif
   d_linear_multistep_method = linear_multistep_method;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setIterationType(
   int iteration_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((iteration_type == CV_FUNCTIONAL) ||
      (iteration_type == CV_NEWTON));
#endif
   d_iteration_type = iteration_type;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setToleranceType(
   int tolerance_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((tolerance_type == CV_SS) ||
      (tolerance_type == CV_SV));
#endif
   d_tolerance_type = tolerance_type;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setRelativeTolerance(
   double relative_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(relative_tolerance >= 0.0);
#endif

   d_relative_tolerance = relative_tolerance;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setAbsoluteTolerance(
   double absolute_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(absolute_tolerance >= 0.0);
#endif
   d_absolute_tolerance_scalar = absolute_tolerance;
   d_use_scalar_absolute_tolerance = true;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setAbsoluteTolerance(
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

void CVODESolver::setSteppingMethod(
   int stepping_method)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((stepping_method == CV_NORMAL) ||
      (stepping_method == CV_ONE_STEP));
#endif
   d_stepping_method = stepping_method;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setInitialValueOfIndependentVariable(
   double t_0)
{
   d_t_0 = t_0;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setFinalValueOfIndependentVariable(
   double t_f,
   bool cvode_needs_initialization)
{
   d_user_t_f = t_f;
   d_CVODE_needs_initialization = cvode_needs_initialization;
}

void CVODESolver::setInitialConditionVector(
   SundialsAbstractVector* ic_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(ic_vector == (SundialsAbstractVector *)NULL));
#endif
   d_ic_vector = ic_vector;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setMaximumLinearMultistepMethodOrder(
   int max_order)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_order >= 0);
#endif
   d_max_order = max_order;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setMaximumNumberOfInternalSteps(
   int max_num_internal_steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_num_internal_steps >= 0);
#endif

   d_max_num_internal_steps = max_num_internal_steps;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setMaximumNumberOfNilStepWarnings(
   int max_num_warnings)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_num_warnings >= 0);
#endif

   d_max_num_warnings = max_num_warnings;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setInitialStepSize(
   double init_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(init_step_size >= 0.0);
#endif
   d_init_step_size = init_step_size;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setMaximumAbsoluteStepSize(
   double max_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_step_size >= 0.0);
#endif
   d_max_step_size = max_step_size;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setMinimumAbsoluteStepSize(
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
 *                                                                       *
 * Accessor functions for CVSpgmr parameters.                            *
 *                                                                       *
 *************************************************************************
 */
void CVODESolver::setPreconditioningType(
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

void CVODESolver::setGramSchmidtType(
   int gs_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((gs_type == CLASSICAL_GS) ||
      (gs_type == MODIFIED_GS));
#endif
   d_gram_schmidt_type = gs_type;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setMaxKrylovDimension(
   int max_krylov_dim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_krylov_dim >= 0);
#endif
   d_max_krylov_dim = max_krylov_dim;
   d_CVODE_needs_initialization = true;
}

void CVODESolver::setCVSpgmrToleranceScaleFactor(
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
 *                                                                       *
 * Accessor functions for results of CVODE integration step.             *
 *                                                                       *
 *************************************************************************
 */
SundialsAbstractVector *CVODESolver::getSolutionVector() const
{
   return d_solution_vector;
}

int CVODESolver::getDkyVector(
   double t,
   int k,
   SundialsAbstractVector* dky) const
{
   int return_code;

   return_code = CVodeGetDky(d_cvode_mem, t, k, dky->getNVector());

   return return_code;
}

double CVODESolver::getActualFinalValueOfIndependentVariable() const
{
   return d_actual_t_f;
}

/*
 *************************************************************************
 *                                                                       *
 * Access methods for CVODE statistics.                                  *
 *                                                                       *
 *************************************************************************
 */

void CVODESolver::printStatistics(
   std::ostream& os) const
{
   printCVODEStatistics(os);
   printCVSpgmrStatistics(os);
}

void CVODESolver::printCVODEStatistics(
   std::ostream& os) const
{

   char buf[STAT_OUTPUT_BUFFER_SIZE];

   os << "\nCVODESolver: CVODE statistics... " << std::endl;

   sprintf(buf, "lenrw           = %5d     leniw            = %5d\n",
      getCVODEMemoryUsageForDoubles(),
      getCVODEMemoryUsageForIntegers());
   os << buf;
   sprintf(buf, "nst             = %5d     nfe              = %5d\n",
      getNumberOfInternalStepsTaken(),
      getNumberOfRHSFunctionCalls());
   os << buf;
   sprintf(buf, "nni             = %5d     nsetups          = %5d\n",
      getNumberOfNewtonIterations(),
      getNumberOfLinearSolverSetupCalls());
   os << buf;
   sprintf(buf, "netf            = %5d     ncfn             = %5d\n",
      getNumberOfLocalErrorTestFailures(),
      getNumberOfNonlinearConvergenceFailures());
   os << buf;
   sprintf(buf, "qu              = %5d     qcur             = %5d\n",
      getOrderUsedDuringLastInternalStep(),
      getOrderToBeUsedDuringNextInternalStep());
   os << buf;
   sprintf(buf, "\nhu              = %e      hcur             = %e\n",
      getStepSizeForLastInternalStep(),
      getStepSizeForNextInternalStep());
   os << buf;
   sprintf(buf, "tcur            = %e      tolsf            = %e\n",
      getCurrentInternalValueOfIndependentVariable(),
      getCVODESuggestedToleranceScalingFactor());
   os << buf;
}

int CVODESolver::getNumberOfInternalStepsTaken() const
{
   long int r;
   int ierr = CVodeGetNumSteps(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfRHSFunctionCalls() const
{
   long int r;
   int ierr = CVodeGetNumRhsEvals(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfLinearSolverSetupCalls() const
{
   long int r;
   int ierr = CVodeGetNumLinSolvSetups(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfNewtonIterations() const
{
   long int r;
   int ierr = CVodeGetNumNonlinSolvIters(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfNonlinearConvergenceFailures() const
{
   long int r;
   int ierr = CVodeGetNumNonlinSolvConvFails(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfLocalErrorTestFailures() const
{
   long int r;
   int ierr = CVodeGetNumErrTestFails(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getOrderUsedDuringLastInternalStep() const
{
   int r;
   int ierr = CVodeGetLastOrder(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getOrderToBeUsedDuringNextInternalStep() const
{
   int r;
   int ierr = CVodeGetCurrentOrder(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getCVODEMemoryUsageForDoubles() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r1);
}

int CVODESolver::getCVODEMemoryUsageForIntegers() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r2);
}

double CVODESolver::getStepSizeForLastInternalStep() const
{
   realtype r;
   int ierr = CVodeGetLastStep(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

double CVODESolver::getStepSizeForNextInternalStep() const
{
   realtype r;
   int ierr = CVodeGetCurrentStep(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

double CVODESolver::getCurrentInternalValueOfIndependentVariable() const
{
   realtype r;
   int ierr = CVodeGetCurrentStep(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

double CVODESolver::getCVODESuggestedToleranceScalingFactor() const
{
   realtype r;
   int ierr = CVodeGetTolScaleFactor(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return r;
}

/*
 *************************************************************************
 *                                                                       *
 * Access methods for CVSpgmr statistics.                                *
 *                                                                       *
 *************************************************************************
 */

void CVODESolver::printCVSpgmrStatistics(
   std::ostream& os) const
{
   if (d_iteration_type == CV_NEWTON) {

      os << "CVODESolver: CVSpgmr statistics... " << std::endl;

      os << "spgmr_lrw       = "
         << tbox::Utilities::intToString(getCVSpgmrMemoryUsageForDoubles(), 5)
         << "     spgmr_liw        = "
         << tbox::Utilities::intToString(getCVSpgmrMemoryUsageForIntegers(),
         5) << std::endl;

      os << "nli             = "
         << tbox::Utilities::intToString(getNumberOfLinearIterations(), 5)
         << "     ncfl             = "
         << tbox::Utilities::intToString(
         getNumberOfLinearConvergenceFailures(), 5) << std::endl;

      os << "npe             = "
         << tbox::Utilities::intToString(
         getNumberOfPreconditionerEvaluations(), 5)
         << "     nps              = "
         << tbox::Utilities::intToString(getNumberOfPrecondSolveCalls(),
         5) << std::endl;
   } else {

      os << "\nCVODESolver not set to use NEWTON iteration . . . \n"
         << std::endl;

   }
}

int CVODESolver::getNumberOfPreconditionerEvaluations() const
{
   long int r;
   int ierr = CVSpilsGetNumPrecEvals(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfLinearIterations() const
{
   long int r;
   int ierr = CVSpilsGetNumLinIters(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfPrecondSolveCalls() const
{
   long int r;
   int ierr = CVSpilsGetNumPrecSolves(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getNumberOfLinearConvergenceFailures() const
{
   long int r;
   int ierr = CVSpilsGetNumConvFails(d_cvode_mem, &r);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r);
}

int CVODESolver::getCVSpgmrMemoryUsageForDoubles() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r1);
}

int CVODESolver::getCVSpgmrMemoryUsageForIntegers() const
{
   long int r1;
   long int r2;
   int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
   CVODE_SAMRAI_ERROR(ierr);
   return static_cast<int>(r2);
}

/*
 *************************************************************************
 *                                                                       *
 * Print CVODESolver object data to given output stream.            *
 *                                                                       *
 *************************************************************************
 */
void CVODESolver::printClassData(
   std::ostream& os) const
{
   os << "\nCVODESolver object data members..." << std::endl;
   os << "Object name = "
      << d_object_name << std::endl;

   os << "this = " << (CVODESolver *)this << std::endl;
   os << "d_solution_vector = "
      << (SundialsAbstractVector *)d_solution_vector << std::endl;

   os << "d_CVODE_functions = "
      << (CVODEAbstractFunctions *)d_cvode_functions << std::endl;

   os << "&d_cvode_mem = " << d_cvode_mem << std::endl;
   os << "d_cvode_log_file = " << (FILE *)d_cvode_log_file << std::endl;
   os << "d_cvode_log_file_name = " << d_cvode_log_file_name << std::endl;

   os << std::endl;
   os << "CVODE parameters..." << std::endl;
   os << "d_t_0 = "
      << d_t_0 << std::endl;
   os << "d_ic_vector = "
      << (SundialsAbstractVector *)d_ic_vector << std::endl;

   os << "d_linear_multistep_method = "
      << d_linear_multistep_method << std::endl;
   os << "d_iteration_type = "
      << d_iteration_type << std::endl;
   os << "d_tolerance_type = "
      << d_tolerance_type << std::endl;
   os << "d_relative_tolerance = "
      << d_relative_tolerance << std::endl;
   os << "d_use_scalar_absolute_tolerance = ";
   if (d_use_scalar_absolute_tolerance) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   }
   os << "d_absolute_tolerance_scalar = "
      << d_absolute_tolerance_scalar << std::endl;
   os << "d_absolute_tolerance_vector= " << std::endl;
   d_absolute_tolerance_vector->printVector();

   os << "Optional CVODE inputs (see CVODE docs for details):"
      << std::endl;

   os << "maximum linear multistep method order = "
      << d_max_order << std::endl;
   os << "maximum number of internal steps = "
      << d_max_num_internal_steps << std::endl;
   os << "maximum number of nil internal step warnings = "
      << d_max_num_warnings << std::endl;

   os << "initial step size = "
      << d_init_step_size << std::endl;
   os << "maximum absolute value of step size = "
      << d_max_step_size << std::endl;
   os << "minimum absolute value of step size = "
      << d_min_step_size << std::endl;
   os << "last step size = "
      << getStepSizeForLastInternalStep() << std::endl;
   os << "...end of CVODE parameters\n" << std::endl;

   os << std::endl;
   os << "CVSpgmr parameters..." << std::endl;
   os << "d_precondition_type = "
      << d_precondition_type << std::endl;
   os << "d_gram_schmidt_type = "
      << d_gram_schmidt_type << std::endl;
   os << "d_max_krylov_dim = "
      << d_max_krylov_dim << std::endl;
   os << "d_tol_scale_factor = "
      << d_tol_scale_factor << std::endl;
   os << "...end of CVSpgmr parameters\n" << std::endl;

   os << "d_CVODE_needs_initialization = ";
   if (d_CVODE_needs_initialization) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   }

   os << "...end of CVODESolver object data members\n" << std::endl;
}

}
}

#endif
