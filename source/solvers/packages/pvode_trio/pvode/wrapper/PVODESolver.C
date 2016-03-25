//
// File:        PVODESolver.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: C++ Wrapper class for PVODE solver package 
//

#include "PVODESolver.h"

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#ifdef HAVE_PVODE

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "tbox/IEEE.h"

// PVODE includes
#ifndef included_cvode_h
#define included_cvode_h
extern "C" {
#include "cvode.h"
}
#endif

#ifndef included_cvspgmr_h
#define included_cvspgmr_h
extern "C" {
#include "cvspgmr.h"
}
#endif

#ifndef STAT_OUTPUT_BUFFER_SIZE
#define STAT_OUTPUT_BUFFER_SIZE 256
#endif


namespace SAMRAI {
   namespace solv {


/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with PVODE package.      *
* See header file for PVODEAbstractFunctions for more information. *
*                                                                       *
*************************************************************************
*/
void PVODESolver::PVODERHSFuncEval(int neq,
                                        double t,
                                        PVodeTrioAbstractVector* y,
                                        PVodeTrioAbstractVector* y_dot,
                                        void* my_solver)
{
   NULL_USE(neq);
   ((PVODESolver*)my_solver)->getPVODEFunctions()->
                                    evaluateRHSFunction(t, y, y_dot);
}

/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with CVSpgmr package.	*
*                                                                       *
*************************************************************************
*/

int PVODESolver::CVSpgmrPrecondSet(int neq,
                      double t,
                      PVodeTrioAbstractVector* y,
                      PVodeTrioAbstractVector* fy,
                      int jok,
                      int *jcurPtr,
                      double gamma,
                      PVodeTrioAbstractVector* ewt,
                      double h,
                      double mach_roundoff,
                      long int *nfePtr,
                      void *my_solver,
                      PVodeTrioAbstractVector* vtemp1,
                      PVodeTrioAbstractVector* vtemp2,
                      PVodeTrioAbstractVector* vtemp3)
{

   int success;
   success = ((PVODESolver*)my_solver)->getPVODEFunctions()->
             CVSpgmrPrecondSet(neq,
                               t,
                               y,
                               fy,
                               jok,
                               jcurPtr,
                               gamma,
                               ewt,
                               h,
                               mach_roundoff,
                               nfePtr,
                               my_solver,
                               vtemp1,
                               vtemp2,
                               vtemp3);

   return (success);
}

int PVODESolver::CVSpgmrPrecondSolve(int neq,
			double t,
			PVodeTrioAbstractVector* y,
			PVodeTrioAbstractVector* fy,
			PVodeTrioAbstractVector* vtemp,
			double gamma,
			PVodeTrioAbstractVector* ewt,
			double delta,
			long int *nfePtr,
			PVodeTrioAbstractVector* r,
			int lr,
			void *my_solver,
			PVodeTrioAbstractVector* z)
{

   int success;
   success = ((PVODESolver*)my_solver)->getPVODEFunctions()->
             CVSpgmrPrecondSolve(neq,
                                 t,
                                 y,
                                 fy,
                                 vtemp,
                                 gamma,
                                 ewt,
                                 delta,
                                 nfePtr,
                                 r,
                                 lr,
                                 my_solver,
                                 z);

   return (success);
}


/*
*************************************************************************
*                                                                       *
* PVODESolver constructor and destructor.                         *
*                                                                       *
*************************************************************************
*/
PVODESolver::PVODESolver(
   const string& object_name,
   PVODEAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!(my_functions == (PVODEAbstractFunctions*)NULL));
#endif

   d_object_name         = object_name;
   d_pvode_functions      = my_functions;
   d_uses_preconditioner = uses_preconditioner;

   d_solution_vector = (PVodeTrioAbstractVector*)NULL;

   /*
    * Set default parameters to safe values or to PVODE/CVSpgmr defaults.
    */

   /*
    * PVODE memory record and log file.
    */
   d_pvode_mem           = NULL;
   d_pvode_log_file      = NULL;
   d_pvode_log_file_name = "pvode.log";
 
   /*
    * ODE parameters.
    */ 
   d_neq = 0;
   d_t_0 = 0.0;
   d_user_t_f = 0.0;
   d_actual_t_f = 0.0;
   d_ic_vector = ((PVodeTrioAbstractVector*)NULL);
 
   /*
    * ODE integration parameters.
    */
   setLinearMultistepMethod(BDF);
   setIterationType(FUNCTIONAL);
   setToleranceType(SS);
   setRelativeTolerance(0.0);
   setAbsoluteTolerance(0.0);
   d_absolute_tolerance_vector = (PVodeTrioAbstractVector*)NULL;
   setSteppingMethod(NORMAL);

   /*
    * Optional PVODE parameters.
    */
   d_use_optional_data = false;
   for (int i = 0; i < OPT_SIZE; i++) {
      d_int_optional_data[i] = 0;
      d_real_optional_data[i] = 0.0;
   }

   /* 
    * CVSpgmr parameters.  
    * 
    * Note that when the maximum krylov dimension and CVSpgmr
    * tolerance scale factor are set to 0, CVSpgmr uses its
    * internal default values.  These are described in the header for 
    * this class.
    */
   setPreconditioningType(NONE);
   setGramSchmidtType(MODIFIED_GS);
   setMaxKrylovDimension(0);
   setCVSpgmrToleranceScaleFactor(0);

   d_PVODE_needs_initialization = true;
}

PVODESolver::~PVODESolver()
{
   if (d_pvode_log_file) fclose(d_pvode_log_file);
   if (d_pvode_mem) CVodeFree(d_pvode_mem);
}

/*
*************************************************************************
*                                                                       *
* Functions to initialize linear solver and reset PVODE structure.      *
*                                                                       *
*************************************************************************
*/

void PVODESolver::initialize(PVodeTrioAbstractVector* solution)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(solution == (PVodeTrioAbstractVector*)NULL));
   assert(d_solution_vector == (PVodeTrioAbstractVector*)NULL);
#endif
   d_solution_vector = solution;
   initializePVODE();
}

void PVODESolver::initializePVODE() 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(d_solution_vector == (PVodeTrioAbstractVector*)NULL));
#endif
   if (d_PVODE_needs_initialization) {
  
      /*
       * Set PVODE log file.
       */ 
      if (d_pvode_log_file) {
         fclose(d_pvode_log_file);
      }
      d_pvode_log_file = fopen(d_pvode_log_file_name.c_str(), "w");

      /*
       * Make sure that the number of equations is positive.
       */
      if (d_neq <= 0 && d_pvode_log_file) {
         fprintf(d_pvode_log_file, 
                 "%s: The number of equations is not positive.",
                 d_object_name.c_str());
      }

      /*
       * Make sure that either the relative tolerance or the
       * absolute tolerance has been set to a nonzero value.
       */
      bool tolerance_error = false;
      if (d_use_scalar_absolute_tolerance) {
         if ( (d_relative_tolerance == 0.0) &&
              (d_absolute_tolerance_scalar == 0.0) ) { 
            tolerance_error = true;
         }
      } else {
         if ( (d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_vector->maxNorm() == 0.0) ) {
            tolerance_error = true;
         }
      }

      if (tolerance_error && d_pvode_log_file) {
         fprintf(d_pvode_log_file, 
                 "%s: Both relative and absolute tolerance have value 0.0", 
                 d_object_name.c_str());
      }

      /*
       * PVODE function pointer.
       */
      RhsFn RHSFunc = PVODESolver::PVODERHSFuncEval;

      /*
       * Set optional data array pointers.
       */
      int optional_data_flag = (d_use_optional_data ? 1 : 0);
      double* real_opt_data_array = 
         (d_use_optional_data ? d_real_optional_data : NULL);
      long int* int_opt_data_array = 
         (d_use_optional_data ? d_int_optional_data : NULL);

      /*
       * Set (void*) absolute tolerance pointer to be passed
       * as an argument to CVodeMalloc().
       */
      void* absolute_tolerance;
      if (d_use_scalar_absolute_tolerance) {
         absolute_tolerance = (void*) &d_absolute_tolerance_scalar;
      } else {
         absolute_tolerance = (void*) d_absolute_tolerance_vector;
      }

      /*
       * Free previously allocated CVode memory.  Note that the
       * CVReInit() function is not used since the d_neq variable
       * might have been changed from the previous initializePVODE()
       * call.
       */
      if (d_pvode_mem) CVodeFree(d_pvode_mem);

      /*
       * Allocate main memory for PVODE package.
       */
      d_pvode_mem = (CVodeMem)CVodeMalloc(d_neq,
                                          RHSFunc,
                                          d_t_0,
                                          d_ic_vector,
                                          d_linear_multistep_method,
                                          d_iteration_type,
                                          d_tolerance_type,
                                          &d_relative_tolerance,
                                          absolute_tolerance,
                                          (void*) this,
                                          d_pvode_log_file,
                                          optional_data_flag,
                                          int_opt_data_array,
                                          real_opt_data_array,
                                          (void*) d_solution_vector);

      /*
       * If the iteration type is set to NEWTON, then initialize
       * the CVSpgmr linear solver.
       */
      if (d_iteration_type == NEWTON) {

         /*
          * Setup CVSpgmr function pointers.
          */
         CVSpgmrPrecondFn precond_set   = NULL;
         CVSpgmrPSolveFn  precond_solve = NULL;

        if (d_uses_preconditioner) {
            precond_set	  = PVODESolver::CVSpgmrPrecondSet;
            precond_solve = PVODESolver::CVSpgmrPrecondSolve;
        }
 
        CVSpgmr(d_pvode_mem,
                d_precondition_type,
                d_gram_schmidt_type,
                d_max_krylov_dim,
                d_tol_scale_factor,
                precond_set,
                precond_solve,
                (void*) this);
      }

   } // if no need to initialize PVODE, function does nothing

   d_PVODE_needs_initialization = false;
}


/*
*************************************************************************
*                                                                       *
* Integrate system of ODEs to d_t_f.  If necessary, re-initialize       *
* PVODE.                                                                *
*                                                                       *
*************************************************************************
*/
int PVODESolver::solve()
{

   int retval = SUCCESS;

   initializePVODE();

   /* 
    * Check to make sure that user specified final value for t
    * is greater than initial value for t.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( d_user_t_f > d_t_0);
#endif

   /*
    * See cvode.h header file for definition of return types.
    */

   retval = CVode((void*) d_pvode_mem,
                     d_user_t_f,
                     d_solution_vector,
                     &d_actual_t_f,
                     d_stepping_method);                  

   return( retval );

}

/*
*************************************************************************
*                                                                       *
* Setting PVODE log file name and print flag for PVODE statistics.    *
*                                                                       *
*************************************************************************
*/

void PVODESolver::setLogFileData(
   const string& log_fname)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!log_fname.empty());
#endif
   if ( !(log_fname == d_pvode_log_file_name) ) {
      d_pvode_log_file_name = log_fname;
      d_PVODE_needs_initialization = true;
   }
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for user-defined function and linear solver.       *
*                                                                       *
*************************************************************************
*/

void PVODESolver::setPVODEFunctions(
   PVODEAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(my_functions == (PVODEAbstractFunctions*)NULL));
#endif

   d_pvode_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_PVODE_needs_initialization = true;
}

PVODEAbstractFunctions* PVODESolver::getPVODEFunctions() const
{
   return(d_pvode_functions);
}

void PVODESolver::setNumberOfEquations(int neq) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( neq > 0 );
#endif
   d_neq = neq;
   d_PVODE_needs_initialization = true;
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for PVODE integration parameters.                  *
*                                                                       *
*************************************************************************
*/
void PVODESolver::setLinearMultistepMethod(int linear_multistep_method)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (linear_multistep_method == ADAMS) || 
           (linear_multistep_method == BDF) );
#endif
   d_linear_multistep_method = linear_multistep_method;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setIterationType(int iteration_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (iteration_type == FUNCTIONAL) || 
           (iteration_type == NEWTON) );
#endif
   d_iteration_type = iteration_type;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setToleranceType(int tolerance_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (tolerance_type == SS) || 
           (tolerance_type == SV) );
#endif
   d_tolerance_type = tolerance_type;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setRelativeTolerance(double relative_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(relative_tolerance >= 0.0);
#endif

   d_relative_tolerance = relative_tolerance;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setAbsoluteTolerance(double absolute_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(absolute_tolerance >= 0.0);
#endif
   d_absolute_tolerance_scalar = absolute_tolerance;
   d_use_scalar_absolute_tolerance = true;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setAbsoluteTolerance(
   PVodeTrioAbstractVector* absolute_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( !(absolute_tolerance == (PVodeTrioAbstractVector*)NULL) );
   assert( absolute_tolerance->vecMin() >= 0.0 ); 
#endif
   d_absolute_tolerance_vector = absolute_tolerance;
   d_use_scalar_absolute_tolerance = false;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setSteppingMethod(int stepping_method)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (stepping_method == NORMAL) || 
           (stepping_method == ONE_STEP) );
#endif
   d_stepping_method = stepping_method;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setInitialValueOfIndependentVariable(double t_0)
{
   d_t_0 = t_0;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setFinalValueOfIndependentVariable(double t_f,
   bool pvode_needs_initialization)
{
   d_user_t_f = t_f;
   d_PVODE_needs_initialization = pvode_needs_initialization;
}

void PVODESolver::setInitialConditionVector(
   PVodeTrioAbstractVector* ic_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( !(ic_vector == (PVodeTrioAbstractVector*)NULL) );
#endif
   d_ic_vector = ic_vector;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setMaximumLinearMultistepMethodOrder(
   int max_order)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( max_order >= 0);
#endif

   d_int_optional_data[MAXORD] = max_order;
   d_use_optional_data = true;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setMaximumNumberOfInternalSteps(
   int max_num_internal_steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( max_num_internal_steps >= 0 );
#endif

   d_int_optional_data[MXSTEP] = max_num_internal_steps;
   d_use_optional_data = true;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setMaximumNumberOfNilStepWarnings(
   int max_num_warnings)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( max_num_warnings >= 0 );
#endif

   d_int_optional_data[MXHNIL] = max_num_warnings;
   d_use_optional_data = true;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setInitialStepSize(double init_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( init_step_size >= 0.0 );
#endif
   d_real_optional_data[H0] = init_step_size;
   d_use_optional_data = true;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setMaximumAbsoluteStepSize(double max_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( max_step_size >= 0.0 );
#endif
   d_real_optional_data[HMAX] = max_step_size;
   d_use_optional_data = true;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setMinimumAbsoluteStepSize(double min_step_size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( min_step_size >= 0.0 );
#endif
   d_real_optional_data[HMIN] = min_step_size;
   d_use_optional_data = true;
   d_PVODE_needs_initialization = true;
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for CVSpgmr parameters.                            *
*                                                                       *
*************************************************************************
*/
void PVODESolver::setPreconditioningType(int precondition_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (precondition_type == NONE) ||
           (precondition_type == LEFT) ||
           (precondition_type == RIGHT) ||
           (precondition_type == BOTH) );
#endif
   d_precondition_type = precondition_type;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setGramSchmidtType(int gs_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( (gs_type == CLASSICAL_GS) ||
           (gs_type == MODIFIED_GS) );
#endif
   d_gram_schmidt_type = gs_type;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setMaxKrylovDimension(int max_krylov_dim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( max_krylov_dim >= 0 );
#endif
   d_max_krylov_dim = max_krylov_dim;
   d_PVODE_needs_initialization = true;
}

void PVODESolver::setCVSpgmrToleranceScaleFactor(double tol_scale_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( tol_scale_factor >= 0 );
#endif
   d_tol_scale_factor = tol_scale_factor;
   d_PVODE_needs_initialization = true;
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for results of PVODE integration step.             *
*                                                                       *
*************************************************************************
*/
PVodeTrioAbstractVector* PVODESolver::getSolutionVector() const
{
   return (d_solution_vector);
}

int PVODESolver::getDkyVector(double t, int k,
   PVodeTrioAbstractVector* dky) const
{
   int return_code;

   return_code = CVodeDky((void*)d_pvode_mem, t, k, dky);

   return (return_code);
}

double PVODESolver::getActualFinalValueOfIndependentVariable() const
{
   return (d_actual_t_f);
}


/*
*************************************************************************
*                                                                       *
* Access methods for PVODE statistics.                                  *
*                                                                       *
*************************************************************************
*/

void PVODESolver::turnOnPVODEStatisticsCollection()
{
   d_use_optional_data = true;
}

void PVODESolver::printStatistics(ostream& os) const
{
   printPVODEStatistics(os);
   printCVSpgmrStatistics(os);
}

void PVODESolver::printPVODEStatistics(ostream& os) const
{

   if (d_use_optional_data) {
      char buf[STAT_OUTPUT_BUFFER_SIZE];

      os << "\nPVODESolver: PVODE statistics... " << endl;

      sprintf(buf, "lenrw           = %5d     leniw            = %5d\n", 
         getPVODEMemoryUsageForDoubles(), 
         getPVODEMemoryUsageForIntegers());
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
         getPVODESuggestedToleranceScalingFactor());
      os << buf;

   } else {

      os << "\nPVODESolver not set to collect statistics . . . \n" << endl;

   }
}

int PVODESolver::getNumberOfInternalStepsTaken() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[NST];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfRHSFunctionCalls() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[NFE];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfLinearSolverSetupCalls() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[NSETUPS];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfNewtonIterations() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[NNI];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfNonlinearConvergenceFailures() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[NCFN];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfLocalErrorTestFailures() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[NETF];
   } else {
      return (-1);
   } 
}

int PVODESolver::getOrderUsedDuringLastInternalStep() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[QU];
   } else {
      return (-1);
   } 
}

int PVODESolver::getOrderToBeUsedDuringNextInternalStep() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[QCUR];
   } else {
      return (-1);
   } 
}

int PVODESolver::getPVODEMemoryUsageForDoubles() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[LENRW];
   } else {
      return (-1);
   } 
}

int PVODESolver::getPVODEMemoryUsageForIntegers() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[LENIW];
   } else {
      return (-1);
   } 
}

double PVODESolver::getStepSizeForLastInternalStep() const
{
   if (d_use_optional_data) {
      return d_real_optional_data[HU];
   } else {
      return (-1);
   } 
}

double PVODESolver::getStepSizeForNextInternalStep() const
{
   if (d_use_optional_data) {
      return d_real_optional_data[HCUR];
   } else {
      return (-1);
   } 
}

double PVODESolver::getCurrentInternalValueOfIndependentVariable() const
{
   if (d_use_optional_data) {
      return d_real_optional_data[TCUR];
   } else {
      return (tbox::IEEE::getDBL_MIN());
   } 
}

double PVODESolver::getPVODESuggestedToleranceScalingFactor() const
{
   if (d_use_optional_data) {
      return d_real_optional_data[TOLSF];
   } else {
      return (-1);
   } 
}

/*
*************************************************************************
*                                                                       *
* Access methods for CVSpgmr statistics.                                *
*                                                                       *
*************************************************************************
*/

void PVODESolver::printCVSpgmrStatistics(ostream& os) const
{
   if (d_iteration_type == NEWTON) {
      if (d_use_optional_data) {
         char buf[STAT_OUTPUT_BUFFER_SIZE];

         os << "PVODESolver: CVSpgmr statistics... " << endl;

         sprintf(buf, "spgmr_lrw       = %5d     spgmr_liw        = %5d\n", 
            getCVSpgmrMemoryUsageForDoubles(),
            getCVSpgmrMemoryUsageForIntegers());
         os << buf;
         sprintf(buf, "nli             = %5d     ncfl             = %5d\n", 
            getNumberOfLinearIterations(),
            getNumberOfLinearConvergenceFailures());
         os << buf;
         sprintf(buf, "npe             = %5d     nps              = %5d\n", 
            getNumberOfPreconditionerEvaluations(),
            getNumberOfPrecondSolveCalls());
         os << buf;
      } else { 

         os << "\nPVODESolver not set to collect statistics . . . \n" 
            << endl;
      }
   } else {
 
      os << "\nPVODESolver not set to use NEWTON iteration . . . \n"
         << endl;

   }
}
   
int PVODESolver::getNumberOfPreconditionerEvaluations() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[SPGMR_NPE];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfLinearIterations() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[SPGMR_NLI];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfPrecondSolveCalls() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[SPGMR_NPS];
   } else {
      return (-1);
   } 
}

int PVODESolver::getNumberOfLinearConvergenceFailures() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[SPGMR_NCFL];
   } else {
      return (-1);
   } 
}

int PVODESolver::getCVSpgmrMemoryUsageForDoubles() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[SPGMR_LRW];
   } else {
      return (-1);
   } 
}

int PVODESolver::getCVSpgmrMemoryUsageForIntegers() const
{
   if (d_use_optional_data) {
      return d_int_optional_data[SPGMR_LIW];
   } else {
      return (-1);
   } 
}

/*
*************************************************************************
*                                                                       *
* Print PVODESolver object data to given output stream.            *
*                                                                       *
*************************************************************************
*/
void PVODESolver::printClassData(ostream& os) const
{
   os << "\nPVODESolver object data members..." << endl;
   os << "Object name = "
      << d_object_name << endl;

   os << "this = " << (PVODESolver*)this << endl;
   os << "d_solution_vector = " 
      << (PVodeTrioAbstractVector*)d_solution_vector << endl;

   os << "d_PVODE_functions = " 
      << (PVODEAbstractFunctions*)d_pvode_functions << endl;

   os << "&d_pvode_mem = " << (CVodeMem*) &d_pvode_mem << endl;
   os << "d_pvode_log_file = " << (FILE*) d_pvode_log_file << endl;
   os << "d_pvode_log_file_name = " << d_pvode_log_file_name << endl;

   os << endl;
   os << "PVODE parameters..." << endl;
   os << "d_neq = "
      << d_neq << endl;
   os << "d_t_0 = "
      << d_t_0 << endl;
   os << "d_ic_vector = "
      << (PVodeTrioAbstractVector*)d_ic_vector << endl;

   os << "d_linear_multistep_method = "
      << d_linear_multistep_method << endl;
   os << "d_iteration_type = "
      << d_iteration_type << endl;
   os << "d_tolerance_type = "
      << d_tolerance_type << endl;
   os << "d_relative_tolerance = "
      << d_relative_tolerance << endl;
   os << "d_use_scalar_absolute_tolerance = ";
   if (d_use_scalar_absolute_tolerance) {
      os << "true" << endl;
   } else {
      os << "false" << endl;
   } 
   os << "d_absolute_tolerance_scalar = "
      << d_absolute_tolerance_scalar << endl;
   os << "d_absolute_tolerance_vector= " << endl;
   d_absolute_tolerance_vector->printVector(); 

   os << "Optional PVODE inputs (see PVODE docs for details):" 
      << endl;

   os << "d_use_optional_data = ";
   if (d_use_optional_data) {
      os << "true" << endl;
   } else {
      os << "false" << endl;
   } 

   os << "maximum linear multistep method order = "
      << d_int_optional_data[MAXORD] << endl;
   os << "maximum number of internal steps = "
      << d_int_optional_data[MXSTEP] << endl;
   os << "maximum number of nil internal step warnings = "
      << d_int_optional_data[MXHNIL] << endl;

   os << "initial step size = "
      << d_real_optional_data[H0] << endl;
   os << "maximum absolute value of step size = "
      << d_real_optional_data[HMAX] << endl;
   os << "minimum absolute value of step size = "
      << d_real_optional_data[HMIN] << endl;
   os << "last step size = "
      << d_real_optional_data[HU] << endl;
   os << "...end of PVODE parameters\n" << endl;

   os << endl;
   os << "CVSpgmr parameters..." << endl;
   os << "d_precondition_type = "
	<< d_precondition_type << endl;
   os << "d_gram_schmidt_type = "
	<< d_gram_schmidt_type << endl;
   os << "d_max_krylov_dim = "
	<< d_max_krylov_dim << endl;
   os << "d_tol_scale_factor = "
	<< d_tol_scale_factor << endl;
   os << "...end of CVSpgmr parameters\n" << endl;

   os << "d_PVODE_needs_initialization = ";
   if (d_PVODE_needs_initialization) {
      os << "true" << endl;
   } else {
      os << "false" << endl;
   } 

   os << "...end of PVODESolver object data members\n" << endl;
}


}
}

#endif
