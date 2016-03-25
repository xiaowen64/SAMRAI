//
// File:        KINSOLSolver.C
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: C++ Wrapper class for KINSOL solver package 
//

#include "KINSOLSolver.h"

#include "tbox/Utilities.h"

#ifdef HAVE_KINSOL

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace solv {

/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with KINSOL package.     *
* See header file for KINSOLAbstractFunctions for more information.*
*                                                                       *
*************************************************************************
*/
void KINSOLSolver::KINSOLFuncEval(int neq,
                                       PVodeTrioAbstractVector* soln,
                                       PVodeTrioAbstractVector* fval,
                                       void* my_solver)
{
   NULL_USE(neq);
   ((KINSOLSolver*)my_solver)->getKINSOLFunctions()->
                                    evaluateNonlinearFunction(soln, fval);
}

int KINSOLSolver::KINSOLPrecondSet(int neq,
                                        PVodeTrioAbstractVector* soln,
                                        PVodeTrioAbstractVector* soln_scale,
                                        PVodeTrioAbstractVector* fval,
                                        PVodeTrioAbstractVector* fval_scale,
                                        PVodeTrioAbstractVector* vtemp1,
                                        PVodeTrioAbstractVector* vtemp2,
                                        SysFn sys_func,
                                        double mach_roundoff,
                                        long int* nfePtr,
                                        void* my_solver)
{
   ((KINSOLSolver*)my_solver)->initializeKINSOL();

   int success = 0;

   NULL_USE(neq);
   NULL_USE(sys_func);
   int num_feval = 0;
   success = ((KINSOLSolver*)my_solver)->getKINSOLFunctions()->
                                                 precondSetup(soln,
                                                              soln_scale,
                                                              fval,
                                                              fval_scale,
                                                              vtemp1,
                                                              vtemp2, 
                                                              mach_roundoff,
                                                              num_feval);
   if (num_feval > 0) {
      *nfePtr += num_feval;
   }

   return(success);
}

int
KINSOLSolver::KINSOLPrecondSolve(int neq,
                                      PVodeTrioAbstractVector* soln,
                                      PVodeTrioAbstractVector* soln_scale,
                                      PVodeTrioAbstractVector* fval,
                                      PVodeTrioAbstractVector* fval_scale,
                                      PVodeTrioAbstractVector* rhs,
                                      PVodeTrioAbstractVector* vtemp,
                                      SysFn sys_func,
                                      double mach_roundoff,
                                      long int* nfePtr,
                                      void* my_solver)
{
   int success = 0;

   NULL_USE(neq);
   NULL_USE(sys_func);
   int num_feval = 0;
   success = ((KINSOLSolver*)my_solver)->getKINSOLFunctions()->
                                              precondSolve(soln,
                                                           soln_scale,
                                                           fval,
                                                           fval_scale,
                                                           rhs,
                                                           vtemp,
                                                           mach_roundoff,
                                                           num_feval);
   if (num_feval > 0) {
      *nfePtr += num_feval;
   }

   return(success);
}

int
KINSOLSolver::KINSOLJacobianTimesVector(void* my_solver,
                                             PVodeTrioAbstractVector* vector,
                                             PVodeTrioAbstractVector* prod,
                                             int* new_soln,
                                             PVodeTrioAbstractVector* soln)
{
   int success = 0;

   bool soln_changed = true; 
   if (*new_soln == 0) {
      soln_changed = false; 
   }

   success = ((KINSOLSolver*)my_solver)->
             getKINSOLFunctions()->jacobianTimesVector(vector,
                                                       prod,
                                                       soln_changed,
                                                       soln);

   return(success);
}


/*
*************************************************************************
*                                                                       *
* KINSOLSolver constructor and destructor.                         *
*                                                                       *
*************************************************************************
*/
KINSOLSolver::KINSOLSolver(
   const string& object_name,
   KINSOLAbstractFunctions* my_functions,
   const int uses_preconditioner,
   const int uses_jac_times_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!(my_functions == (KINSOLAbstractFunctions*)NULL));
#endif

   d_object_name = object_name;
   d_KINSOL_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_uses_jac_times_vector = uses_jac_times_vector;

   d_solution_vector = (PVodeTrioAbstractVector*)NULL;
   d_KINSOL_needs_initialization = true;

   /*
    * Default parameters to safe values or to KINSOL defaults.
    */

   d_kin_mem              = NULL;
   d_kinsol_log_file      = NULL;

   d_soln_scale = (PVodeTrioAbstractVector*)NULL;
   d_my_soln_scale_vector = false;
   d_fval_scale = (PVodeTrioAbstractVector*)NULL;
   d_my_fval_scale_vector = false;

   d_constraints = (PVodeTrioAbstractVector*)NULL;

   d_krylov_dimension   = KINSPGMR_MAXL;
   d_max_restarts       = 0;
   d_max_solves_no_set  = KINSPGMR_MSBPRE;

   d_global_strategy    = INEXACT_NEWTON;
   d_residual_tol       = 0.0;
   d_step_tol           = 0.0;

   d_int_optional_input[ETACHOICE] = ETACONSTANT;

   /*
    * Optional input parameters for KINSOL.
    */
   for (int i = 0; i < OPT_SIZE; i++) {
      d_int_optional_input[i] = 0;
      d_real_optional_input[i] = 0.0;
   }

}

KINSOLSolver::~KINSOLSolver()
{
   if (d_my_soln_scale_vector && d_my_fval_scale_vector && d_soln_scale) {
      d_soln_scale->freeVector();
      d_my_soln_scale_vector = false;
      d_my_fval_scale_vector = false;
   }
   if (d_my_soln_scale_vector && d_soln_scale) d_soln_scale->freeVector();
   if (d_my_fval_scale_vector && d_fval_scale) d_fval_scale->freeVector();

   if ( d_kinsol_log_file ) {
     fclose(d_kinsol_log_file);
   }
   if (d_kin_mem) KINFree(d_kin_mem);
}

/*
*************************************************************************
*                                                                       *
* Functions to initialize nonlinear solver and reset KINSOL structure.  *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::initialize(PVodeTrioAbstractVector* solution)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(solution == (PVodeTrioAbstractVector*)NULL));
   assert(d_solution_vector == (PVodeTrioAbstractVector*)NULL);
#endif
   d_solution_vector = solution;
   initializeKINSOL();
}

void KINSOLSolver::initializeKINSOL() 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(d_solution_vector == (PVodeTrioAbstractVector*)NULL));
#endif
   if (d_KINSOL_needs_initialization) {
   
      if (d_kinsol_log_file) {
         fclose(d_kinsol_log_file);
      }

      d_kinsol_log_file = fopen(d_kinsol_log_file_name.c_str(), "w");

      /*
       * KINSOL function pointers.
       */

      KINSpgmrPrecondFn      precond_set   = NULL;
      KINSpgmrPrecondSolveFn precond_solve = NULL;
      KINSpgmruserAtimesFn   jac_times_vec = NULL;

      if (d_uses_preconditioner) {
         precond_set = KINSOLSolver::KINSOLPrecondSet;
         precond_solve = KINSOLSolver::KINSOLPrecondSolve;
      } 

      if (d_uses_jac_times_vector) {
         jac_times_vec = KINSOLSolver::KINSOLJacobianTimesVector;
      }

      if (d_kin_mem) KINFree(d_kin_mem);

      /*
       * Allocate main memory for KINSOL package.
       */
      d_neq = d_krylov_dimension * (1 + d_max_restarts);
      d_kin_mem = (KINMem)KINMalloc( d_neq,
                                     d_kinsol_log_file,
                                     (void*)d_solution_vector );

      /*
       * Initialize KINSOL memory record.
       */
      KINSpgmr( (void*)d_kin_mem,
                d_krylov_dimension,
                d_max_restarts,
                d_max_solves_no_set,
                precond_set,
                precond_solve,
                jac_times_vec,
                (void*)this ); 

   } // if no need to initialize KINSOL, function does nothing

   d_KINSOL_needs_initialization = false;
}


/*
*************************************************************************
*                                                                       *
* Solve nonlinear system; re-initialize KINSOL solver, if necessary.    *
*                                                                       *
*************************************************************************
*/
int KINSOLSolver::solve()
{

   int retval = KINSOL_SUCCESS;

   initializeKINSOL();

   /*
    * If scaling vectors are not provided, we make defaults here.
    */
   if (!d_soln_scale) {
      d_soln_scale = d_solution_vector->makeNewVector();
      d_soln_scale->setToScalar(1.0);
      d_my_soln_scale_vector = true;

      if (!d_fval_scale) {
         d_fval_scale = d_soln_scale;
         d_my_fval_scale_vector = true;
      }
   }
  
   if (!d_fval_scale) {
      d_fval_scale = d_solution_vector->makeNewVector();
      d_fval_scale->setToScalar(1.0);
      d_my_fval_scale_vector = true;
   }

   /*
    * See kinsol.h header file for definition of return types.
    */

   retval = KINSol( (void*)d_kin_mem,
                     d_neq,
                     d_solution_vector,
                     KINSOLSolver::KINSOLFuncEval,
                     d_global_strategy,
                     d_soln_scale,
                     d_fval_scale,
                     d_residual_tol,
                     d_step_tol,
                     d_constraints,
                     TRUE,                   // Optional inputs
                                             // used by default
                     d_int_optional_input,
                     d_real_optional_input,
                     this);                  // access to user-supplied fcns

   return( retval );

}

/*
*************************************************************************
*                                                                       *
* Setting KINSOL log file name and print flag for KINSOL statistics.    *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setLogFileData(
   const string& log_fname,
   const int flag)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(flag >= 0 && flag <= 3);
#endif
   if ( !(log_fname == d_kinsol_log_file_name) ) {
      if (!log_fname.empty()) {
         d_kinsol_log_file_name = log_fname;
      } else {
         d_kinsol_log_file_name = "kinsol.log";
      }
      d_KINSOL_needs_initialization = true;
   }
   d_int_optional_input[PRINTFL] = flag;
}

/*
*************************************************************************
*                                                                       *
* Accessory functions for setting user-defined function information.    *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setKINSOLFunctions(
   KINSOLAbstractFunctions* my_functions, 
   const int uses_preconditioner, 
   const int uses_jac_times_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(my_functions == (KINSOLAbstractFunctions*)NULL));
#endif

   d_KINSOL_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_uses_jac_times_vector = uses_jac_times_vector;

   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setPreconditioner(const int uses_preconditioner)
{
   d_uses_preconditioner = uses_preconditioner;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setJacobianTimesVector(
   const int uses_jac_times_vector)
{
   d_uses_jac_times_vector = uses_jac_times_vector;
   d_KINSOL_needs_initialization = true;
}

KINSOLAbstractFunctions* KINSOLSolver::getKINSOLFunctions() const
{
   return(d_KINSOL_functions);
}

/*
*************************************************************************
*                                                                       *
* Accessory functions for setting scale vectors.                        *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setSolutionScaleVector(
   PVodeTrioAbstractVector* uscale)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(uscale == (PVodeTrioAbstractVector*)NULL));
#endif
   if (d_my_soln_scale_vector && d_soln_scale) d_soln_scale->freeVector();
   d_soln_scale = uscale;
   d_my_soln_scale_vector = false;
}

void KINSOLSolver::setResidualScaleVector(
   PVodeTrioAbstractVector* fscale)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!(fscale == (PVodeTrioAbstractVector*)NULL));
#endif
   if (d_my_fval_scale_vector && d_fval_scale) d_fval_scale->freeVector();
   d_fval_scale = fscale;
   d_my_fval_scale_vector = false;
}

/*
*************************************************************************
*                                                                       *
* Accessory function for setting constraints for nonlinear system.      *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setConstraintVector(
   PVodeTrioAbstractVector* constraints)
{
   d_constraints = constraints;
}

/*
*************************************************************************
*                                                                       *
* Accessory function for setting nonlinear solver parameters.           *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setResidualStoppingTolerance(const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(tol >= 0.0);
#endif
   d_residual_tol = tol;
}

void KINSOLSolver::setMaxIterations(const int maxits)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(maxits >= 0);
#endif
   d_int_optional_input[MXITER] = maxits;
}

void KINSOLSolver::setMaxKrylovDimension(const int kdim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(kdim >= 0);
#endif
   d_krylov_dimension = kdim;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setGlobalStrategy(const int global)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(global == INEXACT_NEWTON || global == LINESEARCH);
#endif
   d_global_strategy = global;
}

void KINSOLSolver::setMaxNewtonStep(const double maxstep)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(maxstep > 0.0);
#endif
   d_real_optional_input[MXNEWTSTEP] = maxstep;
}

void KINSOLSolver::setNonlinearStepTolerance(const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(tol >= 0.0);
#endif
   d_step_tol = tol;
}

void KINSOLSolver::setRelativeFunctionError(const double reserr)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(reserr > 0.0);
#endif
   d_real_optional_input[RELFUNC] = reserr;
}

void KINSOLSolver::setSolutionUpdateConstraint(const double constraint)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(constraint > 0.0);
#endif
   d_real_optional_input[RELU] = constraint;
}

void KINSOLSolver::setLinearSolverConvergenceTest(const int conv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(conv == ETACONSTANT || conv == ETACHOICE1 || conv == ETACHOICE2);
#endif
   d_int_optional_input[ETACHOICE] = conv;
}

void KINSOLSolver::setEisenstatWalkerParameters(
   const double alpha, const double gamma)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(alpha >= 0.0); 
   assert(gamma >= 0.0);
#endif
   d_real_optional_input[ETAALPHA] = alpha;
   d_real_optional_input[ETAGAMMA] = gamma;
}

void KINSOLSolver::setLinearSolverConstantTolerance(
   const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(tol >= 0.0); 
#endif
   d_real_optional_input[ETACONST] = tol;
}

/*
*************************************************************************
*                                                                       *
* Accessory function for setting preconditioner parameters.             *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setPrecondSetupFlag(const int flag)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(flag == 0 || flag == 1);
#endif
   d_int_optional_input[PRECOND_NO_INIT] = flag;
}

void KINSOLSolver::setMaxStepsWithNoPrecondSetup(const int maxsolv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(maxsolv > 0);
#endif
   d_max_solves_no_set = maxsolv;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxLinearSolveRestarts(const int restarts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(restarts >= 0);
#endif
   d_max_restarts = restarts;
   d_KINSOL_needs_initialization = true;
}


int KINSOLSolver::getTotalNumberOfNonlinearIterations() const
{
   return(d_int_optional_input[NNI]);
}

int KINSOLSolver::getTotalNumberOfFunctionCalls() const
{
   return(d_int_optional_input[NFE]);
}

int KINSOLSolver::getTotalNumberOfBetaConditionFailures() const
{
   return(d_int_optional_input[NBCF]);
}

int KINSOLSolver::getTotalNumberOfBacktracks() const
{
   return(d_int_optional_input[NBKTRK]);
}

double KINSOLSolver::getScaledResidualNorm() const
{
   return(d_real_optional_input[FNORM]);
}

double KINSOLSolver::getNewtonStepLength() const
{
   return(d_real_optional_input[STEPL]);
}

/*
*************************************************************************
*                                                                       *
* Print KINSOLSolver object data to given output stream.        *
*                                                                       *
*************************************************************************
*/
void KINSOLSolver::printClassData(ostream& os) const
{
   os << "\nKINSOLSolver object data members..." << endl;
   os << "this = " << (KINSOLSolver*)this << endl;
   os << "d_solution_vector = " 
      << (PVodeTrioAbstractVector*)d_solution_vector << endl;
   os << "d_soln_scale = " 
      << (PVodeTrioAbstractVector*)d_soln_scale << endl;
   os << "d_fval_scale = " 
      << (PVodeTrioAbstractVector*)d_fval_scale << endl;
   os << "d_my_soln_scale_vector = " << d_my_soln_scale_vector << endl;
   os << "d_my_fval_scale_vector = " << d_my_fval_scale_vector << endl;
   os << "d_constraints = " << (PVodeTrioAbstractVector*)d_constraints 
      << endl;

   os << "d_KINSOL_functions = " 
      << (KINSOLAbstractFunctions*)d_KINSOL_functions << endl;

   os << "d_uses_preconditioner = " << d_uses_preconditioner << endl;
   os << "d_uses_jac_times_vector = " << d_uses_jac_times_vector << endl;

   os << "&d_kin_mem = " << (KINMem*) &d_kin_mem << endl;
   os << "d_kinsol_log_file = " << (FILE*) d_kinsol_log_file << endl;
   os << "d_kinsol_log_file_name = " << d_kinsol_log_file_name << endl;

   os << "d_krylov_dimension = " << d_krylov_dimension << endl;
   os << "d_max_restarts = " << d_max_restarts << endl;
   os << "d_max_solves_no_set = " << d_max_solves_no_set << endl;
   os << "d_neq = " << d_neq << endl;
   os << "d_global_strategy = " << d_global_strategy << endl;
   os << "d_residual_tol = " << d_residual_tol << endl;
   os << "d_step_tol = " << d_step_tol << endl;

   os << "Optional KINSOL inputs/outputs (see KINSOL docs for details):" 
      << endl;
   for (int i = 0; i < OPT_SIZE; i++) {
      os << "d_int_optional_input[" << i << "] = "
         << d_int_optional_input[i] << endl;
   }
   for (int i = 0; i < OPT_SIZE; i++) {
      os << "d_real_optional_input[" << i << "] = "
         << d_real_optional_input[i] << endl;
   }

   os << "...end of KINSOLSolver object data members\n" << endl;
 
}

}
}

#endif
