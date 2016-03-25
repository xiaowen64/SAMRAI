//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/packages/sundials/kinsol/KINSOLSolver.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1996 $
// Modified:    $LastChangedDate: 2008-02-20 08:46:33 -0800 (Wed, 20 Feb 2008) $
// Description: C++ Wrapper class for KINSOL solver package 
//

#include "KINSOLSolver.h"

#include "tbox/Utilities.h"

#ifdef HAVE_SUNDIALS

#include <kinsol/kinsol_impl.h>
#include <kinsol/kinsol_spils.h>

namespace SAMRAI {
    namespace solv {

#define SABSVEC_CAST(v) (static_cast<SAMRAI::solv::SundialsAbstractVector *>(v->content))

/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with KINSOL package.     *
* See header file for KINSOLAbstractFunctions for more information.*
*                                                                       *
*************************************************************************
*/
int KINSOLSolver::KINSOLFuncEval(N_Vector soln,
				  N_Vector fval,
				  void* my_solver)
{
   int success = 0;
   // SGS why no error condition?
   ((KINSOLSolver*)my_solver)->getKINSOLFunctions()->
                                    evaluateNonlinearFunction(SABSVEC_CAST(soln), SABSVEC_CAST(fval));

   return success;
}

int KINSOLSolver::KINSOLPrecondSet(N_Vector uu, 
		     N_Vector uscale,
		     N_Vector fval, 
		     N_Vector fscale,
		     void *my_solver, 
		     N_Vector vtemp1,
		     N_Vector vtemp2)
{
   ((KINSOLSolver*)my_solver)->initializeKINSOL();

   int success = 0;

   int num_feval = 0;
   success = ((KINSOLSolver*)my_solver)->getKINSOLFunctions()->
                                                 precondSetup(SABSVEC_CAST(uu),
                                                              SABSVEC_CAST(uscale),
                                                              SABSVEC_CAST(fval),
                                                              SABSVEC_CAST(fscale),
                                                              SABSVEC_CAST(vtemp1),
                                                              SABSVEC_CAST(vtemp2), 
                                                              num_feval);
   return(success);
}

int KINSOLSolver::KINSOLPrecondSolve(N_Vector uu, 
		       N_Vector uscale, 
		       N_Vector fval, 
		       N_Vector fscale, 
		       N_Vector vv, 
		       void *my_solver,
		       N_Vector vtemp)

{
   int success = 0;
   
   int num_feval = 0;
   success = ((KINSOLSolver*)my_solver)->getKINSOLFunctions()->
                                              precondSolve(SABSVEC_CAST(uu),
                                                           SABSVEC_CAST(uscale),
                                                           SABSVEC_CAST(fval),
                                                           SABSVEC_CAST(fscale),
                                                           SABSVEC_CAST(vv),
                                                           SABSVEC_CAST(vtemp),
                                                           num_feval);
   return(success);
}


int KINSOLSolver::KINSOLJacobianTimesVector(N_Vector v, 
			      N_Vector Jv,
			      N_Vector uu, 
			      int *new_uu, 
			      void *my_solver)
{
   int success = 0;

   bool soln_changed = true; 
   if (*new_uu == 0) {
      soln_changed = false; 
   }

   success = ((KINSOLSolver*)my_solver)->
      getKINSOLFunctions()->jacobianTimesVector(SABSVEC_CAST(v),
						SABSVEC_CAST(Jv),
						soln_changed,
						SABSVEC_CAST(uu));

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
   const std::string& object_name,
   KINSOLAbstractFunctions* my_functions,
   const int uses_preconditioner,
   const int uses_jac_times_vector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!(my_functions == (KINSOLAbstractFunctions*)NULL));
#endif

   d_object_name = object_name;
   d_KINSOL_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_uses_jac_times_vector = uses_jac_times_vector;

   d_KINSOL_needs_initialization = true;

   /*
    * Default parameters to safe values or to KINSOL defaults.
    */

   d_kin_mem              = NULL;
   d_kinsol_log_file      = NULL;
   d_solution_vector      = NULL;
   d_constraints          = NULL;
   d_soln_scale           = NULL;
   d_fval_scale           = NULL;

   d_my_soln_scale_vector = false;
   d_my_fval_scale_vector = false;

   d_krylov_dimension     = KINSPILS_MAXL;

   d_max_restarts         = 0;
   d_max_solves_no_set    = MSBSET_DEFAULT;

   d_no_min_eps           = 0;

   d_max_beta_fails       = MXNBCF_DEFAULT;

   d_no_initial_setup      = 0;
   d_no_residual_monitoring = 0;

   d_global_strategy       = KIN_NONE;
   d_residual_tol          = -1;
   d_step_tol              = -1;

   d_eta_choice            = KIN_ETACONSTANT;
   d_eta_constant          = 0.1;

   d_eta_gamma             = 0.9;
   d_eta_alpha             = 2.0;

   d_omega_min             = 0.00001;
   d_omega_max             = 0.9;

   d_omega                 = 0.0;

   d_max_iter              = MXITER_DEFAULT;
   d_max_newton_step       = -1.0;

   d_maxsub                = MSBSET_SUB_DEFAULT;

   d_relative_function_error = -1;
   
   d_print_level            = 0;

}

void KINSOLSolver::freeInternalVectors(void ) {

   if (d_my_soln_scale_vector && d_my_fval_scale_vector && d_soln_scale) {
      d_soln_scale->freeVector();
      d_soln_scale = NULL;
      d_fval_scale = NULL;
      d_my_soln_scale_vector = false;
      d_my_fval_scale_vector = false;
   }

   if (d_my_soln_scale_vector && d_soln_scale) {
      d_soln_scale->freeVector(); 
      d_soln_scale = NULL;
      d_my_soln_scale_vector = false;
   }

   if (d_my_fval_scale_vector && d_fval_scale) {
      d_fval_scale->freeVector();
      d_fval_scale = NULL;
      d_my_fval_scale_vector = false;
   }
}

KINSOLSolver::~KINSOLSolver()
{

   freeInternalVectors();

   if ( d_kinsol_log_file ) {
     fclose(d_kinsol_log_file);
   }

   if (d_kin_mem) {
      KINFree(&d_kin_mem);
   }
}

/*
*************************************************************************
*                                                                       *
* Functions to initialize nonlinear solver and reset KINSOL structure.  *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::initialize(SundialsAbstractVector* solution,
			      SundialsAbstractVector* uscale, 
			      SundialsAbstractVector* fscale)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(solution == (SundialsAbstractVector*)NULL));
#endif
   
   d_solution_vector = solution;

   // Free previously allocated scaling vectors if
   // KINSOLSolver allocated them.
   freeInternalVectors();

   // If user is providing scaling vectors use them 
   // otherwise allocate them.
   if(uscale) {
      if (d_my_soln_scale_vector && d_soln_scale) {
	 d_soln_scale->freeVector();
      }
      d_soln_scale = uscale;
      d_my_soln_scale_vector = false;
   }

   if(fscale) {
      if (d_my_fval_scale_vector && d_fval_scale) {
	 d_fval_scale->freeVector();
      }
      d_fval_scale = fscale;
      d_my_fval_scale_vector = false;
   }

   // Initialize KINSOL.
   d_KINSOL_needs_initialization = true;

   initializeKINSOL();
}

void KINSOLSolver::initializeKINSOL() 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(d_solution_vector == (SundialsAbstractVector*)NULL));
#endif

   if (d_KINSOL_needs_initialization) {
   
      if (d_kinsol_log_file) {
         fclose(d_kinsol_log_file);
      }

      d_kinsol_log_file = fopen(d_kinsol_log_file_name.c_str(), "w");

      /*
       * KINSOL function pointers.
       */

      KINSpilsPrecSetupFn    precond_set   = NULL;
      KINSpilsPrecSolveFn    precond_solve = NULL;
      KINSpilsJacTimesVecFn  jac_times_vec = NULL;

      if (d_uses_preconditioner) {
         precond_set = KINSOLSolver::KINSOLPrecondSet;
         precond_solve = KINSOLSolver::KINSOLPrecondSolve;
      } else {
         precond_set = NULL;
         precond_solve = NULL;
      }

      if (d_uses_jac_times_vector) {
         jac_times_vec = KINSOLSolver::KINSOLJacobianTimesVector;
      } else {
	 jac_times_vec = NULL;
      }

      if (d_kin_mem) KINFree(&d_kin_mem);

      /*
       * Initialize KINSOL structures and set options
       */
      
      d_kin_mem = KINCreate();


      int ierr = KINMalloc(d_kin_mem, KINSOLSolver::KINSOLFuncEval, d_solution_vector -> getNVector() );
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetFdata(d_kin_mem, this);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetInfoFile(d_kin_mem, d_kinsol_log_file);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetEtaForm(d_kin_mem, d_eta_choice);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetEtaConstValue(d_kin_mem, d_eta_constant);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetEtaParams(d_kin_mem, d_eta_gamma, d_eta_alpha);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetMaxSetupCalls(d_kin_mem, d_max_solves_no_set);
      KINSOL_SAMRAI_ERROR(ierr);


      /*
       * Initialize KINSOL memory record.
       */
      ierr = KINSpgmr( d_kin_mem,
		d_krylov_dimension);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSpilsSetMaxRestarts(d_kin_mem, d_max_restarts);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSpilsSetPreconditioner(d_kin_mem,
				       precond_set,
				       precond_solve,
				       (void*)this);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSpilsSetJacTimesVecFn(d_kin_mem,
				      jac_times_vec,
				      (void*)this);
      KINSOL_SAMRAI_ERROR(ierr);

      if( !(d_residual_tol  < 0) ) {
	 ierr = KINSetFuncNormTol(d_kin_mem, d_residual_tol);
	 KINSOL_SAMRAI_ERROR(ierr);
      }


      if( !(d_step_tol  < 0) ) {
	 ierr = KINSetScaledStepTol(d_kin_mem, d_step_tol);
	 KINSOL_SAMRAI_ERROR(ierr);
      }

      ierr = KINSetConstraints(d_kin_mem, (d_constraints != NULL) ? d_constraints -> getNVector() : NULL );
      KINSOL_SAMRAI_ERROR(ierr);

      // Keep default unless user specifies one.
      if ( !(d_max_newton_step < 0) ) {
	 ierr = KINSetMaxNewtonStep(d_kin_mem, d_max_newton_step);
	 KINSOL_SAMRAI_ERROR(ierr);
      }

      if( !(d_relative_function_error < 0) ) {
	 ierr = KINSetRelErrFunc(d_kin_mem, d_relative_function_error);
	 KINSOL_SAMRAI_ERROR(ierr);
      }

      ierr = KINSetPrintLevel(d_kin_mem, d_print_level);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetConstraints(d_kin_mem, d_constraints == NULL ? NULL : d_constraints -> getNVector());
      KINSOL_SAMRAI_ERROR(ierr);
	
      ierr = KINSetNumMaxIters(d_kin_mem, d_max_iter);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNoInitSetup(d_kin_mem, d_no_initial_setup);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNoResMon(d_kin_mem, d_no_residual_monitoring);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetMaxSubSetupCalls(d_kin_mem, d_maxsub);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetResMonParams(d_kin_mem, d_omega_min, d_omega_max);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetResMonConstValue(d_kin_mem, d_omega);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNoMinEps(d_kin_mem, d_no_min_eps);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetMaxBetaFails(d_kin_mem, d_max_beta_fails);
      KINSOL_SAMRAI_ERROR(ierr);



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

   int retval = KIN_SUCCESS;

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

   retval = KINSol( d_kin_mem,
		    d_solution_vector -> getNVector(),
		    d_global_strategy,
		    d_soln_scale -> getNVector(),
		    d_fval_scale -> getNVector());


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
   const std::string& log_fname,
   const int flag)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(flag >= 0 && flag <= 3);
#endif
   if ( !(log_fname == d_kinsol_log_file_name) ) {
      if (!log_fname.empty()) {
         d_kinsol_log_file_name = log_fname;
      } else {
         d_kinsol_log_file_name = "kinsol.log";
      }
      d_KINSOL_needs_initialization = true;
   }

   d_print_level = flag;
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
   TBOX_ASSERT(!(my_functions == (KINSOLAbstractFunctions*)NULL));
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
* Accessory function for setting constraints for nonlinear system.      *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setConstraintVector(
   SundialsAbstractVector* constraints)
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
   TBOX_ASSERT(tol >= 0.0);
#endif
   d_residual_tol = tol;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxIterations(const int maxits)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxits >= 0);
#endif
   d_max_iter = maxits;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxKrylovDimension(const int kdim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(kdim >= 0);
#endif
   d_krylov_dimension = kdim;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setGlobalStrategy(const int global)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(global == KIN_NONE || global == KIN_LINESEARCH);
#endif
   d_global_strategy = global;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxNewtonStep(const double maxstep)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxstep > 0.0);
#endif
   d_max_newton_step = maxstep;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setNonlinearStepTolerance(const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(tol >= 0.0);
#endif
   d_step_tol = tol;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setRelativeFunctionError(const double reserr)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(reserr > 0.0);
#endif
   d_relative_function_error = reserr;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setLinearSolverConvergenceTest(const int conv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(conv == KIN_ETACONSTANT || conv == KIN_ETACHOICE1 || conv == KIN_ETACHOICE2);
#endif
   d_eta_choice = conv;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setEisenstatWalkerParameters(
   const double alpha, const double gamma)
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

void KINSOLSolver::setLinearSolverConstantTolerance(
   const double tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(tol >= 0.0); 
#endif
   // 
   d_eta_constant = tol;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxSubSetupCalls(const int maxsub) {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxsub >= 0); 
#endif
   d_maxsub = maxsub;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setNoInitialSetup(const bool flag) 
{
   if(flag) {
      d_no_initial_setup = 1;
   } else {
      d_no_initial_setup = 0;
   }
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setNoResidualMonitoring(const bool flag)
{
   if(flag) {
      d_no_residual_monitoring = 1;
   } else {
      d_no_residual_monitoring = 0;
   }
   d_KINSOL_needs_initialization = true;
}


void KINSOLSolver::setResidualMonitoringParams(const double omega_min, const double omega_max)
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

void KINSOLSolver::setResidualMonitoringConstant(const double omega)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(omega >= 0); 
#endif
   d_omega = omega;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setNoMinEps(const bool flag) 
{
   if(flag) {
      d_no_min_eps = 1;
   } else {
      d_no_min_eps = 0;
   }
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxBetaFails(const int max_beta_fails)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_beta_fails >= 0); 
#endif
   d_max_beta_fails = max_beta_fails;
   d_KINSOL_needs_initialization = true;
}


/*
*************************************************************************
*                                                                       *
* Accessory function for setting preconditioner parameters.             *
*                                                                       *
*************************************************************************
*/

void KINSOLSolver::setMaxStepsWithNoPrecondSetup(const int maxsolv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(maxsolv > 0);
#endif
   d_max_solves_no_set = maxsolv;
   d_KINSOL_needs_initialization = true;
}

void KINSOLSolver::setMaxLinearSolveRestarts(const int restarts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(restarts >= 0);
#endif
   d_max_restarts = restarts;
   d_KINSOL_needs_initialization = true;
}


int KINSOLSolver::getTotalNumberOfNonlinearIterations() const
{
   long int num;
   int ierr = KINGetNumNonlinSolvIters(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return num;
}

int KINSOLSolver::getTotalNumberOfFunctionCalls() const
{
   long int num;
   int ierr = KINGetNumFuncEvals(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return num;
}

int KINSOLSolver::getTotalNumberOfBetaConditionFailures() const
{
   long int num;
   int ierr = KINGetNumBetaCondFails(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return num;
}

int KINSOLSolver::getTotalNumberOfBacktracks() const
{
   long int num;
   int ierr = KINGetNumBacktrackOps(d_kin_mem, &num);
   KINSOL_SAMRAI_ERROR(ierr);
   return num;
}

double KINSOLSolver::getScaledResidualNorm() const
{
   realtype norm;
   int ierr = KINGetFuncNorm(d_kin_mem, &norm);
   KINSOL_SAMRAI_ERROR(ierr);
   return norm;
}

double KINSOLSolver::getNewtonStepLength() const
{
   realtype step_length;
   int ierr = KINGetStepLength(d_kin_mem, &step_length);
   KINSOL_SAMRAI_ERROR(ierr);
   return step_length;
}

/*
*************************************************************************
*                                                                       *
* Print KINSOLSolver object data to given output stream.        *
*                                                                       *
*************************************************************************
*/
void KINSOLSolver::printClassData(std::ostream& os) const
{
   os << "\nKINSOLSolver object data members..." << std::endl;
   os << "this = " << (KINSOLSolver*)this << std::endl;
   os << "d_solution_vector = " 
      << (SundialsAbstractVector*)d_solution_vector << std::endl;
   os << "d_soln_scale = " 
      << (SundialsAbstractVector*)d_soln_scale << std::endl;
   os << "d_fval_scale = " 
      << (SundialsAbstractVector*)d_fval_scale << std::endl;
   os << "d_my_soln_scale_vector = " << d_my_soln_scale_vector << std::endl;
   os << "d_my_fval_scale_vector = " << d_my_fval_scale_vector << std::endl;
   os << "d_constraints = " << (SundialsAbstractVector*)d_constraints 
      << std::endl;

   os << "d_KINSOL_functions = " 
      << (KINSOLAbstractFunctions*)d_KINSOL_functions << std::endl;

   os << "d_uses_preconditioner = " << d_uses_preconditioner << std::endl;
   os << "d_uses_jac_times_vector = " << d_uses_jac_times_vector << std::endl;

   os << "d_kin_mem = " << d_kin_mem << std::endl;
   os << "d_kinsol_log_file = " << (FILE*) d_kinsol_log_file << std::endl;
   os << "d_kinsol_log_file_name = " << d_kinsol_log_file_name << std::endl;

   os << "d_krylov_dimension = " << d_krylov_dimension << std::endl;
   os << "d_max_restarts = " << d_max_restarts << std::endl;
   os << "d_max_solves_no_set = " << d_max_solves_no_set << std::endl;
   os << "d_global_strategy = " << d_global_strategy << std::endl;
   os << "d_residual_tol = " << d_residual_tol << std::endl;
   os << "d_step_tol = " << d_step_tol << std::endl;

   // SGS add missing output

   os << "...end of KINSOLSolver object data members\n" << std::endl;
 
}

}
}

#endif
