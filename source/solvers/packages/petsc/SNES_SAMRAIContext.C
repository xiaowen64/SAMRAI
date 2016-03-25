//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/packages/petsc/SNES_SAMRAIContext.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2002 $
// Modified:    $LastChangedDate: 2008-02-21 13:54:25 -0800 (Thu, 21 Feb 2008) $
// Description: Wrapper for SNES solver for use in a SAMRAI-based application.
//

#ifndef included_solv_SNES_SAMRAIContext_C
#define included_solv_SNES_SAMRAIContext_C

#include "SNES_SAMRAIContext.h"

#include "PETSc_SAMRAIVectorReal.h"
#include "SAMRAIVectorReal.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#ifdef HAVE_PETSC

#define SOLV_SNES_SAMRAI_CONTEXT_VERSION (1)

namespace SAMRAI {
namespace solv {

/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with PETSc/SNES package. *
* See header file for SNESAbstractFunctions for more information.  *
*                                                                       *
*************************************************************************
*/
template<int DIM> int SNES_SAMRAIContext<DIM>::SNESFuncEval(SNES snes,
							    Vec x,
							    Vec f,
							    void* ctx)
{
   NULL_USE(snes);
   ((SNES_SAMRAIContext<DIM>*)ctx)->getSNESFunctions()->
      evaluateNonlinearFunction(x, f);
   return(0);
}

template<int DIM> int SNES_SAMRAIContext<DIM>::SNESJacobianSet(SNES snes,
							       Vec x,
							       Mat* A,
							       Mat* B,
							       MatStructure* mstruct,
							       void* ctx)
{
   NULL_USE(snes);
   NULL_USE(B);
   NULL_USE(mstruct);
   int retval = 0;
   if ( ((SNES_SAMRAIContext<DIM>*)ctx)->getUsesExplicitJacobian() ) {
      retval =
	    ((SNES_SAMRAIContext<DIM>*)ctx)->getSNESFunctions()->
	 evaluateJacobian(x);
   } else {
      int ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); 
      PETSC_SAMRAI_ERROR(ierr);
      ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);  
      PETSC_SAMRAI_ERROR(ierr);
   }
   return(retval);
}

template<int DIM> int SNES_SAMRAIContext<DIM>::SNESJacobianTimesVector(Mat M,
								       Vec xin,
								       Vec xout)
{
   void* ctx;
   int ierr = 0;
   if (ierr) NULL_STATEMENT;

   ierr = MatShellGetContext(M, &ctx);         PETSC_SAMRAI_ERROR(ierr);
   return(((SNES_SAMRAIContext<DIM>*)ctx)->
	  getSNESFunctions()->jacobianTimesVector(xin, xout));
}

template<int DIM> int SNES_SAMRAIContext<DIM>::SNESsetupPreconditioner(void* ctx)
{
   int ierr = 0;
   Vec current_solution;
   ierr = SNESGetSolution(((SNES_SAMRAIContext<DIM>*)ctx)->getSNESSolver(), 
			  &current_solution);      PETSC_SAMRAI_ERROR(ierr);
   return( ((SNES_SAMRAIContext<DIM>*)ctx)->
	   getSNESFunctions()->setupPreconditioner(current_solution) );
}

template<int DIM> int SNES_SAMRAIContext<DIM>::SNESapplyPreconditioner(void* ctx,
								       Vec r,
								       Vec z)
{
   return( ((SNES_SAMRAIContext<DIM>*)ctx)->
	   getSNESFunctions()->applyPreconditioner(r, z) );
}

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for SNES_SAMRAIContext<DIM>.  The         *
* constructor sets default values for data members, then overrides      *
* them with values read from input or restart.  The destructor destroys *
* the SNES object.                                                      *
*                                                                       *
*************************************************************************
*/
template<int DIM>  SNES_SAMRAIContext<DIM>::SNES_SAMRAIContext(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   SNESAbstractFunctions* my_functions)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!(my_functions == (SNESAbstractFunctions*)NULL));
#endif

   d_object_name = object_name;
   d_context_needs_initialization = true;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, 
							   this);

   /* 
    * Set default state.
    */

   d_SNES_solver    = ((SNES)NULL);
   d_krylov_solver  = ((KSP)NULL);
   d_jacobian       = ((Mat)NULL);
   d_preconditioner = ((PC)NULL);
   d_solution_vector = ((Vec)NULL);
   d_residual_vector = ((Vec)NULL);

   d_SNES_functions = my_functions;

   /*
    * Default nonlinear solver parameters.
    */

   d_absolute_tolerance = PETSC_DEFAULT;
   d_relative_tolerance = PETSC_DEFAULT;
   d_step_tolerance = PETSC_DEFAULT;
   d_maximum_nonlinear_iterations = PETSC_DEFAULT;
   d_maximum_function_evals = PETSC_DEFAULT;

   d_forcing_term_strategy = "CONSTANT";
   d_forcing_term_flag = PETSC_DEFAULT;

   d_constant_forcing_term = PETSC_DEFAULT;
   d_initial_forcing_term = PETSC_DEFAULT;
   d_maximum_forcing_term = PETSC_DEFAULT;
   d_EW_choice2_alpha = PETSC_DEFAULT;
   d_EW_choice2_gamma = PETSC_DEFAULT;
   d_EW_safeguard_exponent = PETSC_DEFAULT;
   d_EW_safeguard_disable_threshold = PETSC_DEFAULT;

   d_SNES_completion_code = SNES_CONVERGED_ITERATING;

   /*
    * Default linear solver parameters.
    */

   d_linear_solver_absolute_tolerance = PETSC_DEFAULT;
   d_linear_solver_divergence_tolerance = PETSC_DEFAULT;
   d_maximum_linear_iterations = PETSC_DEFAULT;

   d_maximum_gmres_krylov_dimension = PETSC_DEFAULT;
   d_gmres_orthogonalization_algorithm = PETSC_DEFAULT;

   /*
    * Default "Matrix-free" parameters.
    */

   d_function_evaluation_error = PETSC_DEFAULT;
   d_differencing_parameter_strategy = MATMFFD_WP;

   /*
    * Default output parameters.
    */

   d_nonlinear_iterations = 0;

   /*
    * Initialize members with data read from the input and restart
    * databases.  Note that PETSc object parameters are set in
    * initialize().
    */

   if  (tbox::RestartManager::getManager()->isFromRestart()) {
      getFromRestart();
   }
   getFromInput(input_db);

}
 
template<int DIM>  SNES_SAMRAIContext<DIM>::~SNES_SAMRAIContext()
{
   if (d_solution_vector) {
      PETSc_SAMRAIVectorReal<DIM,double>::destroyPETScVector(
	 d_solution_vector);
   }

   if (d_residual_vector) {
      PETSc_SAMRAIVectorReal<DIM,double>::destroyPETScVector(
	 d_residual_vector);
   }  

   destroyPetscObjects();
}

/*
*************************************************************************
*                                                                       *
* Access functions for PETSc objects, and user-supplied functions.      *
*                                                                       *
*************************************************************************
*/
template<int DIM> SNES SNES_SAMRAIContext<DIM>::getSNESSolver() const
{
   return(d_SNES_solver);
}

template<int DIM> SNESAbstractFunctions* SNES_SAMRAIContext<DIM>::getSNESFunctions() const
{
   return(d_SNES_functions);
}

template<int DIM> KSP SNES_SAMRAIContext<DIM>::getKrylovSolver() const
{
   return(d_krylov_solver);
}

template<int DIM> Mat SNES_SAMRAIContext<DIM>::getJacobianMatrix() const
{
   return(d_jacobian);
}

/*
*************************************************************************
*                                                                       *
*  Access functions for parameters that control solver behavior.        *
*                                                                       *
*************************************************************************
*/
template<int DIM> double SNES_SAMRAIContext<DIM>::getAbsoluteTolerance() const
{
   return(d_absolute_tolerance);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setAbsoluteTolerance(double abs_tol)
{
   d_absolute_tolerance = abs_tol;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getRelativeTolerance() const
{
   return(d_relative_tolerance);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setRelativeTolerance(double rel_tol)
{
   d_relative_tolerance = rel_tol;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getStepTolerance() const
{
   return(d_step_tolerance);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setStepTolerance(double step_tol)
{
   d_step_tolerance = step_tol;
   d_context_needs_initialization = true;
}

template<int DIM> int SNES_SAMRAIContext<DIM>::getMaxNonlinearIterations() const
{
   return(d_maximum_nonlinear_iterations);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setMaxNonlinearIterations(int max_nli)
{
   d_maximum_nonlinear_iterations = max_nli;
   d_context_needs_initialization = true;
}

template<int DIM> int SNES_SAMRAIContext<DIM>::getMaxFunctionEvaluations() const
{
   return(d_maximum_function_evals);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setMaxFunctionEvaluations(int max_feval)
{
   d_maximum_function_evals = max_feval;
   d_context_needs_initialization = true;
}

template<int DIM> std::string SNES_SAMRAIContext<DIM>::getForcingTermStrategy() const
{
   return(d_forcing_term_strategy);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setForcingTermStrategy(std::string& strategy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( strategy == "CONSTANT"  || 
		strategy == "EWCHOICE1" || 
		strategy == "EWCHOICE2" );
#endif
   d_forcing_term_strategy = strategy;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getConstantForcingTerm() const
{
   return(d_constant_forcing_term);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setConstantForcingTerm(double fixed_eta)
{
   d_constant_forcing_term = fixed_eta;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getInitialForcingTerm() const
{
   return(d_initial_forcing_term);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setInitialForcingTerm(double initial_eta)
{
   d_initial_forcing_term = initial_eta;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getMaximumForcingTerm() const
{
   return(d_maximum_forcing_term);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setMaximumForcingTerm(double max_eta)
{
   d_maximum_forcing_term = max_eta;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getEWChoice2Exponent() const
{
   return(d_EW_choice2_alpha);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setEWChoice2Exponent(double alpha)
{
   d_EW_choice2_alpha = alpha;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getEWChoice2SafeguardExponent() const
{
   return(d_EW_safeguard_exponent);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setEWChoice2SafeguardExponent(double beta)
{
   d_EW_safeguard_exponent = beta;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getEWChoice2ScaleFactor() const
{
   return(d_EW_choice2_gamma);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setEWChoice2ScaleFactor(double gamma)
{
   d_EW_choice2_gamma = gamma;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getEWSafeguardThreshold() const
{
   return(d_EW_safeguard_disable_threshold);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setEWSafeguardThreshold(double threshold)
{
   d_EW_safeguard_disable_threshold = threshold;
   d_context_needs_initialization = true;
}

template<int DIM> std::string SNES_SAMRAIContext<DIM>::getLinearSolverType() const
{
   return(d_linear_solver_type);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setLinearSolverType(std::string& type)
{
   d_linear_solver_type = type;
   d_context_needs_initialization = true;
}

template<int DIM> bool SNES_SAMRAIContext<DIM>::getUsesPreconditioner() const
{
   return(d_uses_preconditioner);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setUsesPreconditioner(bool uses_preconditioner)
{
   d_uses_preconditioner = uses_preconditioner;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getLinearSolverAbsoluteTolerance() const
{
   return(d_linear_solver_absolute_tolerance);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setLinearSolverAbsoluteTolerance(double abs_tol)
{
   d_linear_solver_absolute_tolerance = abs_tol;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getLinearSolverDivergenceTolerance() const
{
   return(d_linear_solver_divergence_tolerance);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setLinearSolverDivergenceTolerance(
   double div_tol)
{
   d_linear_solver_divergence_tolerance = div_tol;
   d_context_needs_initialization = true;
}

template<int DIM> int SNES_SAMRAIContext<DIM>::getMaximumLinearIterations() const
{
   return(d_maximum_linear_iterations);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setMaximumLinearIterations(int max_li)
{
   d_maximum_linear_iterations = max_li;
   d_context_needs_initialization = true;
}

template<int DIM> int SNES_SAMRAIContext<DIM>::getMaximumGMRESKrylovDimension() const
{
   return(d_maximum_gmres_krylov_dimension);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setMaximumGMRESKrylovDimension(int d)
{
   d_maximum_gmres_krylov_dimension = d;
   d_context_needs_initialization = true;
}

template<int DIM> std::string SNES_SAMRAIContext<DIM>::getGMRESOrthogonalizationMethod() const
{
   return(d_gmres_orthogonalization_algorithm);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setGMRESOrthogonalizationMethod(std::string& method)
{
   d_gmres_orthogonalization_algorithm = method;
   d_context_needs_initialization = true;
}

template<int DIM> bool SNES_SAMRAIContext<DIM>::getUsesExplicitJacobian() const
{
   return(d_uses_explicit_jacobian);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setUsesExplicitJacobian(bool use_jac)
{
   d_uses_explicit_jacobian = use_jac;
   d_context_needs_initialization = true;
}

template<int DIM> std::string SNES_SAMRAIContext<DIM>::getDifferencingParameterMethod() const
{
   return(d_differencing_parameter_strategy);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setDifferencingParameterMethod(std::string& method)
{
   d_differencing_parameter_strategy = method;
   d_context_needs_initialization = true;
}

template<int DIM> double SNES_SAMRAIContext<DIM>::getFunctionEvaluationError() const
{
   return(d_function_evaluation_error);
}

template<int DIM> void SNES_SAMRAIContext<DIM>::setFunctionEvaluationError(
   double evaluation_error)
{
   d_function_evaluation_error = evaluation_error;
   d_context_needs_initialization = true;
}

/*
*************************************************************************
*                                                                       *
* Routines to initialize PETSc/SNES solver and solve nonlinear system.  *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SNES_SAMRAIContext<DIM>::initialize(
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > solution)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!solution.isNull());
#endif

   /*
    * Set up vectors for solution and nonlinear residual.
    */

   d_solution_vector =
      PETSc_SAMRAIVectorReal<DIM,double>::createPETScVector(solution);

   tbox::Pointer< SAMRAIVectorReal<DIM,double> > residual =
      solution->cloneVector("residual");
   residual->allocateVectorData();
   d_residual_vector = 
      PETSc_SAMRAIVectorReal<DIM,double>::createPETScVector(residual);

   createPetscObjects();
   initializePetscObjects();
}

/*
*************************************************************************
*                                                                       *
* Reset the state of the nonlinear solver.                              *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SNES_SAMRAIContext<DIM>::resetSolver(const int coarsest_level,
							    const int finest_level)
{
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > solution_vector =
      PETSc_SAMRAIVectorReal<DIM,double>::getSAMRAIVector(
	 d_solution_vector);
   solution_vector->deallocateVectorData();
   solution_vector->resetLevels(coarsest_level, finest_level);
   solution_vector->allocateVectorData();

   tbox::Pointer< SAMRAIVectorReal<DIM,double> > residual_vector =
      PETSc_SAMRAIVectorReal<DIM,double>::getSAMRAIVector(
	 d_residual_vector);
   residual_vector->deallocateVectorData();
   residual_vector->resetLevels(coarsest_level, finest_level);
   residual_vector->allocateVectorData();

   destroyPetscObjects();
   createPetscObjects();
   initializePetscObjects();
}

/*
*************************************************************************
*                                                                       *
* Solve the nonlinear system.                                           *
*                                                                       *
*************************************************************************
*/
template<int DIM> int SNES_SAMRAIContext<DIM>::solve()
{
   int ierr;

   if ( d_context_needs_initialization ) initializePetscObjects();

   Vec initial_guess;

   ierr = VecDuplicate(d_solution_vector, &initial_guess);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = VecSet( initial_guess, 0.0); 
   PETSC_SAMRAI_ERROR(ierr);
   
   ierr = SNESSolve(d_SNES_solver,
		    initial_guess,
		    d_solution_vector); 
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetLinearSolveIterations(d_SNES_solver,
				       &d_nonlinear_iterations);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetConvergedReason(d_SNES_solver,
				 &d_SNES_completion_code);  
   PETSC_SAMRAI_ERROR(ierr);

   ierr = VecDestroy(initial_guess);
   PETSC_SAMRAI_ERROR(ierr);

   return(((int)d_SNES_completion_code > 0) ? 1 : 0);
}

/*
*************************************************************************
*                                                                       *
* Get the number of nonlinear iterations used in last solve.            *
*                                                                       *
*************************************************************************
*/
template<int DIM> int SNES_SAMRAIContext<DIM>::getNonlinearIterationCount() const
{
   return(d_nonlinear_iterations);
}

/*
*************************************************************************
*                                                                       *
* Get the total number of linear iterations used in last solve.         *
*                                                                       *
*************************************************************************
*/
template<int DIM> int SNES_SAMRAIContext<DIM>::getTotalLinearIterationCount() const
{
   int linear_itns;
   int ierr = SNESGetLinearSolveIterations(d_SNES_solver, 
					   &linear_itns); 
   PETSC_SAMRAI_ERROR(ierr);
   return(linear_itns);
}

/*
*************************************************************************
*                                                                       *
*  Report the reason for termination of nonlinear iterations.  SNES     *
*  return codes are translated here, and a message is placed in the     *
*  specified output stream.  Test only on relevant completion codes.    *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SNES_SAMRAIContext<DIM>::reportCompletionCode(std::ostream & os) const
{
   switch ((int) d_SNES_completion_code) {
      case SNES_CONVERGED_FNORM_ABS:
	 os << " Fnorm less than specified absolute tolerance.\n";
	 break;
      case SNES_CONVERGED_FNORM_RELATIVE:
	 os << " Fnorm less than specified relative tolerance.\n";
	 break;
      case SNES_CONVERGED_PNORM_RELATIVE:
	 os << " Step size less than specified tolerance.\n";
	 break;
      case SNES_DIVERGED_FUNCTION_COUNT:
	 os << " Maximum function evaluation count exceeded.\n";
	 break;
      case SNES_DIVERGED_FNORM_NAN:
	 os << " Norm of F is NAN.\n";
	 break;
      case SNES_DIVERGED_MAX_IT	:
	 os << " Maximum nonlinear iteration count exceeded.\n";
	 break;
      case SNES_DIVERGED_LS_FAILURE:
	 os << " Failure in linesearch procedure.\n";
	 break;
      default:
	 os << " Inappropriate completion code reported.\n";
	 break;
   }
}

/*
*************************************************************************
*                                                                       *
* Create needed Petsc objects and cache a pointer to them.              *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SNES_SAMRAIContext<DIM>::createPetscObjects()
{
   int ierr = 0;
   if (ierr) NULL_STATEMENT;

   /*
    * Create the nonlinear solver, specify linesearch backtracking,
    * and register method for nonlinear residual evaluation.
    */
   ierr = SNESCreate(PETSC_COMM_SELF,
		     &d_SNES_solver); 
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESSetType(d_SNES_solver, 
		      SNESLS); 
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESSetFunction(d_SNES_solver,
			  d_residual_vector,
			  SNES_SAMRAIContext<DIM>::SNESFuncEval,
			  (void*) this);
   PETSC_SAMRAI_ERROR(ierr);
   /*
    * Cache the linear solver object, as well as the wrapped Krylov
    * solver and preconditioner.
    */
//   ierr = SNESGetSLES(d_SNES_solver,
//                      &d_SLES_solver);
//                      PETSC_SAMRAI_ERROR(ierr);

//   ierr = SLESGetKSP(d_SLES_solver,
//                     &d_krylov_solver);
//                     PETSC_SAMRAI_ERROR(ierr);



   ierr = SNESGetKSP(d_SNES_solver,
		     &d_krylov_solver);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = KSPSetPreconditionerSide(d_krylov_solver,
				   PC_RIGHT);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = KSPGetPC(d_krylov_solver,
		   &d_preconditioner);
   PETSC_SAMRAI_ERROR(ierr);


}

/*
*************************************************************************
*                                                                       *
* Initialize the state of cached Petsc objects from cached information. *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SNES_SAMRAIContext<DIM>::initializePetscObjects()
{
   int ierr = 0;
   if (ierr) NULL_STATEMENT;

   /*
    * Set tolerances in nonlinear solver.  Also set parameters if 
    * the Jacobian-free option has been selected.
    */
   ierr = SNESSetTolerances(d_SNES_solver,
			    d_absolute_tolerance,
			    d_relative_tolerance,
			    d_step_tolerance,
			    d_maximum_nonlinear_iterations,
			    d_maximum_function_evals);
   PETSC_SAMRAI_ERROR(ierr);

   if (!(d_forcing_term_strategy == "CONSTANT") ) {

      ierr = SNESKSPSetUseEW(d_SNES_solver, PETSC_TRUE);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = SNESKSPSetParametersEW(d_SNES_solver,
				    d_forcing_term_flag,
				    d_initial_forcing_term,
				    d_maximum_forcing_term,
				    d_EW_choice2_gamma,
				    d_EW_choice2_alpha,
				    d_EW_safeguard_exponent,
				    d_EW_safeguard_disable_threshold);
      PETSC_SAMRAI_ERROR(ierr);
   }

   /*
    * Create data structures needed for Jacobian.  This is done
    * here in case an application toggles use of an explicit
    * Jacobian within a run.  
    *
    * First delete any Jacobian object that already has been created.
    */
   if (d_jacobian) MatDestroy(d_jacobian);
   if (d_uses_explicit_jacobian) {

      ierr = MatCreateShell(PETSC_COMM_SELF,
			    0,   // dummy number of local rows
			    0,   // dummy number of local columns
			    PETSC_DETERMINE,
			    PETSC_DETERMINE,
			    (void*) this,
			    &d_jacobian);
      PETSC_SAMRAI_ERROR(ierr);
    
      ierr = MatShellSetOperation(d_jacobian, 
				  MATOP_MULT,
				  (void(*)()) SNES_SAMRAIContext<DIM>::
				  SNESJacobianTimesVector);
      PETSC_SAMRAI_ERROR(ierr); 

   } else {

      ierr = MatCreateSNESMF(d_SNES_solver,
			     &d_jacobian);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = MatMFFDSetType(
	 d_jacobian,
	 (MatMFFDType)d_differencing_parameter_strategy.c_str() );

      ierr = MatMFFDSetFunctionError(d_jacobian, 
				     d_function_evaluation_error);
   }

   /*
    * Register method for setting up Jacobian; this is the same
    * for both options.
    *
    * N.B.  In principle, the second Mat argument should not 
    * be the same as the first Mat argument.  However we 
    * restrict to either no preconditioner, or a shell 
    * preconditioner; in these circumstances that seems to 
    * cause no problem, since the shell preconditioner provides 
    * its own setup method.
    */
   ierr = SNESSetJacobian(d_SNES_solver,
			  d_jacobian,
			  d_jacobian,
			  SNES_SAMRAIContext<DIM>::SNESJacobianSet,
			  (void*) this);
   PETSC_SAMRAI_ERROR(ierr);

   /*
    * Initialize the Krylov solver object.  This includes setting the
    * type of Krylov method that is used and tolerances used by the 
    * method.
    */
   ierr = KSPSetType(d_krylov_solver,
		     (KSPType) d_linear_solver_type.c_str());
   PETSC_SAMRAI_ERROR(ierr);

   if ( d_linear_solver_type == "gmres" ) {

      ierr = KSPGMRESSetRestart(
	 d_krylov_solver,
	 d_maximum_gmres_krylov_dimension);
      PETSC_SAMRAI_ERROR(ierr);

      if (d_gmres_orthogonalization_algorithm == "modifiedgramschmidt") {

	 ierr = KSPGMRESSetOrthogonalization(
	    d_krylov_solver,
	    KSPGMRESModifiedGramSchmidtOrthogonalization); 
	 PETSC_SAMRAI_ERROR(ierr);

      } else if (d_gmres_orthogonalization_algorithm == "gmres_cgs_refine_ifneeded") {

	 ierr = KSPGMRESSetCGSRefinementType(
	    d_krylov_solver,
	    KSP_GMRES_CGS_REFINE_IFNEEDED);
	 PETSC_SAMRAI_ERROR(ierr);
      } else if (d_gmres_orthogonalization_algorithm == "gmres_cgs_refine_always") {

	 ierr = KSPGMRESSetCGSRefinementType(
	    d_krylov_solver,
	    KSP_GMRES_CGS_REFINE_ALWAYS);
	 PETSC_SAMRAI_ERROR(ierr);
      }
   }

   if ( d_forcing_term_strategy == "CONSTANT" ) {

      ierr = KSPSetTolerances(d_krylov_solver,
			      d_constant_forcing_term,
			      d_linear_solver_absolute_tolerance,
			      d_linear_solver_divergence_tolerance,
			      d_maximum_linear_iterations);
      PETSC_SAMRAI_ERROR(ierr);
   } 

   /*
    * Initialize the precondtioner.  Only shell PCs are supported.
    * For these, register the methods used to set up and apply
    * the preconditioner.
    */
   if (d_uses_preconditioner) {

      std::string pc_type = "shell";
      ierr = PCSetType(d_preconditioner,
		       (PCType) pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetSetUp(
	 d_preconditioner,
	 SNES_SAMRAIContext<DIM>::SNESsetupPreconditioner);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetContext(d_preconditioner, this);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetApply(d_preconditioner,
			     SNES_SAMRAIContext<DIM>::SNESapplyPreconditioner); 
      PETSC_SAMRAI_ERROR(ierr); 

   } else {
    
      std::string pc_type = "none";
      ierr = PCSetType(d_preconditioner,
		       (PCType) pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

   }

   d_context_needs_initialization = false;
}

/*
*************************************************************************
*                                                                       *
* Destroy cached Petsc objects.                                         *
*                                                                       *
*************************************************************************
*/
template<int DIM> void SNES_SAMRAIContext<DIM>::destroyPetscObjects()
{
   if (d_jacobian) {
      MatDestroy(d_jacobian);
      d_jacobian = ((Mat)NULL);
   }

   if (d_SNES_solver) {
      SNESDestroy(d_SNES_solver);
//     if (d_SLES_solver) d_SLES_solver = ((SLES)NULL);
      if (d_preconditioner) d_preconditioner = ((PC)NULL);
      if (d_krylov_solver) d_krylov_solver = ((KSP)NULL);
   }
}

/*
*************************************************************************
*                                                                       *
* Read parameters from input that are cached in this object.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void SNES_SAMRAIContext<DIM>::getFromInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   if (db->keyExists("maximum_nonlinear_iterations")) {
      d_maximum_nonlinear_iterations = 
	 db->getInteger("maximum_nonlinear_iterations");
   }
   if (db->keyExists("maximum_function_evals")) {
      d_maximum_function_evals = db->getInteger("maximum_function_evals");
   }

   if (db->keyExists("uses_preconditioner")) {
      d_uses_preconditioner = db->getBool("uses_preconditioner");
   }
   if (db->keyExists("uses_explicit_jacobian")) {
      d_uses_explicit_jacobian = db->getBool("uses_explicit_jacobian");
   }
   if (db->keyExists("absolute_tolerance")) {
      d_absolute_tolerance = db->getDouble("absolute_tolerance");
   }
   if (db->keyExists("relative_tolerance")) {
      d_relative_tolerance = db->getDouble("relative_tolerance");
   }
   if (db->keyExists("step_tolerance")) {
      d_step_tolerance = db->getDouble("step_tolerance");
   }

   if (db->keyExists("forcing_term_strategy")) {
      d_forcing_term_strategy = db->getString("forcing_term_strategy");
      if (d_forcing_term_strategy == "EWCHOICE1") {
	 d_forcing_term_flag = 1; 
      } else if (d_forcing_term_strategy == "EWCHOICE2") {
	 d_forcing_term_flag = 2;
      } else if ( !(d_forcing_term_strategy == "CONSTANT") ) {
	 TBOX_ERROR(d_object_name << ": "
		    << "Key data `forcing_term_strategy' = "
		    << d_forcing_term_strategy << " in input not recognized.");
      }
   }

   if (db->keyExists("constant_forcing_term")) {
      d_constant_forcing_term = db->getDouble("constant_forcing_term");
   }
   if (db->keyExists("initial_forcing_term")) {
      d_initial_forcing_term = db->getDouble("initial_forcing_term");
   }
   if (db->keyExists("maximum_forcing_term")) {
      d_maximum_forcing_term = db->getDouble("maximum_forcing_term");
   }
   if (db->keyExists("EW_choice2_alpha")) {
      d_EW_choice2_alpha = db->getDouble("EW_choice2_alpha");
   }
   if (db->keyExists("EW_choice2_gamma")) {
      d_EW_choice2_gamma = db->getDouble("EW_choice2_gamma");
   }
   if (db->keyExists("EW_safeguard_exponent")) {
      d_EW_safeguard_exponent = db->getDouble("EW_safeguard_exponent");
   }
   if (db->keyExists("EW_safeguard_disable_threshold")) {
      d_EW_safeguard_disable_threshold =
	 db->getDouble("EW_safeguard_disable_threshold");
   }

   if (db->keyExists("linear_solver_type")) {
      d_linear_solver_type = db->getString("linear_solver_type");
   }
   if (db->keyExists("linear_solver_absolute_tolerance")) {
      d_linear_solver_absolute_tolerance =
	 db->getDouble("linear_solver_absolute_tolerance");
   }
   if (db->keyExists("linear_solver_divergence_tolerance")) {
      d_linear_solver_divergence_tolerance =
	 db->getDouble("linear_solver_divergence_tolerance");
   }
   if (db->keyExists("maximum_linear_iterations")) {
      d_maximum_linear_iterations = 
	 db->getInteger("maximum_linear_iterations");
   }

   if (db->keyExists("maximum_gmres_krylov_dimension")) {
      d_maximum_gmres_krylov_dimension = 
	 db->getInteger("maximum_gmres_krylov_dimension");		
   }
   if (db->keyExists("gmres_orthogonalization_algorithm")) {
      d_gmres_orthogonalization_algorithm = 
	 db->getString("gmres_orthogonalization_algorithm");
   }

   if (db->keyExists("differencing_parameter_strategy")) {
      d_differencing_parameter_strategy =
	 db->getString("differencing_parameter_strategy");
   }
   if (db->keyExists("function_evaluation_error")) {
      d_function_evaluation_error = 
	 db->getDouble("function_evaluation_error");
   }

}

/*
*************************************************************************
*                                                                       *
* Routines to read/write from/to restart/database.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void SNES_SAMRAIContext<DIM>::getFromRestart()
{

   tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
		 << d_object_name << " not found in restart file");
   }

   int ver = db->getInteger("SOLV_SNES_SAMRAI_CONTEXT_VERSION");
   if (ver != SOLV_SNES_SAMRAI_CONTEXT_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
		 << "Restart file version different "
		 << "than class version.");
   }

   d_uses_preconditioner = db->getBool("d_uses_preconditioner");
   d_uses_explicit_jacobian = db->getBool("d_uses_explicit_jacobian");

   d_maximum_nonlinear_iterations =
      db->getInteger("d_maximum_nonlinear_iterations");
   d_maximum_function_evals = db->getInteger("d_maximum_function_evals");

   d_absolute_tolerance = db->getDouble("d_absolute_tolerance");
   d_relative_tolerance = db->getDouble("d_relative_tolerance");
   d_step_tolerance = db->getDouble("d_step_tolerance");

   d_forcing_term_strategy = db->getString("d_forcing_term_strategy");
   d_forcing_term_flag = db->getInteger("d_forcing_term_flag");

   d_constant_forcing_term = db->getDouble("d_constant_forcing_term");
   d_initial_forcing_term = db->getDouble("d_initial_forcing_term");
   d_maximum_forcing_term = db->getDouble("d_maximum_forcing_term");
   d_EW_choice2_alpha = db->getDouble("d_EW_choice2_alpha");
   d_EW_choice2_gamma = db->getDouble("d_EW_choice2_gamma");
   d_EW_safeguard_exponent = db->getDouble("d_EW_safeguard_exponent");
   d_EW_safeguard_disable_threshold =
      db->getDouble("d_EW_safeguard_disable_threshold");

   d_linear_solver_type = db->getString("d_linear_solver_type");
   d_linear_solver_absolute_tolerance = 
      db->getDouble("d_linear_solver_absolute_tolerance");
   d_linear_solver_divergence_tolerance = 
      db->getDouble("d_linear_solver_divergence_tolerance");
   d_maximum_linear_iterations = 
      db->getInteger("d_maximum_linear_iterations");

   d_maximum_gmres_krylov_dimension = 
      db->getInteger("d_maximum_gmres_krylov_dimension");
   d_gmres_orthogonalization_algorithm = 
      db->getString("d_gmres_orthogonalization_algorithm");

   d_function_evaluation_error = db->getDouble("d_function_evaluation_error");
   d_differencing_parameter_strategy = 
      db->getString("d_differencing_parameter_strategy");

}

template<int DIM> void SNES_SAMRAIContext<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("SOLV_SNES_SAMRAI_CONTEXT_VERSION",
		  SOLV_SNES_SAMRAI_CONTEXT_VERSION);

   db->putBool("d_uses_preconditioner", d_uses_preconditioner);
   db->putBool("d_uses_explicit_jacobian", d_uses_explicit_jacobian);

   db->putInteger("d_maximum_nonlinear_iterations",
		  d_maximum_nonlinear_iterations);
   db->putInteger("d_maximum_function_evals", d_maximum_function_evals);

   db->putDouble("d_absolute_tolerance", d_absolute_tolerance);
   db->putDouble("d_relative_tolerance", d_relative_tolerance);
   db->putDouble("d_step_tolerance", d_step_tolerance);

   db->putString("d_forcing_term_strategy", d_forcing_term_strategy);
   db->putInteger("d_forcing_term_flag", d_forcing_term_flag);

   db->putDouble("d_constant_forcing_term", d_constant_forcing_term);
   db->putDouble("d_initial_forcing_term", d_initial_forcing_term);
   db->putDouble("d_maximum_forcing_term", d_maximum_forcing_term);
   db->putDouble("d_EW_choice2_alpha", d_EW_choice2_alpha);
   db->putDouble("d_EW_choice2_gamma", d_EW_choice2_gamma);
   db->putDouble("d_EW_safeguard_exponent", d_EW_safeguard_exponent);
   db->putDouble("d_EW_safeguard_disable_threshold",
		 d_EW_safeguard_disable_threshold);

   db->putString("d_linear_solver_type", d_linear_solver_type);
   db->putDouble("d_linear_solver_absolute_tolerance", 
		 d_linear_solver_absolute_tolerance);
   db->putDouble("d_linear_solver_divergence_tolerance",
		 d_linear_solver_divergence_tolerance);
   db->putInteger("d_maximum_linear_iterations", 
		  d_maximum_linear_iterations);

   db->putInteger("d_maximum_gmres_krylov_dimension",
		  d_maximum_gmres_krylov_dimension);
   db->putString("d_gmres_orthogonalization_algorithm",
		 d_gmres_orthogonalization_algorithm);

   db->putDouble("d_function_evaluation_error",
		 d_function_evaluation_error);
   db->putString("d_differencing_parameter_strategy",
		 d_differencing_parameter_strategy);

}

/*
*************************************************************************
*                                                                       *
* Write all class data members to specified output stream.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void SNES_SAMRAIContext<DIM>::printClassData(std::ostream& os) const
{
   os << "\nSNES_SAMRAIContext<DIM>::printClassData..." << std::endl;
   os << "SNES_SAMRAIContext<DIM>: this = "
      << (SNES_SAMRAIContext<DIM>*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_SNES_functions = " 
      << (SNESAbstractFunctions*)d_SNES_functions << std::endl; 
   os << "d_SNES_solver = " << (SNES)d_SNES_solver << std::endl; 
//   os << "d_SLES_solver = " << (SLES)d_SLES_solver << std::endl; 
   os << "d_krylov_solver = " << (KSP)d_krylov_solver << std::endl; 
   os << "d_jacobian = " << (Mat*)&d_jacobian << std::endl; 
   os << "d_preconditioner = " << (PC*)&d_preconditioner << std::endl; 

   os << "d_solution_vector = " << (Vec*)&d_solution_vector << std::endl; 
   os << "d_residual_vector = " << (Vec*)&d_residual_vector << std::endl; 

   os << "d_uses_preconditioner = " << d_uses_preconditioner << std::endl;
   os << "d_uses_explicit_jacobian = " << d_uses_explicit_jacobian << std::endl;

   os << "d_maximum_nonlinear_iterations = " 
      << d_maximum_nonlinear_iterations << std::endl;
   os << "d_maximum_function_evals = " << d_maximum_function_evals << std::endl;

   os << "d_absolute_tolerance = " << d_absolute_tolerance << std::endl;
   os << "d_relative_tolerance = " << d_relative_tolerance << std::endl;
   os << "d_step_tolerance = " << d_step_tolerance << std::endl;

   os << "d_forcing_term_strategy = " << d_forcing_term_strategy << std::endl;
   os << "d_forcing_term_flag = " << d_forcing_term_flag << std::endl;
 
   os << "d_constant_forcing_term = " << d_constant_forcing_term << std::endl;
   os << "d_initial_forcing_term = " << d_initial_forcing_term << std::endl;
   os << "d_EW_choice2_alpha = " << d_EW_choice2_alpha << std::endl;
   os << "d_EW_choice2_gamma = " << d_EW_choice2_gamma << std::endl;
   os << "d_EW_safeguard_exponent = " << d_EW_safeguard_exponent << std::endl;
   os << "d_EW_safeguard_disable_threshold = " 
      << d_EW_safeguard_disable_threshold << std::endl;

   os << "d_linear_solver_type = " << d_linear_solver_type << std::endl;
   os << "d_linear_solver_absolute_tolerance = " 
      << d_linear_solver_absolute_tolerance << std::endl;
   os << "d_linear_solver_divergence_tolerance = " 
      << d_linear_solver_divergence_tolerance << std::endl;
   os << "d_maximum_linear_iterations = " 
      << d_maximum_linear_iterations << std::endl;

   os << "d_maximum_gmres_krylov_dimension = " 
      << d_maximum_gmres_krylov_dimension << std::endl;
   os << "d_gmres_orthogonalization_algorithm = " 
      << d_gmres_orthogonalization_algorithm << std::endl;

   os << "d_differencing_parameter_strategy = " 
      << d_differencing_parameter_strategy << std::endl;
   os << "d_function_evaluation_error = " 
      << d_function_evaluation_error << std::endl;
}

}
}

#endif
#endif
