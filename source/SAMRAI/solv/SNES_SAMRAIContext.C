/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Wrapper for SNES solver for use in a SAMRAI-based application.
 *
 ************************************************************************/

#ifndef included_solv_SNES_SAMRAIContext_C
#define included_solv_SNES_SAMRAIContext_C

#include "SAMRAI/solv/SNES_SAMRAIContext.h"

#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/tbox/RestartManager.h"

#ifdef HAVE_PETSC

namespace SAMRAI {
namespace solv {

const int SNES_SAMRAIContext::SOLV_SNES_SAMRAI_CONTEXT_VERSION = 1;

/*
 *************************************************************************
 *
 * Static member functions that provide linkage with PETSc/SNES package.
 * See header file for SNESAbstractFunctions for more information.
 *
 *************************************************************************
 */

int
SNES_SAMRAIContext::SNESJacobianSet(
   SNES snes,
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
   if (((SNES_SAMRAIContext *)ctx)->getUsesExplicitJacobian()) {
      retval =
         ((SNES_SAMRAIContext *)ctx)->getSNESFunctions()->
         evaluateJacobian(x);
   } else {
      int ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
      PETSC_SAMRAI_ERROR(ierr);
      ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
      PETSC_SAMRAI_ERROR(ierr);
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Constructor and destructor for SNES_SAMRAIContext.  The
 * constructor sets default values for data members, then overrides
 * them with values read from input or restart.  The destructor destroys
 * the SNES object.
 *
 *************************************************************************
 */
SNES_SAMRAIContext::SNES_SAMRAIContext(
   const std::string& object_name,
   SNESAbstractFunctions* my_functions,
   const boost::shared_ptr<tbox::Database>& input_db) :
   d_object_name(object_name),
   d_context_needs_initialization(true),
   d_SNES_solver(((SNES)NULL)),
   d_krylov_solver(((KSP)NULL)),
   d_jacobian(((Mat)NULL)),
   d_preconditioner(((PC)NULL)),
   d_solution_vector(((Vec)NULL)),
   d_residual_vector(((Vec)NULL)),
   d_SNES_functions(my_functions),
   d_uses_preconditioner(true),
   d_uses_explicit_jacobian(true),
   d_maximum_nonlinear_iterations(PETSC_DEFAULT),
   d_maximum_function_evals(PETSC_DEFAULT),
   d_absolute_tolerance(PETSC_DEFAULT),
   d_relative_tolerance(PETSC_DEFAULT),
   d_step_tolerance(PETSC_DEFAULT),
   d_forcing_term_strategy("CONSTANT"),
   d_forcing_term_flag(PETSC_DEFAULT),
   d_constant_forcing_term(PETSC_DEFAULT),
   d_initial_forcing_term(PETSC_DEFAULT),
   d_maximum_forcing_term(PETSC_DEFAULT),
   d_EW_choice2_alpha(PETSC_DEFAULT),
   d_EW_choice2_gamma(PETSC_DEFAULT),
   d_EW_safeguard_exponent(PETSC_DEFAULT),
   d_EW_safeguard_disable_threshold(PETSC_DEFAULT),
   d_SNES_completion_code(SNES_CONVERGED_ITERATING),
   d_linear_solver_absolute_tolerance(PETSC_DEFAULT),
   d_linear_solver_divergence_tolerance(PETSC_DEFAULT),
   d_maximum_linear_iterations(PETSC_DEFAULT),
   d_maximum_gmres_krylov_dimension(PETSC_DEFAULT),
   d_differencing_parameter_strategy(MATMFFD_WP),
   d_function_evaluation_error(PETSC_DEFAULT),
   d_nonlinear_iterations(0)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!(my_functions == (SNESAbstractFunctions *)NULL));

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name,
      this);

   /*
    * Initialize members with data read from the input and restart
    * databases.  Note that PETSc object parameters are set in
    * initialize().
    */
   if (tbox::RestartManager::getManager()->isFromRestart()) {
      getFromRestart();
   }
   getFromInput(input_db);

}

SNES_SAMRAIContext::~SNES_SAMRAIContext()
{
   if (d_solution_vector) {
      PETSc_SAMRAIVectorReal<double>::destroyPETScVector(
         d_solution_vector);
   }

   if (d_residual_vector) {
      PETSc_SAMRAIVectorReal<double>::destroyPETScVector(
         d_residual_vector);
   }

   destroyPetscObjects();
}

/*
 *************************************************************************
 *
 * Routines to initialize PETSc/SNES solver and solve nonlinear system.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::initialize(
   const boost::shared_ptr<SAMRAIVectorReal<double> >& solution)
{
   TBOX_ASSERT(solution);

   /*
    * Set up vectors for solution and nonlinear residual.
    */

   d_solution_vector =
      PETSc_SAMRAIVectorReal<double>::createPETScVector(solution);

   boost::shared_ptr<SAMRAIVectorReal<double> > residual(
      solution->cloneVector("residual"));
   residual->allocateVectorData();
   d_residual_vector =
      PETSc_SAMRAIVectorReal<double>::createPETScVector(residual);

   createPetscObjects();
   initializePetscObjects();
}

/*
 *************************************************************************
 *
 * Reset the state of the nonlinear solver.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::resetSolver(
   const int coarsest_level,
   const int finest_level)
{
   boost::shared_ptr<SAMRAIVectorReal<double> > solution_vector(
      PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(d_solution_vector));
   solution_vector->deallocateVectorData();
   solution_vector->resetLevels(coarsest_level, finest_level);
   solution_vector->allocateVectorData();

   boost::shared_ptr<SAMRAIVectorReal<double> > residual_vector(
      PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(d_residual_vector));
   residual_vector->deallocateVectorData();
   residual_vector->resetLevels(coarsest_level, finest_level);
   residual_vector->allocateVectorData();

   destroyPetscObjects();
   createPetscObjects();
   initializePetscObjects();
}

/*
 *************************************************************************
 *
 * Solve the nonlinear system.
 *
 *************************************************************************
 */
int
SNES_SAMRAIContext::solve()
{
   int ierr;

   if (d_context_needs_initialization) initializePetscObjects();

   Vec initial_guess;

   ierr = VecDuplicate(d_solution_vector, &initial_guess);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = VecSet(initial_guess, 0.0);
   PETSC_SAMRAI_ERROR(ierr);
   ierr = SNESSolve(d_SNES_solver,
         initial_guess,
         d_solution_vector);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetIterationNumber(d_SNES_solver,
         &d_nonlinear_iterations);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetConvergedReason(d_SNES_solver,
         &d_SNES_completion_code);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = VecDestroy(initial_guess);
   PETSC_SAMRAI_ERROR(ierr);

   return ((int)d_SNES_completion_code > 0) ? 1 : 0;
}

/*
 *************************************************************************
 *
 *  Report the reason for termination of nonlinear iterations.  SNES
 *  return codes are translated here, and a message is placed in the
 *  specified output stream.  Test only on relevant completion codes.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::reportCompletionCode(
   std::ostream& os) const
{
   switch ((int)d_SNES_completion_code) {
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
      case SNES_DIVERGED_MAX_IT:
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
 *
 * Create needed Petsc objects and cache a pointer to them.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::createPetscObjects()
{
   int ierr = 0;
   NULL_USE(ierr);

   /*
    * Create the nonlinear solver, specify linesearch backtracking,
    * and register method for nonlinear residual evaluation.
    */
   ierr = SNESCreate(PETSC_COMM_SELF, &d_SNES_solver);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESSetType(d_SNES_solver, SNESLS);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESSetFunction(d_SNES_solver,
         d_residual_vector,
         SNES_SAMRAIContext::SNESFuncEval,
         (void *)this);
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

   ierr = SNESGetKSP(d_SNES_solver, &d_krylov_solver);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = KSPSetPreconditionerSide(d_krylov_solver, PC_RIGHT);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = KSPGetPC(d_krylov_solver, &d_preconditioner);
   PETSC_SAMRAI_ERROR(ierr);

}

/*
 *************************************************************************
 *
 * Initialize the state of cached Petsc objects from cached information.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::initializePetscObjects()
{
   int ierr = 0;
   NULL_USE(ierr);

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

   if (!(d_forcing_term_strategy == "CONSTANT")) {

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
   if (d_jacobian) {
      MatDestroy(d_jacobian);
   }
   if (d_uses_explicit_jacobian) {

      ierr = MatCreateShell(PETSC_COMM_SELF,
            0,                   // dummy number of local rows
            0,                   // dummy number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            (void *)this,
            &d_jacobian);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = MatShellSetOperation(d_jacobian,
            MATOP_MULT,
            (void (*)()) SNES_SAMRAIContext::
            SNESJacobianTimesVector);
      PETSC_SAMRAI_ERROR(ierr);

   } else {

      ierr = MatCreateSNESMF(d_SNES_solver,
            &d_jacobian);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = MatMFFDSetType(
            d_jacobian,
            (MatMFFDType)d_differencing_parameter_strategy.c_str());

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
         SNES_SAMRAIContext::SNESJacobianSet,
         (void *)this);
   PETSC_SAMRAI_ERROR(ierr);

   /*
    * Initialize the Krylov solver object.  This includes setting the
    * type of Krylov method that is used and tolerances used by the
    * method.
    */
   ierr = KSPSetType(d_krylov_solver, (KSPType)d_linear_solver_type.c_str());
   PETSC_SAMRAI_ERROR(ierr);

   if (d_linear_solver_type == "gmres") {

      ierr = KSPGMRESSetRestart(
            d_krylov_solver,
            d_maximum_gmres_krylov_dimension);
      PETSC_SAMRAI_ERROR(ierr);

      if (d_gmres_orthogonalization_algorithm == "modifiedgramschmidt") {

         ierr = KSPGMRESSetOrthogonalization(
               d_krylov_solver,
               KSPGMRESModifiedGramSchmidtOrthogonalization);
         PETSC_SAMRAI_ERROR(ierr);

      } else if (d_gmres_orthogonalization_algorithm ==
                 "gmres_cgs_refine_ifneeded") {

         ierr = KSPGMRESSetCGSRefinementType(
               d_krylov_solver,
               KSP_GMRES_CGS_REFINE_IFNEEDED);
         PETSC_SAMRAI_ERROR(ierr);
      } else if (d_gmres_orthogonalization_algorithm ==
                 "gmres_cgs_refine_always") {

         ierr = KSPGMRESSetCGSRefinementType(
               d_krylov_solver,
               KSP_GMRES_CGS_REFINE_ALWAYS);
         PETSC_SAMRAI_ERROR(ierr);
      }
   }

   if (d_forcing_term_strategy == "CONSTANT") {

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
      ierr = PCSetType(d_preconditioner, (PCType)pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetSetUp(d_preconditioner,
            SNES_SAMRAIContext::SNESsetupPreconditioner);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetContext(d_preconditioner, this);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetApply(d_preconditioner,
            SNES_SAMRAIContext::SNESapplyPreconditioner);
      PETSC_SAMRAI_ERROR(ierr);

   } else {

      std::string pc_type = "none";
      ierr = PCSetType(d_preconditioner, (PCType)pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

   }

   d_context_needs_initialization = false;
}

/*
 *************************************************************************
 *
 * Destroy cached Petsc objects.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::destroyPetscObjects()
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
 *
 * Read parameters from input that are cached in this object.
 *
 *************************************************************************
 */

void
SNES_SAMRAIContext::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{
   if (input_db) {
      if (input_db->keyExists("maximum_nonlinear_iterations")) {
         d_maximum_nonlinear_iterations =
            input_db->getInteger("maximum_nonlinear_iterations");
      }
      if (input_db->keyExists("maximum_function_evals")) {
         d_maximum_function_evals = input_db->getInteger(
            "maximum_function_evals");
      }

      if (input_db->keyExists("uses_preconditioner")) {
         d_uses_preconditioner = input_db->getBool("uses_preconditioner");
      }
      if (input_db->keyExists("uses_explicit_jacobian")) {
         d_uses_explicit_jacobian = input_db->getBool("uses_explicit_jacobian");
      }
      if (input_db->keyExists("absolute_tolerance")) {
         d_absolute_tolerance = input_db->getDouble("absolute_tolerance");
      }
      if (input_db->keyExists("relative_tolerance")) {
         d_relative_tolerance = input_db->getDouble("relative_tolerance");
      }
      if (input_db->keyExists("step_tolerance")) {
         d_step_tolerance = input_db->getDouble("step_tolerance");
      }

      if (input_db->keyExists("forcing_term_strategy")) {
         d_forcing_term_strategy = input_db->getString(
            "forcing_term_strategy");
         if (d_forcing_term_strategy == "EWCHOICE1") {
            d_forcing_term_flag = 1;
         } else if (d_forcing_term_strategy == "EWCHOICE2") {
            d_forcing_term_flag = 2;
         } else if (!(d_forcing_term_strategy == "CONSTANT")) {
            TBOX_ERROR(
               d_object_name << ": "
                             << "Key data `forcing_term_strategy' = "
                             << d_forcing_term_strategy
                             << " in input not recognized.");
         }
      }

      if (input_db->keyExists("constant_forcing_term")) {
         d_constant_forcing_term = input_db->getDouble(
            "constant_forcing_term");
      }
      if (input_db->keyExists("initial_forcing_term")) {
         d_initial_forcing_term = input_db->getDouble("initial_forcing_term");
      }
      if (input_db->keyExists("maximum_forcing_term")) {
         d_maximum_forcing_term = input_db->getDouble("maximum_forcing_term");
      }
      if (input_db->keyExists("EW_choice2_alpha")) {
         d_EW_choice2_alpha = input_db->getDouble("EW_choice2_alpha");
      }
      if (input_db->keyExists("EW_choice2_gamma")) {
         d_EW_choice2_gamma = input_db->getDouble("EW_choice2_gamma");
      }
      if (input_db->keyExists("EW_safeguard_exponent")) {
         d_EW_safeguard_exponent = input_db->getDouble(
            "EW_safeguard_exponent");
      }
      if (input_db->keyExists("EW_safeguard_disable_threshold")) {
         d_EW_safeguard_disable_threshold =
            input_db->getDouble("EW_safeguard_disable_threshold");
      }

      if (input_db->keyExists("linear_solver_type")) {
         d_linear_solver_type = input_db->getString("linear_solver_type");
      }
      if (input_db->keyExists("linear_solver_absolute_tolerance")) {
         d_linear_solver_absolute_tolerance =
            input_db->getDouble("linear_solver_absolute_tolerance");
      }
      if (input_db->keyExists("linear_solver_divergence_tolerance")) {
         d_linear_solver_divergence_tolerance =
            input_db->getDouble("linear_solver_divergence_tolerance");
      }
      if (input_db->keyExists("maximum_linear_iterations")) {
         d_maximum_linear_iterations =
            input_db->getInteger("maximum_linear_iterations");
      }

      if (input_db->keyExists("maximum_gmres_krylov_dimension")) {
         d_maximum_gmres_krylov_dimension =
            input_db->getInteger("maximum_gmres_krylov_dimension");
      }
      if (input_db->keyExists("gmres_orthogonalization_algorithm")) {
         d_gmres_orthogonalization_algorithm =
            input_db->getString("gmres_orthogonalization_algorithm");
      }

      if (input_db->keyExists("differencing_parameter_strategy")) {
         d_differencing_parameter_strategy =
            input_db->getString("differencing_parameter_strategy");
      }
      if (input_db->keyExists("function_evaluation_error")) {
         d_function_evaluation_error =
            input_db->getDouble("function_evaluation_error");
      }
   }

}

/*
 *************************************************************************
 *
 * Routines to read/write from/to restart/database.
 *
 *************************************************************************
 */

void
SNES_SAMRAIContext::getFromRestart()
{

   boost::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file");
   }
   boost::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("SOLV_SNES_SAMRAI_CONTEXT_VERSION");
   if (ver != SOLV_SNES_SAMRAI_CONTEXT_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
                               << "Restart file version different "
                               << "than class version.");
   }

   d_uses_preconditioner = db->getBool("uses_preconditioner");
   d_uses_explicit_jacobian = db->getBool("uses_explicit_jacobian");

   d_maximum_nonlinear_iterations =
      db->getInteger("maximum_nonlinear_iterations");
   d_maximum_function_evals = db->getInteger("maximum_function_evals");

   d_absolute_tolerance = db->getDouble("absolute_tolerance");
   d_relative_tolerance = db->getDouble("relative_tolerance");
   d_step_tolerance = db->getDouble("step_tolerance");

   d_forcing_term_strategy = db->getString("forcing_term_strategy");
   d_forcing_term_flag = db->getInteger("d_forcing_term_flag");

   d_constant_forcing_term = db->getDouble("constant_forcing_term");
   d_initial_forcing_term = db->getDouble("initial_forcing_term");
   d_maximum_forcing_term = db->getDouble("maximum_forcing_term");
   d_EW_choice2_alpha = db->getDouble("EW_choice2_alpha");
   d_EW_choice2_gamma = db->getDouble("EW_choice2_gamma");
   d_EW_safeguard_exponent = db->getDouble("EW_safeguard_exponent");
   d_EW_safeguard_disable_threshold =
      db->getDouble("EW_safeguard_disable_threshold");

   d_linear_solver_type = db->getString("linear_solver_type");
   d_linear_solver_absolute_tolerance =
      db->getDouble("linear_solver_absolute_tolerance");
   d_linear_solver_divergence_tolerance =
      db->getDouble("linear_solver_divergence_tolerance");
   d_maximum_linear_iterations =
      db->getInteger("maximum_linear_iterations");

   d_maximum_gmres_krylov_dimension =
      db->getInteger("maximum_gmres_krylov_dimension");
   d_gmres_orthogonalization_algorithm =
      db->getString("gmres_orthogonalization_algorithm");

   d_function_evaluation_error = db->getDouble("function_evaluation_error");
   d_differencing_parameter_strategy =
      db->getString("differencing_parameter_strategy");

}

void
SNES_SAMRAIContext::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("SOLV_SNES_SAMRAI_CONTEXT_VERSION",
      SOLV_SNES_SAMRAI_CONTEXT_VERSION);

   restart_db->putBool("uses_preconditioner", d_uses_preconditioner);
   restart_db->putBool("uses_explicit_jacobian", d_uses_explicit_jacobian);

   restart_db->putInteger("maximum_nonlinear_iterations",
      d_maximum_nonlinear_iterations);
   restart_db->putInteger("maximum_function_evals",
      d_maximum_function_evals);

   restart_db->putDouble("absolute_tolerance", d_absolute_tolerance);
   restart_db->putDouble("relative_tolerance", d_relative_tolerance);
   restart_db->putDouble("step_tolerance", d_step_tolerance);

   restart_db->putString("forcing_term_strategy", d_forcing_term_strategy);
   restart_db->putInteger("d_forcing_term_flag", d_forcing_term_flag);

   restart_db->putDouble("constant_forcing_term", d_constant_forcing_term);
   restart_db->putDouble("initial_forcing_term", d_initial_forcing_term);
   restart_db->putDouble("maximum_forcing_term", d_maximum_forcing_term);
   restart_db->putDouble("EW_choice2_alpha", d_EW_choice2_alpha);
   restart_db->putDouble("EW_choice2_gamma", d_EW_choice2_gamma);
   restart_db->putDouble("EW_safeguard_exponent", d_EW_safeguard_exponent);
   restart_db->putDouble("EW_safeguard_disable_threshold",
      d_EW_safeguard_disable_threshold);

   restart_db->putString("linear_solver_type", d_linear_solver_type);
   restart_db->putDouble("linear_solver_absolute_tolerance",
      d_linear_solver_absolute_tolerance);
   restart_db->putDouble("linear_solver_divergence_tolerance",
      d_linear_solver_divergence_tolerance);
   restart_db->putInteger("maximum_linear_iterations",
      d_maximum_linear_iterations);

   restart_db->putInteger("maximum_gmres_krylov_dimension",
      d_maximum_gmres_krylov_dimension);
   restart_db->putString("gmres_orthogonalization_algorithm",
      d_gmres_orthogonalization_algorithm);

   restart_db->putDouble("function_evaluation_error",
      d_function_evaluation_error);
   restart_db->putString("differencing_parameter_strategy",
      d_differencing_parameter_strategy);

}

/*
 *************************************************************************
 *
 * Write all class data members to specified output stream.
 *
 *************************************************************************
 */

void
SNES_SAMRAIContext::printClassData(
   std::ostream& os) const
{
   os << "\nSNES_SAMRAIContext::printClassData..." << std::endl;
   os << "SNES_SAMRAIContext: this = "
      << (SNES_SAMRAIContext *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_SNES_functions = "
      << (SNESAbstractFunctions *)d_SNES_functions << std::endl;
   os << "d_SNES_solver = " << (SNES)d_SNES_solver << std::endl;
//   os << "d_SLES_solver = " << (SLES)d_SLES_solver << std::endl;
   os << "d_krylov_solver = " << (KSP)d_krylov_solver << std::endl;
   os << "d_jacobian = " << (Mat *)&d_jacobian << std::endl;
   os << "d_preconditioner = " << (PC *)&d_preconditioner << std::endl;

   os << "d_solution_vector = " << (Vec *)&d_solution_vector << std::endl;
   os << "d_residual_vector = " << (Vec *)&d_residual_vector << std::endl;

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
