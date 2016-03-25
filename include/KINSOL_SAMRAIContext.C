//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/packages/sundials/kinsol/KINSOL_SAMRAIContext.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1941 $
// Modified:    $LastChangedDate: 2008-01-31 08:25:47 -0800 (Thu, 31 Jan 2008) $
// Description: KINSOL solver for use within a SAMRAI-based application.
//
 
#ifndef included_solv_KINSOL_SAMRAIContext_C
#define included_solv_KINSOL_SAMRAIContext_C

#include "KINSOL_SAMRAIContext.h"

#ifdef HAVE_SUNDIALS

#include "Sundials_SAMRAIVector.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#define SOLV_KINSOL_SAMRAI_CONTEXT_VERSION (1)

namespace SAMRAI {
    namespace solv {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for KINSOL_SAMRAIContext<DIM>.  The       *
* constructor sets default values for data members, then overrides      *
* them with values read from input or restart.  The C++ wrapper for     *
* KINSOL is also created in the constructor.  The destructor destroys   *
* the wrappers for KINSOL and the solution vector.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM>  KINSOL_SAMRAIContext<DIM>::KINSOL_SAMRAIContext(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   KINSOLAbstractFunctions* my_functions)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!(my_functions == (KINSOLAbstractFunctions*)NULL));
#endif

   d_object_name = object_name;
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   /* 
    * Set default state.
    */

   d_KINSOL_solver = new KINSOLSolver(object_name,
                                           my_functions,
                                           0, 0);

   d_solution_vector = ((SundialsAbstractVector*)NULL);

   d_residual_stop_tolerance = 0.0;
   d_max_nonlinear_iterations = 0;
   d_max_krylov_dimension = 0;
   d_global_newton_strategy = 0;
   d_max_newton_step = 0.0;
   d_nonlinear_step_tolerance = 0.0;
   d_relative_function_error = 0.0;
   d_solution_update_constraint = 0.0;
   d_linear_convergence_test = 0;
   d_eisenstat_walker_params[0] = 0.0;
   d_eisenstat_walker_params[1] = 0.0;
   d_linear_solver_constant_tolerance = 0.0;
   d_precond_setup_flag = 0;
   d_max_solves_no_precond_setup = 0;
   d_max_linear_solve_restarts = 0;
   d_KINSOL_print_flag = 0;
   d_uses_preconditioner = false;
   d_uses_jac_times_vector = false;

   /*
    * Initialize object with data read from the input and restart databases.
    */

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart ) {
      getFromRestart();
   }
   getFromInput(input_db);

   d_KINSOL_solver->setPreconditioner( ((d_uses_preconditioner == false)
                                       ? 0 : 1 ) );
   d_KINSOL_solver->setJacobianTimesVector( ((d_uses_jac_times_vector == false)
                                            ? 0 : 1) );
}
 
template<int DIM>  KINSOL_SAMRAIContext<DIM>::~KINSOL_SAMRAIContext()
{
   if (d_solution_vector) {
      Sundials_SAMRAIVector<DIM>::destroySundialsVector(d_solution_vector);
   }
   if (d_KINSOL_solver) {
      delete d_KINSOL_solver;
   }
}

/*
*************************************************************************
*                                                                       *
* Routines to initialize KINSOL solver and solve nonlinear system.      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void KINSOL_SAMRAIContext<DIM>::initialize(
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > solution)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!solution.isNull());
#endif
   
   d_solution_vector = Sundials_SAMRAIVector<DIM>::createSundialsVector(solution);
   d_KINSOL_solver->initialize(d_solution_vector);
}

template<int DIM> int KINSOL_SAMRAIContext<DIM>::solve()
{
   return(d_KINSOL_solver->solve());
}

template<int DIM> KINSOLSolver* KINSOL_SAMRAIContext<DIM>::getKINSOLSolver()
{
   return(d_KINSOL_solver);
}

/*
*************************************************************************
*                                                                       *
* Initialize KINSOL solver and solve nonlinear system from input.       *
* Note that all restart values for parameters may be overridden with    *
* input values.                                                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void KINSOL_SAMRAIContext<DIM>::getFromInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   if (db->keyExists("residual_stop_tolerance")) {
      d_residual_stop_tolerance = db->getDouble("residual_stop_tolerance");
      d_KINSOL_solver->setResidualStoppingTolerance(d_residual_stop_tolerance);
   }

   if (db->keyExists("max_nonlinear_iterations")) {
      d_max_nonlinear_iterations = db->getInteger("max_nonlinear_iterations");
      d_KINSOL_solver->setMaxIterations(d_max_nonlinear_iterations);
   }

   if (db->keyExists("max_krylov_dimension")) {
      d_max_krylov_dimension = db->getInteger("max_krylov_dimension");
      d_KINSOL_solver->setMaxKrylovDimension(d_max_krylov_dimension);
   }

   if (db->keyExists("global_newton_strategy")) {
      d_global_newton_strategy = db->getInteger("global_newton_strategy");
      d_KINSOL_solver->setGlobalStrategy(d_global_newton_strategy);
   }
  
   if (db->keyExists("max_newton_step")) {
      d_max_newton_step = db->getDouble("max_newton_step");
      d_KINSOL_solver->setMaxNewtonStep(d_max_newton_step);
   }

   if (db->keyExists("nonlinear_step_tolerance")) {
      d_nonlinear_step_tolerance = db->getDouble("nonlinear_step_tolerance");
      d_KINSOL_solver->setNonlinearStepTolerance(d_nonlinear_step_tolerance);
   }

   if (db->keyExists("relative_function_error")) {
      d_relative_function_error = db->getDouble("relative_function_error");
      d_KINSOL_solver->setRelativeFunctionError(d_relative_function_error);
   }

   if (db->keyExists("linear_convergence_test")) {
      d_linear_convergence_test = db->getInteger("linear_convergence_test");
      d_KINSOL_solver->
         setLinearSolverConvergenceTest(d_linear_convergence_test);
   }

   if (db->keyExists("max_subsetup_calls")) {
      d_max_subsetup_calls = db->getInteger("max_subsetup_calls");
      d_KINSOL_solver->
         setMaxSubSetupCalls(d_max_subsetup_calls);
   }

   if (db->keyExists("residual_monitoring_params")) {
      db->getDoubleArray("residual_monitoring_params",
                         d_residual_monitoring_params, 2);
      d_KINSOL_solver->
         setResidualMonitoringParams(d_residual_monitoring_params[0], 
				     d_residual_monitoring_params[1]);
   }

   if (db->keyExists("residual_monitoring_constant")) {
      d_residual_monitoring_constant = 
	 db->getDouble("residual_monitoring_constant");
      d_KINSOL_solver->
         setResidualMonitoringConstant(d_residual_monitoring_constant);
   }


   if(db->keyExists("no_min_eps")) {
      d_no_min_eps = db->getBool("no_min_eps");
      d_KINSOL_solver-> setNoMinEps(d_no_min_eps);
   }

   if (db->keyExists("eisenstat_walker_params")) {
      db->getDoubleArray("eisenstat_walker_params",
                         d_eisenstat_walker_params, 2);
      d_KINSOL_solver->
         setEisenstatWalkerParameters(d_eisenstat_walker_params[0], 
                                      d_eisenstat_walker_params[1]);
   }

   if (db->keyExists("linear_solver_constant_tolerance")) {
      d_linear_solver_constant_tolerance =
         db->getDouble("linear_solver_constant_tolerance");
      d_KINSOL_solver->
         setLinearSolverConstantTolerance(d_linear_solver_constant_tolerance);
   }

   if (db->keyExists("max_solves_no_precond_setup")) {
      d_max_solves_no_precond_setup = 
         db->getInteger("max_solves_no_precond_setup");
      d_KINSOL_solver->
         setMaxStepsWithNoPrecondSetup(d_max_solves_no_precond_setup);
   }

   if (db->keyExists("max_linear_solve_restarts")) {
      d_max_linear_solve_restarts = 
         db->getInteger("max_linear_solve_restarts");
      d_KINSOL_solver->setMaxLinearSolveRestarts(d_max_linear_solve_restarts);
   }

   bool set_print_options = false;
   if (db->keyExists("KINSOL_log_filename")) {
      d_KINSOL_log_filename = db->getString("KINSOL_log_filename");
      set_print_options = true;
   }

   if (db->keyExists("KINSOL_print_flag")) {
      d_KINSOL_print_flag = db->getInteger("KINSOL_print_flag");
      set_print_options = true;
   }

   if (set_print_options) {
      d_KINSOL_solver->setLogFileData(d_KINSOL_log_filename,
                                      d_KINSOL_print_flag);
   }

   if(db->keyExists("uses_preconditioner")) {
     d_uses_preconditioner = db->getBool("uses_preconditioner");
   }

   if(db->keyExists("uses_jac_times_vector")) {
     d_uses_jac_times_vector = db->getBool("uses_jac_times_vector");
   }

}

/*
*************************************************************************
*                                                                       *
* Read data members from restart database.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void KINSOL_SAMRAIContext<DIM>::getFromRestart()
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

   int ver = db->getInteger("SOLV_KINSOL_SAMRAI_CONTEXT_VERSION");
   if (ver != SOLV_KINSOL_SAMRAI_CONTEXT_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different "
              << "than class version.");
   }

   d_residual_stop_tolerance = db->getDouble("d_residual_stop_tolerance");
   d_max_nonlinear_iterations = db->getInteger("d_max_nonlinear_iterations");
   d_max_krylov_dimension = db->getInteger("d_max_krylov_dimension");
   d_global_newton_strategy = db->getInteger("d_global_newton_strategy");
   d_max_newton_step = db->getDouble("d_max_newton_step");
   d_nonlinear_step_tolerance = db->getDouble("d_nonlinear_step_tolerance");
   d_relative_function_error = db->getDouble("d_relative_function_error");
   d_solution_update_constraint = db->getDouble("d_solution_update_constraint");
   d_linear_convergence_test = db->getInteger("d_linear_convergence_test");
   db->getDoubleArray("d_eisenstat_walker_params",
                      d_eisenstat_walker_params, 2); 
   d_linear_solver_constant_tolerance = 
      db->getDouble("d_linear_solver_constant_tolerance");
   d_precond_setup_flag = db->getInteger("d_precond_setup_flag");
   d_max_solves_no_precond_setup = 
      db->getInteger("d_max_solves_no_precond_setup");
   d_max_linear_solve_restarts = db->getInteger("d_max_linear_solve_restarts");
   d_KINSOL_log_filename = db->getString("d_KINSOL_log_filename");
   d_KINSOL_print_flag = db->getInteger("d_KINSOL_print_flag");

}

/*
*************************************************************************
*                                                                       *
* Write data members to database.                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM> void KINSOL_SAMRAIContext<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("SOLV_KINSOL_SAMRAI_CONTEXT_VERSION",
                  SOLV_KINSOL_SAMRAI_CONTEXT_VERSION);

   db->putDouble("d_residual_stop_tolerance", d_residual_stop_tolerance);
   db->putInteger("d_max_nonlinear_iterations", d_max_nonlinear_iterations);
   db->putInteger("d_max_krylov_dimension", d_max_krylov_dimension);
   db->putInteger("d_global_newton_strategy", d_global_newton_strategy);
   db->putDouble("d_max_newton_step", d_max_newton_step);
   db->putDouble("d_nonlinear_step_tolerance", d_nonlinear_step_tolerance);
   db->putDouble("d_relative_function_error", d_relative_function_error);
   db->putDouble("d_solution_update_constraint", d_solution_update_constraint);
   db->putInteger("d_linear_convergence_test", d_linear_convergence_test);
   db->putDoubleArray("d_eisenstat_walker_params",
                      d_eisenstat_walker_params, 2);
   db->putDouble("d_linear_solver_constant_tolerance",
                 d_linear_solver_constant_tolerance);
   db->putInteger("d_precond_setup_flag", d_precond_setup_flag);
   db->putInteger("d_max_solves_no_precond_setup",
                  d_max_solves_no_precond_setup);
   db->putInteger("d_max_linear_solve_restarts", d_max_linear_solve_restarts);
   db->putString("d_KINSOL_log_filename", d_KINSOL_log_filename);  
   db->putInteger("d_KINSOL_print_flag", d_KINSOL_print_flag);
}

/*
*************************************************************************
*                                                                       *
* Write all class data members to specified output stream.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void KINSOL_SAMRAIContext<DIM>::printClassData(std::ostream& os) const
{
   os << "\nKINSOL_SAMRAIContext<DIM>::printClassData..." << std::endl;
   os << "KINSOL_SAMRAIContext<DIM>: this = "
      << (KINSOL_SAMRAIContext<DIM>*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_KINSOL_solver = " << (KINSOLSolver*)d_KINSOL_solver << std::endl;
   os << "d_solution_vector = " 
      << (SundialsAbstractVector*)d_solution_vector;
   os << "\nd_residual_stop_tolerance = " << d_residual_stop_tolerance << std::endl;
   os << "d_max_nonlinear_iterations = " << d_max_nonlinear_iterations << std::endl;
   os << "d_max_krylov_dimension = " << d_max_krylov_dimension << std::endl;
   os << "d_global_newton_strategy = " << d_global_newton_strategy << std::endl;
   os << "d_max_newton_step = " << d_max_newton_step << std::endl;
   os << "d_nonlinear_step_tolerance = " << d_nonlinear_step_tolerance << std::endl;
   os << "d_relative_function_error = " << d_relative_function_error << std::endl;
   os << "d_solution_update_constraint = " << d_solution_update_constraint;
   os << "\nd_linear_convergence_test = " << d_linear_convergence_test << std::endl;
   os << "d_eisenstat_walker_params = "
      << d_eisenstat_walker_params[0] << " , " 
      << d_eisenstat_walker_params[1] << std::endl;
   os << "d_linear_solver_constant_tolerance = " 
      << d_linear_solver_constant_tolerance << std::endl;
   os << "d_precond_setup_flag = " << d_precond_setup_flag << std::endl;
   os << "d_max_solves_no_precond_setup = " << d_max_solves_no_precond_setup;
   os << "\nd_max_linear_solve_restarts = " << d_max_linear_solve_restarts;
   os << "\nd_KINSOL_log_filename = " << d_KINSOL_log_filename << std::endl;
   os << "d_KINSOL_print_flag = " << d_KINSOL_print_flag << std::endl;

}

}
}

#endif
#endif
