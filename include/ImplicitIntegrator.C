//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/implicit/ImplicitIntegrator.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2122 $
// Modified:    $LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description: Implicit time integration manager class for nonlinear problems.
//

#ifndef included_algs_ImplicitIntegrator_C
#define included_algs_ImplicitIntegrator_C

#include "ImplicitIntegrator.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "tbox/SAMRAI_MPI.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define ALGS_IMPLICIT_INTEGRATOR_VERSION (1)

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#include "ImplicitIntegrator.I"
#endif
namespace SAMRAI {
    namespace algs {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for ImplicitIntegrator<DIM>.  The         *
* constructor sets default values for data members, then overrides      *
* them with values read from input or restart.  The destructor does     *
* nothing interesting.                                                  * 
*                                                                       *
*************************************************************************
*/

template<int DIM>  ImplicitIntegrator<DIM>::ImplicitIntegrator(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   ImplicitEquationStrategy<DIM>* implicit_equations,
   solv::NonlinearSolverStrategy<DIM>* nonlinear_solver,
   const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(implicit_equations != ((ImplicitEquationStrategy<DIM>*)NULL));
   TBOX_ASSERT(nonlinear_solver != ((solv::NonlinearSolverStrategy<DIM>*)NULL));
   TBOX_ASSERT(!hierarchy.isNull());
#endif

   d_object_name         = object_name;
   d_implicit_equations  = implicit_equations;
   d_nonlinear_solver    = nonlinear_solver;
   d_patch_hierarchy     = hierarchy;

   d_solution_vector.setNull();

   d_initial_time =
   d_final_time =
   d_current_time =
   d_current_dt =
   d_old_dt = tbox::MathUtilities<double>::getSignalingNaN();

   d_integrator_step      = 0;
   d_max_integrator_steps = 0;

   d_finest_level = -1;

   /*
    * Initialize object with data read from input and restart databases.
    */

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart ) {
      getFromRestart();
   }

   getFromInput(input_db);

   d_current_dt    = d_implicit_equations->getInitialDt();
   d_old_dt        = 0.0;
   d_current_time  = d_initial_time;

}

template<int DIM>  ImplicitIntegrator<DIM>::~ImplicitIntegrator()
{
}

/*
*************************************************************************
*                                                                       *
* Initialize integrator and nonlinear solver:                           *
*                                                                       *
* (1) Create vector containing solution state advanced in time.         *
*                                                                       *
* (2) Equation class registers data components with solution vector.    *
*                                                                       *
* (3) Initialize nonlinear solver.                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void ImplicitIntegrator<DIM>::initialize()
{
   d_finest_level = d_patch_hierarchy->getFinestLevelNumber();

   d_solution_vector = new solv::SAMRAIVectorReal<DIM,double>("solution_vector",
                                                          d_patch_hierarchy,
                                                          0, d_finest_level);

   d_implicit_equations->setupSolutionVector(d_solution_vector);

   if (d_solution_vector->getNumberOfComponents() == 0) {
      TBOX_ERROR("Solution vector has zero components.");
   }

   d_nonlinear_solver->initialize(d_solution_vector);
}

/*
*************************************************************************
*                                                                       *
* Integrate solution through given time increment:                      *
*                                                                       *
* (1) If number of levels in hierarchy has changed since last advance   *
*     due to regridding, the range of levels in the vectors is reset.   * 
*                                                                       *
* (2) Construct initial guess at new solution by extrapolation.         *
*                                                                       *
* (3) Call the equation advance set up routine.                         *
*                                                                       *
* (4) Compute the new solution using the nonlinear solver.              *
*                                                                       *
* (5) Return integer return code define by nonlinear solver.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> int ImplicitIntegrator<DIM>::advanceSolution(const double dt,
                                              const bool first_step)
{
   int retcode = tbox::MathUtilities<int>::getMax();

   if (stepsRemaining() && (d_current_time < d_final_time)) {

      d_current_dt = dt;

      const int finest_now = d_patch_hierarchy->getFinestLevelNumber();

      if (first_step && (finest_now != d_finest_level)) {
   
         d_finest_level = finest_now;

         d_solution_vector->resetLevels(0, d_finest_level); 

         d_nonlinear_solver->initialize(d_solution_vector);

      }

      d_implicit_equations->setInitialGuess(first_step,
                                            d_current_time,
                                            d_current_dt,
                                            d_old_dt);

      retcode = d_nonlinear_solver->solve();

   }
   
   return(retcode);
}

/*
*************************************************************************
*                                                                       *
* Get next dt from user-supplied equation class.  Timestep selection    *
* is generally based on whether the nonlinear solution iteration        *
* converged and, if so, whether the solution meets some user-defined    *
* criteria.  It is assumed that, before this routine is called, the     *
* routine checkNewSolution() has been called.  The boolean argument     *
* is the return value from that call.  The integer argument is          *
* that which is returned by the particular nonlinear solver package     *
* that generated the solution.                                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> double ImplicitIntegrator<DIM>::getNextDt(
   const bool good_solution,
   const int solver_retcode)
{
   double dt_next = d_implicit_equations->getNextDt(good_solution,
                                                    solver_retcode);

   double global_dt_next = tbox::SAMRAI_MPI::minReduction(dt_next);

   global_dt_next = 
      tbox::MathUtilities<double>::Min(global_dt_next,
                                       d_final_time - d_current_time); 

   return(global_dt_next);
}

/*
*************************************************************************
*                                                                       *
* Check whether time advanced solution is acceptable.  Note that the    *
* user-supplied solution checking routine must interpret the integer    *
* return code generated by the nonlinear solver correctly.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool ImplicitIntegrator<DIM>::checkNewSolution(
   const int solver_retcode) const
{
   bool good_solution = 
      d_implicit_equations->checkNewSolution(solver_retcode);

   int good = (good_solution ? 1 : 0);
   int global_good = tbox::SAMRAI_MPI::minReduction(good);

   return (global_good == 0 ? false : true);
}

/*
*************************************************************************
*                                                                       *
* Assuming an acceptable time advanced solution is found, update        *
* solution quantities and time information state of integrator.         *
* Return the current simulation time.                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> double ImplicitIntegrator<DIM>::updateSolution()
{
   d_current_time += d_current_dt;
   d_old_dt        = d_current_dt;
   d_integrator_step++;

   d_implicit_equations->updateSolution(d_current_time);

   return(d_current_time);
}

/*
*************************************************************************
*                                                                       *
* If simulation is not from restart, read data from input database.     *
* Otherwise, override restart values for a subset of the data members   *
* with those found in input.                                            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void ImplicitIntegrator<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif
 
  if ( tbox::RestartManager::getManager()->isFromRestart() ) {

     if (db->keyExists("final_time")) {
        d_final_time = db->getDouble("final_time");
        if (d_final_time < d_initial_time) {
           TBOX_ERROR(d_object_name << " -- Error in input data "
                      << "final_time < initial_time.");
        }
     }

     if (db->keyExists("max_integrator_steps")) {
        d_max_integrator_steps = db->getInteger("max_integrator_steps");
        if (d_max_integrator_steps < 0) {
           TBOX_ERROR(d_object_name << " -- Error in input data "
                      << "max_integrator_steps < 0.");
        } else {
           if (d_max_integrator_steps < d_integrator_step) {
              TBOX_ERROR(d_object_name << " -- Error in input data "
                         << "max_integrator_steps < current integrator step."); 
           }
        }
  
     }

  } else {

     if (db->keyExists("initial_time")) {
        d_initial_time = db->getDouble("initial_time");
     } else {
        TBOX_ERROR(d_object_name << " -- Key data `initial_time'"
                                 << " missing in input.");
     }

     if (db->keyExists("final_time")) {
        d_final_time = db->getDouble("final_time");
        if (d_final_time < d_initial_time) {
           TBOX_ERROR(d_object_name << " -- Error in input data "
                      << "final_time < initial_time.");
        }
     } else {
        TBOX_ERROR(d_object_name << " -- Key data `final_time'"
                                 << " missing in input.");
     }

     if (db->keyExists("max_integrator_steps")) {
        d_max_integrator_steps = db->getInteger("max_integrator_steps");
        if (d_max_integrator_steps < 0) {
           TBOX_ERROR(d_object_name << " -- Error in input data "
                      << "max_integrator_steps < 0.");
        }
     } else {
        TBOX_ERROR(d_object_name << " -- Key data `max_integrator_steps'"
                                 << " missing in input.");
     }

  }
   
}

/*
*************************************************************************
*                                                                       *
* Write out class version number and data members to database.          *
*                                                                       *
*************************************************************************
*/

template<int DIM> void ImplicitIntegrator<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   db->putInteger("ALGS_IMPLICIT_INTEGRATOR_VERSION",
                  ALGS_IMPLICIT_INTEGRATOR_VERSION);

   db->putDouble("d_initial_time", d_initial_time);
   db->putDouble("d_final_time", d_final_time);
   db->putDouble("d_current_time", d_current_time);
   db->putDouble("d_current_dt", d_current_dt);
   db->putDouble("d_old_dt", d_old_dt);

   db->putInteger("d_integrator_step", d_integrator_step);
   db->putInteger("d_max_integrator_steps", d_max_integrator_steps);

}

/*
*************************************************************************
*                                                                       *
* Check to make sure that the version number of the class is that same  *
* as the version number in the restart file.  If these values are equal *
* then read values for data members from the restart file.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void ImplicitIntegrator<DIM>::getFromRestart()
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

   int ver = db->getInteger("ALGS_IMPLICIT_INTEGRATOR_VERSION");
   if (ver != ALGS_IMPLICIT_INTEGRATOR_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version different "
              << "than class version.");
   }

   d_initial_time = db->getDouble("d_initial_time");
   d_final_time = db->getDouble("d_final_time");
   d_current_time = db->getDouble("d_current_time");
   d_current_dt = db->getDouble("d_current_dt");
   d_old_dt = db->getDouble("d_old_dt");

   d_integrator_step = db->getInteger("d_integrator_step");
   d_max_integrator_steps = db->getInteger("d_max_integrator_steps");

}

/*
*************************************************************************
*                                                                       *
* Print class data members to given output stream.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void ImplicitIntegrator<DIM>::printClassData(std::ostream& os) const
{
   os << "\nImplicitIntegrator<DIM>::printClassData..." << std::endl;
   os << "ImplicitIntegrator<DIM>: this = "
      << (ImplicitIntegrator<DIM>*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl; 
   os << "d_implicit_equations = " 
      << (ImplicitEquationStrategy<DIM>*)d_implicit_equations << std::endl;
   os << "d_nonlinear_solver = " 
      << (solv::NonlinearSolverStrategy<DIM>*)d_nonlinear_solver << std::endl;
   os << "d_patch_hierarchy = " 
      << (hier::PatchHierarchy<DIM>*)d_patch_hierarchy << std::endl;
   os << "d_solution_vector = " 
      << (solv::SAMRAIVectorReal<DIM,double>*)d_solution_vector << std::endl;

   os << "d_finest_level = " << d_finest_level << std::endl;
   os << "d_initial_time = " << d_initial_time << std::endl;
   os << "d_final_time = " << d_final_time << std::endl;
   os << "d_current_time = " << d_current_time << std::endl;
   os << "d_current_dt = " << d_current_dt << std::endl;
   os << "d_old_dt = " << d_old_dt << std::endl;
   os << "d_integrator_step = " << d_integrator_step << std::endl;
   os << "d_max_integrator_steps = " << d_max_integrator_steps << std::endl;
}

}
}
#endif
