//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/method_of_lines/MethodOfLinesIntegrator.h $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Basic method-of-lines time integration algorithm
//

#ifndef included_algs_MethodOfLinesIntegrator
#define included_algs_MethodOfLinesIntegrator

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h" 
#include "CoarsenAlgorithm.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "tbox/Database.h" 
#include "Geometry.h"
#include "GriddingAlgorithm.h"
#include "tbox/List.h" 
#include "MethodOfLinesPatchStrategy.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h" 
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "tbox/Serializable.h" 
#ifndef included_String
#include <string>
#define included_String
#endif 
#include "StandardTagAndInitStrategy.h"
#include "Variable.h"
#include "VariableContext.h"


#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace algs {

/*!
 * \brief Class MethodOfLinesIntegrator<DIM> implements a spatially
 * adaptive version of the Strong Stability Preserving (SSP) Runge-Kutta 
 * time integration algorithm.
 *
 * The original non-adaptive version of the algorithm is described in
 * S. Gottlieb, C.W. Shu, E. Tadmor, SIAM Review, Vol. 43, No. 1, pp. 89-112.
 * The advanceHierarchy() method integrates all levels of an AMR hierarchy 
 * through a specified timestep.  See this method for details of the 
 * time-stepping process.  Application-specific numerical routines that 
 * are necessary for these operations are provided by the 
 * MethodOfLinesPatchStrategy<DIM> data member.  The collaboration between 
 * this class and the patch strategy follows the the Strategy design pattern.
 * A concrete patch strategy object is derived from the base class to 
 * provide those routines for a specific problem.
 * 
 * This class is derived from the mesh::StandardTagAndInitStrategy<DIM> abstract
 * base class which defines an interface for routines required by the
 * dynamic adaptive mesh refinement routines in the mesh::GriddingAlgorithm<DIM>
 * class.  This collaboration also follows the Strategy design pattern.
 *
 * Initialization of an MethodOfLinesIntegrator<DIM> object is performed 
 * by first setting default values, then reading from input.  All input 
 * values may override values read from restart.  Data read from input is 
 * summarized as follows:
 *
 * Required input keys and data types: NONE
 *
 * Optional input keys, data types, and defaults:
 * 
 *    - \b    order   
 *       integer value specifying order of Runge-Kutta scheme.  If no input
 *       value is given, third order (i.e. order = 3) is used.
 *
 *    - \b    alpha_1   
 *    - \b    alpha_2
 *    - \b    beta 
 *       arrays of double values (length = order) specifying the coeffients
 *       used in the multi-step Strong Stability Preserving (SSP) Runge-Kutta 
 *       algorithm.  If no input is supplied, the default alpha_1, alpha_2,
 *       and beta values are automatically set to correspond to the 
 *       specified order.  
 *
 *
 * The following represents a sample input entry: 
 * 
 * \verbatim
 *  MethodOfLinesIntegrator{
 *     order                 = 3         
 *     alpha_1               = 1., 0.75, 0.33333
 *     alpha_2               = 0., 0.25, 0.66666
 *     beta                  = 1., 0.25, 0.66666
 *  }
 *  \endverbatim
 *
 * @see mesh::StandardTagAndInitStrategy 
 */

template<int DIM> class MethodOfLinesIntegrator : 
   public tbox::Serializable,
   public mesh::StandardTagAndInitStrategy<DIM>
{
public:
   /*!
    * Enumerated type for the different categories of variable 
    * quantities allowed by the method of lines integration algorithm.
    * See registerVariable(...) function for more details.
    * 
    * - \b SOLN           {Solution quantity for time-dependent ODE problem
    *                      solved by RK time-stepping algorithm.}
    * - \b RHS            {Right-hand-side of ODE problem solved;
    *                      i.e., du/dt = RHS.}
    * 


    */
   enum MOL_VAR_TYPE{ SOLN = 0,
                      RHS = 1 };

   /*!
    * The constructor for MethodOfLinesIntegrator<DIM> configures the method
    * of lines integration algorithm with the concrete patch strategy object 
    * (containing problem-specific numerical routines) and initializes 
    * integration algorithm parameters provided in the specified input 
    * database and in the restart database corresponding to the 
    * specified object_name.  The constructor also registers this object
    * for restart using the specified object name when the boolean
    * argument is true.  Whether object will write its state to
    * restart files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * When assertion checking is active, passing in any null pointer
    * or an empty string will result in an unrecoverable exception.
    */
   MethodOfLinesIntegrator(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      MethodOfLinesPatchStrategy<DIM>* patch_strategy,
      bool register_for_restart = true);

   /*!
    * The destructor for MethodOfLinesIntegrator<DIM> unregisters
    * the integrator object with the restart manager when so registered.
    */
   virtual ~MethodOfLinesIntegrator<DIM>();

   /*!
    * Initialize integrator by setting the number of time levels
    * of data needed based on specifications of the gridding algorithm.
    *
    * This routine also invokes variable registration in the patch strategy.
    */
   void initializeIntegrator(
      tbox::Pointer< mesh::GriddingAlgorithm<DIM> > gridding_alg);

   /*!
    * Return a suitable time increment over which to integrate the ODE
    * problem.  A minimum is taken over the increment computed on
    * each patch in the hierarchy.
    */
   double getTimestep(
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const double time ) const;

   /*!
    * Advance the solution through the specified dt, which is assumed
    * for the problem and state of the solution.  Advances all patches
    * in the hierarchy passed in.
    */
   void advanceHierarchy(const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
                         const double time,
                         const double dt);

   /*!
    * Register variable quantity defined in the patch strategy with the
    * method of lines integrator which manipulates its storage.
    */
   void registerVariable(const tbox::Pointer< hier::Variable<DIM> > variable,
                         const hier::IntVector<DIM>& ghosts,
                         const MOL_VAR_TYPE m_v_type,
                         const tbox::Pointer< xfer::Geometry<DIM> >& transfer_geom,
                         const std::string &coarsen_name = std::string(),
                         const std::string &refine_name = std::string());

   /*!
    * Print all data members of MethodOfLinesIntegrator<DIM> object.
    */
   virtual void printClassData(std::ostream& os) const;

   /*!
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.  The old_level pointer corresponds to
    * the level that resided in the hierarchy before the level with the
    * specified number was introduced.  If the pointer is NULL, there was
    * no level in the hierarchy prior to the call and the level data is set
    * based on the user routines and the simulation time.  Otherwise, the 
    * specified level replaces the old level and the new level receives data 
    * from the old level appropriately before it is destroyed.
    *
    * Typically, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the data is to be set, the level number must
    * match that of the old level, if non-NULL.  If the old level is
    * non-NULL, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.  Then, user-supplied patch routines are called to further
    * initialize the data if needed.  The boolean argument after_regrid
    * is passed into the user's routines.
    *
    * The boolean argument initial_time indicates whether the integration
    * time corresponds to the initial simulation time.  If true, the level
    * should be initialized with initial simulation values.  Otherwise, it
    * should be assumed that the simulation time is at some point after the
    * start of the simulation.  This information is provided since the
    * initialization of the data on a patch may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether 
    * the level is the finest allowable level in the hierarchy.
    * 
    * Note: This function is overloaded from the base class 
    *       mesh::StandardTagAndInitStrategy<DIM>.
    */
   void
   initializeLevelData(const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                       const int level_number,
                       const double init_time,
                       const bool can_be_refined,
                       const bool initial_time,
                       const tbox::Pointer< hier::BasePatchLevel<DIM> > old_level =
                             tbox::Pointer< hier::BasePatchLevel<DIM> >(NULL), 
                       const bool allocate_data = true);

   /*!
    * Reset cached communication schedules after the hierarchy has changed
    * (due to regridding, for example) and the data has been initialized on
    * the new levels.  The intent is that the cost of data movement on the
    * hierarchy will be amortized across multiple communication cycles,
    * if possible.  Note, that whenever this routine is called, communication
    * schedules are updated for every level finer than and including that
    * indexed by coarsest_level.
    * 
    * Note: This function is overloaded from the base class 
    *       mesh::StandardTagAndInitStrategy<DIM>.
    */
   void resetHierarchyConfiguration(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    * Set integer tags to "one" on the given level where refinement
    * of that level should occur using the user-supplied gradient detector.
    * The boolean argument initial_time is true when the level is being 
    * subject to error estimation at initialization time.  If it is false, 
    * the error estimation process is being invoked at some later time
    * after the AMR hierarchy was initially constructed.  The boolean argument
    * uses_richardson_extrapolation_too is true when Richardson
    * extrapolation error estimation is used in addition to the gradient
    * detector, and false otherwise.  This argument helps the user to
    * manage multiple regridding criteria.  This information
    * is passed along to the user's patch data tagging routines since the
    * application of the error estimator may be different in each of those
    * circumstances.
    * 
    * Note: This function is overloaded from the base class 
    *       mesh::StandardTagAndInitStrategy<DIM>.
    */
   virtual void
   applyGradientDetector(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /*!
    * Writes object state out to the given database.
    *  
    * When assertion checking is enabled, the database pointer must be non-null.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

private:

   /*
    * Copy all solution data from current context to scratch context.
    */
   void copyCurrentToScratch(
      const tbox::Pointer< hier::PatchLevel<DIM> > level) const;

   /*
    * Copy all solution data from scratch context to current context.
    */
   void copyScratchToCurrent(
      const tbox::Pointer< hier::PatchLevel<DIM> > level) const;

   /*
    * Reads in parameters from the input database.  All
    * values from the input file take precedence over values from the
    * restart file.
    *
    * When assertion checking enabled: db must not be a non-NULL pointer.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     bool is_from_restart);

   /*
    * Read object state from the restart file and initialize class data
    * members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not
    *     match.
    * 
    *    -Data is missing from the restart database or is inconsistent.
    *
    */
   void getFromRestart();

   /*
    * The object name is used as a handle to the database stored in 
    * restart files and for error reporting purposes.  The boolean
    * is used to control restart file writing operations.
    */
   std::string d_object_name;
   bool d_registered_for_restart;

   /*
    * Order of the Runge-Kutta method, and array of alpha values used in
    * updating solution during multi-step process.
    */
   int d_order;
   tbox::Array<double> d_alpha_1;
   tbox::Array<double> d_alpha_2;
   tbox::Array<double> d_beta;

   /*
    * A pointer to the method of lines patch model that will perform
    * the patch-based numerical operations.
    */
   MethodOfLinesPatchStrategy<DIM>* d_patch_strategy;

   /*
    * The communication algorithms and schedules are created and
    * maintained to manage inter-patch communication during AMR integration.
    * The algorithms are created in the class constructor.  They are 
    * initialized when variables are registered with the integrator.
    */

   /*
    * The "advance" schedule is used prior to advancing a level and
    * prior to computing dt at initialization.  It must be reset each
    * time a level is regridded.  All ghosts are filled with current
    * data at specified time.
    */
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_bdry_fill_advance;
   tbox::Array< tbox::Pointer< xfer::RefineSchedule<DIM> > > d_bdry_sched_advance;

   /*
    * Algorithm for transferring data from coarse patch to fine patch 
    * after a regrid.
    */
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_fill_after_regrid;

   /*
    * Algorithm for copying data from current context to scratch context,
    * on the same level.
    */
   tbox::Pointer< xfer::RefineAlgorithm<DIM> > d_fill_before_tagging;

   /*
    * Algorithm and communication schedule for transferring data from 
    * fine to coarse grid.
    */
   tbox::Pointer< xfer::CoarsenAlgorithm<DIM> > d_coarsen_algorithm;
   tbox::Array< tbox::Pointer< xfer::CoarsenSchedule<DIM> > > d_coarsen_schedule;

   /*
    * This algorithm has two variable contexts.  The current context is the
    * solution state at the current simulation time.  These data values do
    * not need ghost cells.  The scratch context is the temporary solution
    * state during the time integration process.  These variables will require
    * ghost cell widths that depend on the spatial discretization.
    */
   tbox::Pointer<hier::VariableContext> d_current;
   tbox::Pointer<hier::VariableContext> d_scratch;

   tbox::List< tbox::Pointer< hier::Variable<DIM> > > d_soln_variables;
   tbox::List< tbox::Pointer< hier::Variable<DIM> > > d_rhs_variables;

   /*
    * The component selectors for current and scratch are used by the
    * algorithm to collectively allocate/deallocate data for variables.
    */
   hier::ComponentSelector d_current_data;
   hier::ComponentSelector d_scratch_data;
   hier::ComponentSelector d_rhs_data;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MethodOfLinesIntegrator.C"
#endif
