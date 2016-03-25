//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/time_refinement/TimeRefinementIntegrator.h $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Time integration manager for AMR with local time stepping.
//

#ifndef included_algs_TimeRefinementIntegrator
#define included_algs_TimeRefinementIntegrator

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_BasePatchHierarchy
#include "BasePatchHierarchy.h"
#endif
#ifndef included_mesh_BaseGriddingAlgorithm
#include "BaseGriddingAlgorithm.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif
#ifndef included_algs_TimeRefinementLevelStrategy
#include "TimeRefinementLevelStrategy.h"
#endif

namespace SAMRAI {
    namespace algs {

/**
 * Class TimeRefinementIntegrator<DIM> manages time integration over an
 * AMR patch hierarchy using local time refinement for finer hierarchy levels.
 * This class orchestrates hierarchy construction, data advancement and 
 * synchronization, and the dynamic grid refinement processes.  The basic 
 * ideas behind these algorithms are described in several sources on 
 * structured adaptive mesh refinement.  See Berger and Colella, J. Comp. 
 * Phys. (82)1:64-84, 1989 for an introduction to algorithm.  See Trangenstein,
 * SIAM J. Sci. Comput. 16(4):819-839, 1995, or Hornung, PhD thesis, Dept. 
 * of Mathematics, Duke University, 1994 for further discussion.  
 *
 * This class can be used in two different modes:  refined timestepping,
 * which divides they hierarchy's timestep into smaller timesteps on finer
 * levels, or synchronized timestepping, which advances all levels in a
 * hierarchy by the same timestep.  The mode that is used is determined
 * by querying the level integrator that is passed into the constructor.
 * One and only one mode can be used for each instantiation of
 * this class.
 *
 * The algorithm requires that integration steps on different levels are 
 * interleaved since the time increment used on each level is determined
 * by the spatial resolution of the mesh on that level (e.g., CFL condition).  
 * Generally, when using refined timestepping, coarser levels use larger time
 * increments than finer levels.  Thus, data must be synchronized between
 * levels and the dynamic regridding process must be coordinated with
 * the time stepping sequence.
 *
 * The routines in this class are implemented in a manner that is generic 
 * with respect to the details of the level integration and regridding 
 * methods, and the discrete equations being solved.  Thus, the class may
 * be employed for a variety of applications.  Upon construction, an object
 * of this class is configured with routines suitable for a given problem.
 * The TimeRefinementLevelStrategy<DIM> data member supplies routines
 * for advancing levels and synchronizing data on different levels during
 * time integration.  The mesh::BaseGriddingAlgorithm<DIM> data member provides
 * operations that construct and dynamically reconfigure the patch hierarchy.
 * The collaboration between this class and each of those objects follows 
 * the Strategy design pattern.
 *
 * Initialization begins by setting data on the coarsest AMR hierarchy
 * level.  Then, each successively finer level is created and initialized
 * after invoking the regridding procedures on the most recently initialized 
 * level.  This process is performed until either a maximum number of levels 
 * is reached or no further refinement is needed.
 *
 * Time integration is performed by invoking a recursive advance procedure on 
 * the coarsest AMR hierarchy level.  On a level, data is integrated to
 * a given point using a sequence of integration steps, where the size of 
 * each time increment depends on the problem being solved.  After each
 * step on a level, the next finer level (if it exists) is integrated to
 * the same time using a sequence of time steps appropriate for the level.
 * This class may dynamically adjust the time step sequence used on each 
 * level during the data advance process depending on requirements of the
 * integrator and information about stable time step size.  Dynamic 
 * mesh regridding is invoked during the integration process so that 
 * time integration, data synchronization, and mesh movement are coordinated
 * properly.
 *
 * An object of this class requires numerous parameters to be read from 
 * input.  Also, data must be written to and read from files for restart.
 * The input data are summarized as follows.
 *
 * Required input keys and data types:
 * 
 *    - \b    start_time   
 *        double value representing the start time for the simulation.
 * 
 *    - \b    end_time   
 *        double value representing the end time for the simulation.
 * 
 *    - \b    grow_dt   
 *        double value representing the maximum factor by which each
 *        succesive time increment may grow (typically >= 1.0).
 * 
 *    - \b    max_integrator_steps   
 *        integer value representing the maximum number of timesteps
 *        performed on the coarsest hierarchy level during the simulation.
 *
 *
 * Optional input keys, data types, and defaults:
 * 
 *    - \b    tag_buffer   
 *       array of integer values (one for each level that may be refined)
 *       representing the number of cells by which tagged cells are buffered 
 *       before clustering into boxes.  If no input is given, a default value 
 *       equal to the number of steps taken on the level before the next 
 *       regrid is used.
 *
 *
 * Note that the input values for end_time, grow_dt, max_integrator_step, 
 * and tag_buffer override values read in from restart.
 * 
 * A sample input file entry might look like:
 * 
 * \verbatim
 *
 *    start_time            = 0.e0      // initial simulation time
 *    end_time              = 10.e0     // final simulation time
 *    grow_dt               = 1.1e0     // growth factor for timesteps
 *    max_integrator_steps  = 50        // max number of simulation timesteps
 *    tag_buffer            = 1,1,1,1   // a max of 4 finer levels in hierarchy
 *
 * \endverbatim
 *
 * When running in synchronized timestepping mode, an additional input
 * key 'regrid_interval' can be added to specify the number of timesteps
 * between each regrid of the hierarchy.
 *
 * @see algs::TimeRefinementLevelStrategy
 * @see mesh::BaseGriddingAlgorithm
 */

template<int DIM> class TimeRefinementIntegrator 
:
   public virtual tbox::DescribedClass, 
   public tbox::Serializable
{
public:
   /**
    * The constructor for TimeRefinementIntegrator<DIM> initializes the
    * time stepping parameters needed to integrate the levels in the AMR
    * hierarchy.   Some data is set to default values; others are read 
    * from the specified input database and the restart database 
    * corresponding to the specified object_name.  Consult top of
    * this header file for further details.  The constructor also
    * registers this object for restart using the specified object name
    * when the boolean argument is true.  Whether object will write its state to
    * restart files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * Note that this object also invokes the variable creation and 
    * registration process in the level strategy.
    * 
    * If assertion checking is turned on, an unrecoverable assertion will
    * result if any of the input database, patch hierarchy,
    * level strategy, or regridding algorithm pointers is null.  Exceptions
    * may also be thrown if any checks for consistency between parameters 
    * in the gridding algorithm, level strategy, and this object fail.
    */
   TimeRefinementIntegrator(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      TimeRefinementLevelStrategy<DIM>* level_integrator, 
      tbox::Pointer< mesh::BaseGriddingAlgorithm<DIM> > gridding_algorithm,
      bool register_for_restart = true);

   /**
    * The destructor for TimeRefinementIntegrator<DIM> unregisters 
    * the integrator object with the restart manager when so registered.
    */
   ~TimeRefinementIntegrator<DIM>();

   /**
    * Set AMR patch hierarchy configuration and data at start of simulation.
    * If the run is begun from a restart file, the hierarchy and data 
    * are read from the hierarchy database.  Otherwise, the hierarchy 
    * and data are initialized by the gridding algorithm data member.
    * In this case, the coarsest level is constructed and initialized.
    * Then, error estimation is performed to determine if and where it 
    * should be refined.  Successively finer levels are created and 
    * initialized until the maximum allowable number of levels is achieved 
    * or no further refinement is needed.  The double return value is the 
    * time increment for the first data advance step on the coarsest 
    * hierarchy level (i.e., level 0).
    *
    * This function assumes that the hierarchy exists, but that it contains
    * no patch levels, when it is called.  On return from this function, the 
    * initial hierarchy configuration and simulation data is set properly for 
    * the advanceHierarchy() function to be called.  In particular, on each
    * level constructed only the data needed for initialization exists.
    * 
    * When assertion checking is active, the hierachy database pointer
    * must be non-null.
    */
   double initializeHierarchy();
 
   /**
    * Advance each level in the hierarchy through the given time increment 
    * and return an appropriate time increment for subsequent advances of the 
    * coarsest hierarchy level (level 0).  The boolean argument indicates 
    * whether the coarsest hierarchy level (i.e., level 0) should be load 
    * balanced before the levels are advanced.  In general, the problem
    * domain (determined by the union of patches on level 0) does not change 
    * once set.  However, the boolean flag here allows one to reconfigure
    * the patches on the coarsest level which constitute this union.  This
    * may be required depending on a dynamic change of the work load.  
    * By default, the level will not be subject to load balancing.
    *
    * This function assumes that all data on each level in the hierarchy
    * has been set and that only the data need for initialization exists 
    * on each level (as opposed to both current and new data, for example).
    * Upon return from this function, the simulation data on each hierarchy 
    * levels is advanced through the time increment dt.  In addition, data on 
    * all hierarchy levels has been synchronized so that it is consistent at 
    * the new simulation time (where this synchronization process is defined 
    * by the level strategy).  Thus, the data is set properly for any 
    * subsequent calls to this function.
    */
   double advanceHierarchy(
      const double dt, 
      const bool rebalance_coarsest = false);
 
   /**
    * Return true if the current step count for the level indicates 
    * that regridding should occur.  In particular, true is returned 
    * if both the level allows refinement and the step count is an 
    * integer multiple of the regrid step interval.  
    * Otherwise, false is returned.
    */
   bool atRegridPoint(const int level_number) const;
 
   /**
    * Return current integration time for coarsest hierarchy level.
    */
   double getIntegratorTime() const;
 
   /**
    * Return initial integration time.
    */
   double getStartTime() const;
 
   /**
    * Return final integration time.
    */
   double getEndTime() const;

   /**
    * Return integration step count for entire hierarchy 
    * (i.e., number of steps taken by the coarsest level).
    */
   int getIntegratorStep() const; 

   /**
    * Return maximum number of integration steps allowed for entire 
    * hierarchy (i.e., steps allowed on coarsest level).
    */
   int getMaxIntegratorSteps() const; 

   /**
    * Return true if any steps remain in current step sequence on level 
    * (i.e., before it will synchronize with some coarser level).
    * Return false otherwise. 
    */
   bool stepsRemaining(const int level_number) const; 

   /**
    * Return true if any integration steps remain, false otherwise. 
    */
   bool stepsRemaining() const; 

   /**
    * Return current time increment used to advance level. 
    */
   double getLevelDtActual(const int level_number) const; 

   /**
    * Return maximum time increment currently allowed on level. 
    */
   double getLevelDtMax(const int level_number) const; 

   /**
    * Return current simulation time for level.
    */
   double getLevelSimTime(const int level_number) const; 

   /**
    * Return step count for current integration sequence on level.
    */
   int getLevelStep(const int level_number) const; 

   /**
    * Return maximum number of time steps allowed on level in 
    * current integration step sequence.
    */
   int getLevelMaxSteps(const int level_number) const; 

   /**
    * Return const pointer to patch hierarchy managed by integrator. 
    */
   const tbox::Pointer< hier::BasePatchHierarchy<DIM> > getPatchHierarchy() const;

   /**
    * Return pointer to level integrator.
    */
   tbox::Pointer<TimeRefinementLevelStrategy<DIM> >
   getLevelIntegrator() const;

   /**
    * Return pointer to gridding algorithm object.
    */
   tbox::Pointer< mesh::BaseGriddingAlgorithm<DIM> > getGriddingAlgorithm() const;

   /**
    * Return true if current step on level is first in current step
    * sequence; otherwise return false.
    */
   bool firstLevelStep(const int level_number) const; 

   /**
    * Return true if current step on level is last in current step
    * sequence; otherwise return false.
    */
   bool lastLevelStep(const int level_number) const; 

   /**
    * set the regrid interval to a new value.  This may only be used
    * when using synchronized timestepping.
    */
   void setRegridInterval(const int regrid_interval);

   /**
    * Print data representation of this object to given output stream.
    */
   virtual void printClassData(std::ostream& os) const; 
  
   /**
    * Print time stepping data for a single level to given output stream.
    */
   void printDataForLevel(std::ostream& os, const int level_number) const;

   /**
    * Write object state out to the given database.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void putToDatabase( tbox::Pointer<tbox::Database> db);

private:
   /*
    * Initialize data on given level.  If the level can be refined, a problem-
    * dependent error estimation procedure is invoked to determine whether
    * further refinement is needed.  If a new level is created, this function 
    * is called recursively to initialize the next finest level.
    */
   void initializeRefinedTimesteppingLevelData(const int level_number);
   void initializeSynchronizedTimesteppingLevelData(const int level_number);

   /*
    * Advance the data on the level to the specified time using a dynamically 
    * adjusted sequence of time increments.  If any finer levels exist 
    * in the hierarchy when this function is called or are generated during
    * the regridding process, they will be advanced to the specified time
    * as well.  This function is recursive.  After a single timestep is 
    * performed on each level, all finer levels are advanced to the new
    * integration time before another timestep is taken on the original level.
    */
   void advanceRecursivelyForRefinedTimestepping(const int level_number,
                                                 const double end_time);

   double advanceForSynchronizedTimestepping(const double dt);

   /*
    * Determine the next stable time increment (dt) and adjust the step 
    * sequence if necessary for the given level.  The computed dt will
    * be less than or equal to the specified bound and the time remaining.
    * In adjusting the step sequence, an attempt is made to partition the
    * remaining time interval into a sequence of equal time increments.
    * However, the total number of time steps in the sequence for the 
    * level must satisfy any constraints placed on it by the regridding 
    * procedure.
    */
   bool findNextDtAndStepsRemaining(const int level_number,
                                    const double time_remaining,
                                    const double dt_bound);

   /*
    * Return true if the this level can be regridded at the current step
    * and the next coarser level can be regridded too.  Otherwise, 
    * false is returned.
    */
   bool coarserLevelRegridsToo(const int level_number) const;

   /*
    * Read input data from specified database and initialize class members.
    * The argument is_from_restart should be set to true if the simulation
    * is from restart.  Otherwise, it should be set to false.
    *
    * If the simulation is not from restart, read in start_time, end_time,
    * grow_dt, max_integrator_step, and possibly tag_buffer
    * from the database. 
    * 
    * If the simulation is from restart, then only read in end_time,
    * grow_dt, max_integrator_step and tag_buffer if they are 
    * found in the input database.
    *
    * When assertion checking is active, the databse pointer must be non-null. 
    */
   virtual void getFromInput(tbox::Pointer<tbox::Database> db,
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
    */
   virtual void getFromRestart();

   /* 
    * The object name is used as a handle to databases stored in 
    * restart files and for error reporting purposes.  The boolean
    * is used to control restart file writing operations.
    */
   std::string d_object_name;
   bool d_registered_for_restart;

   /*
    * Pointers to the patch hierarchy, level integration and gridding 
    * algorithm objects associated with this time integration object.
    * The level integrator defines operations for advancing data on
    * individual levels in the AMR patch hierarchy.  The gridding algorithm
    * provides grid generation and regridding routines for the AMR hierarchy. 
    */
   tbox::Pointer< hier::BasePatchHierarchy<DIM> > d_patch_hierarchy;
   TimeRefinementLevelStrategy<DIM>* d_refine_level_integrator;
   tbox::Pointer< mesh::BaseGriddingAlgorithm<DIM> > d_gridding_algorithm;

   /*
    */
   bool d_use_refined_timestepping;

   /*
    * Integrator data read from input or set at initialization.
    */
   double d_start_time;
   double d_end_time;
   double d_grow_dt;
   int d_max_integrator_steps;

   /*
    * The regrid interval indicates the number of integration steps taken
    * on a level between successive invocations of the regridding process
    * on that level.  In general, this class enforces the constraint that 
    * each synchronization time between two successive hierarchy levels 
    * will always be a potential regrid point for the coarser of the two 
    * levels.  Specifically, it sets the regrid interval for each level 
    * to be the greatest common divisor between the entries in the grid 
    * refinement ratio vector between the level and the next coarsest level 
    * in the hierarchy.  The regrid interval for level 0 is set equal to 
    * that for level 1.  In the future, users may be able to specify 
    * this value in the input file.
    */
   tbox::Array<int> d_regrid_interval;

   /*
    * The tag buffer indicates the number of cells on each level by which
    * tagged cells will be buffered after they have selected for refinement.
    * These values are passed into the gridding algorithm routines during
    * hierarchy construction and regridding.  The tag buffer helps to 
    * guarantee that refined cells near important features in the solution 
    * will remain refined until the level is regridded next.  
    *
    * Important note: these values may be specified in the input file.
    * If not, default values are set based on the regrid intervals.  
    * However, if the user attempts to specify these values, care must
    * be taken to assure that improper tag buffering will not degrade the
    * calculation.
    */
   tbox::Array<int> d_tag_buffer;

   /*
    * Integrator data that evolves during time integration and maintains
    * the state of the timestep sequence over the levels in the AMR hierarchy. 
    */
   double d_integrator_time; 
   int d_integrator_step;
   bool d_just_regridded;
   int  d_last_finest_level;
   tbox::Array<double> d_level_old_old_time;
   tbox::Array<double> d_level_old_time;
   tbox::Array<double> d_level_sim_time;
   tbox::Array<double> d_dt_max_level;
   tbox::Array<double> d_dt_actual_level;
   tbox::Array<int> d_step_level;
   tbox::Array<int> d_max_steps_level;

   double d_dt;

   /*
    * tbox::Timer objects for performance measurement.
    */
   tbox::Pointer<tbox::Timer> t_initialize_hier;
   tbox::Pointer<tbox::Timer> t_advance_hier;

   // The following are not implemented:
   TimeRefinementIntegrator(const TimeRefinementIntegrator<DIM>&);
   void operator=(const TimeRefinementIntegrator<DIM>&);

};

}
}
#ifndef DEBUG_NO_INLINE
#include "TimeRefinementIntegrator.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "TimeRefinementIntegrator.C"
#endif
