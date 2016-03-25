//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/time_refinement/TimeRefinementLevelStrategy.h $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Interface to level routines for time-refinement integrator.
//
 
#ifndef included_algs_TimeRefinementLevelStrategy
#define included_algs_TimeRefinementLevelStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_BasePatchHierarchy
#include "BasePatchHierarchy.h"
#endif
#ifndef included_hier_BasePatchLevel
#include "BasePatchLevel.h"
#endif
#ifndef included_mesh_BaseGriddingAlgorithm
#include "BaseGriddingAlgorithm.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace algs {

/**
 * Class TimeRefinementLevelStrategy<DIM> is an abstract base class that
 * defines the interface to level integration and synchronization routines 
 * needed by the hierarchy integration class TimeRefinementIntegrator<DIM>.
 * In particular, this class insulates the hierarchy integrator from the 
 * routines that manipulate data on the hierarchy levels in a problem-specific 
 * fashion.  When the AMR hierarchy integration and regridding sequence 
 * provided by the class TimeRefinementIntegrator<DIM> are appropriate 
 * for some computational problem, a subclass of this base class can be 
 * used to provide the necessary operations to the hierarchy integrator. 
 * That is, a TimeRefinementIntegrator<DIM> object may be configured with 
 * a concrete implementation of this base class by passing the concrete 
 * object into the to the time refinement integrator constructor.
 *
 * Three functions in this class are defined as virtual rather than
 * pure virtual, and default implementations are provided.  This is
 * done to allow for an implementation of this strategy class that is used
 * by the time refinement integrator only for refined timestepping or only
 * for synchronized timestepping.  If an implementation is only used for
 * refined timestepping, getMaxLevelDt() and one version of
 * standardLevelSynchronization() must be overloaded in the concrete class.
 * If an implementation is only used for synchronized timestepping,
 * the other version of standardLevelSynchronization must be overloaded.
 * If an implementation supports both refined and synchronized timestepping,
 * all three virtual functions must be overloaded.
 *
 * @see algs::TimeRefinementIntegrator
 */

template<int DIM> class TimeRefinementLevelStrategy
: 
public virtual tbox::DescribedClass
{
public:
   /**
    * Default constructor for TimeRefinementLevelStrategy<DIM>. 
    */     
   TimeRefinementLevelStrategy();

   /**
    * Virtual destructor for TimeRefinementLevelStrategy<DIM>.
    */
   virtual ~TimeRefinementLevelStrategy<DIM>();

   /**
    * Initialize the state of the integrator that performs time 
    * advances on the levels of the hierarchy.  Typically, this
    * involves setting up information, such as communication algorithms,
    * to manage variable storage. The pointer to the gridding algorithm 
    * is provided so that the integrator may access information about 
    * regridding procedures or the structure the hierarchy, which is not 
    * yet created.
    */
   virtual void 
   initializeLevelIntegrator(
      tbox::Pointer< mesh::BaseGriddingAlgorithm<DIM> > gridding_alg) = 0;

   /**
    * Return appropriate time increment for given level in the patch 
    * hierarchy.  This routine is called during the initial 
    * generation of the AMR patch hierarchy and possibly during regridding
    * if time regridding uses a time advance.  It should be assumed that 
    * the only data that exists on the level when this routine is called 
    * is that which is needed to initialize the level.  The initial_time
    * boolean flag is true if this routine is called at the initial
    * simulation time (i.e., when hierarchy is generated for first time);
    * otherwise (e.g., at an advance step) it is false.
    *
    * The recompute_dt option specifies whether to compute 
    * the timestep using the current level data or to return the value
    * stored by the time integrator. The default true setting means
    * the timestep will be computed if no value is supplied.  
    */
   virtual double
   getLevelDt(const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
              const double dt_time,
              const bool initial_time) = 0;

   /**
    * Return the maximum allowable time increment for level in the hierarchy
    * with the given level number.  The point of this routine is to determine
    * the increment for that level based on the current time increment used
    * on the next coarser level.  The coarse dt is the current time step 
    * size for the next coarser level.  The ratio is the mesh refinement 
    * ratio between the two levels.
    *
    * This is defined as a virtual function, and it is only used by the
    * time refinement integrator for refined timestepping.  It should be
    * overloaded in any concrete implementation of this class that is
    * used for refined timestepping, but need not be overloaded in
    * an implementation that is used only for synchronized timestepping.
    */
   virtual double
   getMaxFinerLevelDt(const int finer_level_number,
                      const double coarse_dt,
                      const hier::IntVector<DIM>& ratio) {

      NULL_USE(finer_level_number);
      NULL_USE(coarse_dt);
      NULL_USE(ratio);

      return (0.0);
   }

   /**
    * Advance data on all patches on specified patch level from current time
    * (current_time) to new time (new_time).   The boolean value first_step
    * indicates whether the advance step is the first in a time sequence
    * on the level.  The boolean value last_step indicates whether the 
    * advance step is the last in a time sequence on the level.   Usually,
    * the timestep sequence refers to the steps taken to advance the solution
    * through the time interval of the most recent advance on the next coarser
    * level, if such a level exists.  The boolean regrid_advance is false
    * when the advance is part of the actual hierarchy integration process 
    * and true when the advance is called during time-dependent regridding
    * (e.g., when using Richardson extrapolation).  The default value is false.
    * The last boolean argument is true when the level is in the hierarchy,
    * and false otherwise.  The default value is true.  Usually, this value
    * is false only during time-dependent regridding operations performed
    * on some temporary level; thus, a schedule must be generated for the 
    * level before the advance can occur, for example.
    * 
    * When this function is called, the level data required to begin the 
    * advance must be allocated and be defined appropriately.  Typically, 
    * this is equivalent to what is needed to initialize a new level after
    * regridding.  Upon exiting this routine, both current and new data may 
    * exist on the level.  This data is needed until level synchronization 
    * occurs, in general. Current and new data may be reset by calling 
    * the member function resetTimeDependentData(). 
    * 
    * This routine is called from two different points within the 
    * TimeRefinementIntegrator<DIM> class: during the regular time
    * advance sequence, and at the initial simulation time.  The second 
    * call is made to advance the solution on a coarser level ahead in time to 
    * provide time-dependent boundary values for some finer level when 
    * time-dependent regridding is used.  In the first case, the values of
    * the boolean flags are:
    * 


    *    - \b  first_step 
    *        = true for first step in level time step sequence; else, false. 
    *    - \b  last_step 
    *        = true for last step in level time step sequence; else, false. 
    *    - \b  regrid_advance 
    *        = false. 
    * 


    * In the second case, the values of the boolean flags are:
    * 


    *    - \b  first_step 
    *        = true.
    *    - \b  last_step 
    *        = false.
    *    - \b  regrid_advance 
    *        = true. 
    * 

 
    */ 
   virtual double 
   advanceLevel(const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
                const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                const double current_time, 
                const double new_time, 
                const bool first_step,
                const bool last_step,
                const bool regrid_advance = false) = 0;

   /**
    * Synchronize data on specified patch levels in AMR hierarchy at the
    * given synchronization time.  The array of time values provides the
    * previous integration time for each level involved in the synchronization.
    * In other words, (sync_time - old_times[ln]) is the most recent time
    * increment used to advance data on level ln.  These times are used when
    * the synchronization process requires re-integration of the data.  Note
    * that other synchronization routines are defined below for other points
    * in the hierarchy integration sequence.  
    *
    * When this routine is called, both current and new data may exist on each
    * level involved in the synchronization.  The new data on each level
    * corresponds to the synchronization time.   Each entry in the array
    * of time values specifies the time to which the current data on each 
    * level corresponds.  It is assumed that this routine will reset the 
    * synchronized data on each level so that only the current data will 
    * exist on each level when done.
    *
    * Note that this routine is distinct from the synchronizeNewLevels()
    * function below.  This routine is used to synchronize levels during
    * the time integration process.  The other routine is used to synchronize
    * new levels in the hierarchy, either at initialization time or after
    * regridding.
    *
    * This function must be overloaded in a concrete implementation
    * that is used by the time refinement integrator for refined
    * timestepping
    */
   virtual void 
   standardLevelSynchronization(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy, 
      const int coarsest_level, 
      const int finest_level, 
      const double sync_time, 
      const tbox::Array<double>& old_times) {

      NULL_USE(hierarchy);
      NULL_USE(coarsest_level);
      NULL_USE(finest_level);
      NULL_USE(sync_time);
      NULL_USE(old_times);
   }

   /**
    * This version of standardLevelSynchronization must be overloaded
    * for implementations that support synchronized timestepping.
    * The interface is the same as the other version of this function,
    * except for the old_time argument.  In synchronized timestepping
    * all levels advance the same timestep, so only a single old time
    * is required, rather than an array.
    */
   virtual void
   standardLevelSynchronization(
      const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int coarsest_level,
      const int finest_level,
      const double sync_time,
      const double old_time) {
     
      NULL_USE(hierarchy);
      NULL_USE(coarsest_level);
      NULL_USE(finest_level);
      NULL_USE(sync_time);
      NULL_USE(old_time);
   }

   /**
    * Synchronize specified levels after regridding has occurred or during
    * initial construction of the AMR patch hierarchy.  Note that this 
    * synchronization may be different than the standard time-dependent
    * synchronization (above in standardLevelSynchronization()) depending 
    * on the level integration algorithm.
    *
    * Before this routine is called, all time-dependent data on all levels
    * involved in the synchronization has been reset.  Thus, this routine
    * must only synchronize the current data on each level.  On return from
    * this function, only the current data on each level must be present.
    */
   virtual void
   synchronizeNewLevels(const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
                        const int coarsest_level,
                        const int finest_level,
                        const double sync_time,
                        const bool initial_time) = 0;

   /**
    * Reset time-dependent data storage for the specified patch level.  This
    * routine is called when the current level data is no longer needed
    * and it is appropriate to replace the current data with the new data 
    * on the level, if such data exists.
    */
   virtual void 
   resetTimeDependentData(const tbox::Pointer< hier::BasePatchLevel<DIM> > level,
                          const double new_time,
                          const bool can_be_refined) = 0;

   /**
    * Reset data on the patch level to state before time advance.
    * This is needed whenever the solution on a level is advanced beyond
    * the new time on a level.  For example, during time-dependent regridding 
    * (Richardson extrapolation) or initialization, it is necessary to 
    * do such an advance.  This routine is called to discard the new solution
    * data so that subsequent calls to advance are provided proper data at the 
    * correct time.
    */
   virtual void 
   resetDataToPreadvanceState(const tbox::Pointer< hier::BasePatchLevel<DIM> > level) = 0;

   /**
    * Return true if the implementation of this class is constructed
    * to use refined timestepping, and false otherwise.
    */
   virtual bool
   usingRefinedTimestepping() const = 0;

private:

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "TimeRefinementLevelStrategy.C"
#endif
