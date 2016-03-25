//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/gridding/BaseGriddingAlgorithm.h $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2003 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1818 $
// Modified:    $LastChangedDate: 2007-12-20 15:50:44 -0800 (Thu, 20 Dec 2007) $
// Description: AMR hierarchy generation and regridding routines.
//

#ifndef included_mesh_BaseGriddingAlgorithm
#define included_mesh_BaseGriddingAlgorithm

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_BasePatchHierarchy
#include "BasePatchHierarchy.h"
#endif
#ifndef included_mesh_TagAndInitializeStrategy
#include "TagAndInitializeStrategy.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif

namespace SAMRAI {
   namespace mesh {

/*!
 * @brief Virtual base class providing interface for gridding algorithm.
 *
 * Class BaseGriddingAlgorithm in a virtual base class that provides
 * an abstract interface for a gridding algorithm.  This allows higher-level
 * classes in SAMRAI and in applications that use SAMRAI to interface
 * with a gridding algorithm in an abstract manner without knowing whether
 * it is a mesh::GriddingAlgorithm that handles gridding operations on a
 * rectangular domain or if it is an algorithm for gridding on, for example,
 * a multiblock domain.
 *
 * @see mesh::GriddingAlgorithm
 */

template <int DIM> class BaseGriddingAlgorithm
:
public tbox::Serializable
{
public:
   /*!
    * @brief Default constructor
    */
   BaseGriddingAlgorithm();

   /*!
    * @brief Virtual destructor for BaseGriddingAlgorithm<DIM>.
    */
   virtual ~BaseGriddingAlgorithm<DIM>();

   /*!
    * @brief Construct the coarsest level in a hierarchy (i.e., level 0).
    *
    * If level 0 does not already exist, then the domain specification is
    * checked against the constraints of the grid generation procedures.
    * The level gridding strategy data member defines these constraints.
    * Recall that the domain specification is maintained by the grid geometry
    * associated with the hierarchy.  Generally, an unrecoverable
    * exception will result if the constraints are not satisfied.
    * 
    * If level 0 already exists in the hierarchy, then the routine will 
    * generate a new level by re-applying the load balancing procedure to 
    * the existing level.  Data will be moved from the old level to the 
    * new level and the pre-existing level 0 will be discarded.  Note that 
    * this routine is different than the routine makeFinerLevel() below, 
    * which is used to construct levels 1 and finer.  In particular, this
    * routine does not select cells for refinement, whereas the other 
    * routine does.
    *
    * Important note: If assertion checking is turned on, then an 
    * unrecoverable assertion will result if either the patch hierarchy 
    * or its grid geometry is NULL.
    *
    * @param hierarchy The hierarchy on which coarse level is constructed.
    * @param level_time Simulation time when level is constructed
    */
   virtual void makeCoarsestLevel(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const double level_time) = 0;

   /*!
    * @brief This routine attempts to create a new level in a hierarchy 
    * finer than the finest level currently residing in the hierarchy.  
    *
    * It will select cells for refinement on the finest level and construct 
    * a new finest level, if necessary.  If no cells are selected for 
    * refinement, no new level will be added to the hierarchy.   The boolean 
    * argument initial_time indicates whether the routine is called at the 
    * initial simulation time.  If true, this routine is used to build 
    * individual levels during the construction of the AMR hierarchy at the 
    * initial simulation time.  If false, the routine is being used to add 
    * new levels to the hierarchy at some later point.  In either case, the 
    * time value is the current simulation time.  Note that this routine 
    * cannot be used to construct the coarsest level in the hierarchy
    * (i.e., level 0).  The routine makeCoarsestLevel() above must be used
    * for that purpose.
    *
    * The tag buffer indicates the number of cells by which cells selected
    * for refinement will be buffered before new finer level boxes are 
    * constructed.  The buffer is important to keep phenomena of interest 
    * on refined regions of the mesh until adaptive regridding occurs next.
    * Thus, the buffer size should take into account how the simulation may 
    * evolve before regridding occurs (e.g., number of timesteps taken).
    *
    * Important note: If assertion checking is activated, several checks
    * are applied to the functions arguments.  If any check is violated,
    * an unrecoverable assertion will result.  In particular, the hierarchy 
    * pointer must be non-NULL and the given level number must match that of 
    * the finest level currently residing in the hierarchy.  Also, the 
    * the tag buffer must be positive.
    *
    * @param hierarchy The hierarchy on which coarse level is constructed.
    * @param level_time Simulation time when level is constructed
    * @param initial_time Must be true if level_time is the initial time
    *                     of the simulation, false otherwise
    * @param tag_buffer Size of buffer around tagged cells that will be
    *                   covered by the fine level
    * @param regrid_start_time The simulation time when the regridding
    *                          operation began (this parameter is ignored
    *                          except when using Richardson extrapolation)
    */
   virtual void makeFinerLevel(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const double level_time,
      const bool initial_time,
      const int tag_buffer,
      const double regrid_start_time = 0.) = 0;

   /*!
    * @brief This routine attempts to reconfigure the patches on each level in
    * an AMR patch hierarchy which is finer than the specified level.
    *
    * The given level number is that of the coarsest level on which cells
    * will be selected for refinement.  In other words, that level is the 
    * finest level that will not be subject to a change in its patch 
    * configuration during the regridding process.  Generally, this routine
    * should be used to alter a pre-existing AMR patch hierarchy based on
    * the need to adapt the computational mesh around some phenomenon of 
    * interest.  The routine makeFinerLevel() above should be used to 
    * construct an initial hierarchy configuration or to add more than one
    * new level into the hierarchy.  Also, this routine will not reconfigure
    * the patches on level 0 (i.e., the coarsest in any hierarchy).  The
    * routine makeCoarsestLevel() above is provided for that purpose.
    *
    * Note that the current algorithm permits at most one new finest level 
    * to be added to the hierarchy with each invocation of the regridding 
    * process.  This constraint, though seemingly restrictive makes the 
    * process of maintaining properly nested levels much easier.
    *
    * The tag buffer array indicates the number of cells by which cells 
    * selected for refinement on a level will be buffered before new finer
    * level boxes are constructed.  The buffer is important to keep phenomena 
    * of interest on refined regions of the mesh until adaptive regridding 
    * occurs next.  Thus, the buffer size should take into account how the 
    * simulation may evolve before regridding occurs (e.g., number of 
    * timesteps taken on each level).
    *
    * The boolean argument is used for regridding in time-dependent 
    * problems.  When true, it indicates that the specified level is
    * the coarsest level to synchronize at the current regrid time 
    * before this regridding method is called.  This is a pretty
    * idiosyncratic argument but allows some flexibility in the way 
    * memory is managed during time-dependent regridding operations.
    *
    * Important note: If assertion checking is activated, several checks
    * are applied to the functions arguments.  If any check is violated,
    * an unrecoverable assertion will result.  In particular, the hierarchy
    * pointer must be non-NULL and the given level number must match that of
    * of some level in the hierarchy.  Also, the tag buffer array must 
    * contain a positive value for each level in the hierarchy.
    *
    * @param hierarchy The hierarchy on which coarse level is constructed.
    * @param level_number Coarsest level on which cells will be tagged for
    *                     refinement
    * @param regrid_time Simulaition time when regridding occurs
    * @param tag_buffer Size of buffer on each level around tagged cells
    *                   that will be covered by the next finer level
    * @param regrid_start_time The simulation time when the regridding
    *                          operation began on each level (this parameter is
    *                          ignored except when using Richardson
    *                          extrapolation)
    * @param level_is_coarsest_to_sync Level is the coarsest to sync
    */
   virtual void regridAllFinerLevels(
      tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const double regrid_time,
      const tbox::Array<int>& tag_buffer,
      tbox::Array<double> regrid_start_time = tbox::Array<double>(),
      const bool level_is_coarsest_to_sync = true) = 0;

   /*!
    * @brief Return true if error estimation process uses time integration;
    * otherwise, return false.
    */
   virtual bool errorEstimationUsesTimeIntegration() const = 0;

   /*!
    * @brief Return the error coarsen ratio.
    *
    * This is needed for cases where an error estimation scheme uses time
    * integration (e.g. Richardson extrapolation) to determine how many time
    * levels to maintain to properly apply the estimtion scheme. In general,
    * an even refine ratio (e.g. 2, 4, 8) will maintain two time levels, while
    * while an odd refine ratio (e.g. 3) will maintain three. 
    */
   virtual int getErrorCoarsenRatio() const = 0;

   /*!
    * @brief Return true if level associated with the specified level number
    * can be refined; i.e., the level number is less than that of the finest
    * level allowed in the hierarchy.  Otherwise, false is returned.
    */
   virtual bool levelCanBeRefined(const int level_number) const = 0;

   /*!
    * @brief Return pointer to level gridding strategy data member.
    */
   virtual
   tbox::Pointer< TagAndInitializeStrategy<DIM> > 
      getTagAndInitializeStrategy() const = 0;

   /*!
    * @brief Return maximum number of levels allowed in hierarchy.
    */
   virtual int getMaxLevels() const = 0;

   /*!
    * @brief Return const reference to ratio between specified level and next
    * coarser.
    */
   virtual const hier::IntVector<DIM>& 
      getRatioToCoarserLevel(const int level_number) const = 0;

   /*!
    * @brief Return efficiency tolerance for clustering tags on level.
    */
   virtual double getEfficiencyTolerance(const int level_number) const = 0;

   /*!
    * @brief Return combine efficiency for clustering tags on level.
    */
   virtual double getCombineEfficiency(const int level_number) const = 0;

   /*!
    * @brief Return proper nesting buffer width for level.
    */
   virtual int getProperNestingBuffer(const int level_number) const = 0;

   /*!
    * @brief Write object state out to the given database.
    * 
    * When assertion checking is active, the database pointer must be non-null.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db) = 0;

private:

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BaseGriddingAlgorithm.C"
#endif
