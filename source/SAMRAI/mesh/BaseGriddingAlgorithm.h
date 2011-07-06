/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines. 
 *
 ************************************************************************/

#ifndef included_mesh_BaseGriddingAlgorithm
#define included_mesh_BaseGriddingAlgorithm

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Serializable.h"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Virtual base class providing interface for gridding
 * algorithms.
 *
 * Class BaseGriddingAlgorithm in a virtual base class that provides
 * an abstract interface for a gridding algorithm.  This allows
 * higher-level classes in SAMRAI and in applications that use SAMRAI
 * to interface with a gridding algorithm in an abstract manner
 * without knowing whether it is a mesh::GriddingAlgorithm that
 * handles gridding operations on a rectangular domain or if it is an
 * algorithm for gridding on, for example, a multiblock domain.
 *
 * QUESTION: Should this be called GriddingAlgorithmStrategy?  There
 * are only pure virtual interfaces.
 *
 * Each BaseGriddingAlgorithm is constructed a PatchHierarchy.  All
 * hierarchy operations refer to this hierarchy.
 *
 * The BaseGriddingAlgorithm constructor requires a PatchHierarchy.
 * The implementation is then responsible for manipulating levels in
 * that hierarchy using the other interfaces defined here.
 *
 * @see mesh::GriddingAlgorithm
 */

class BaseGriddingAlgorithm:
   public tbox::Serializable
{
public:
   /*!
    * @brief Constructor
    *
    * @param hierarchy The hierarchy that this gridding algorithm
    * will operate on.  All hierarchy operations and data refers to
    * this hierarchy.
    */
   BaseGriddingAlgorithm(
      const tbox::Pointer<hier::PatchHierarchy> &hierarchy );

   /*!
    * @brief Virtual destructor for BaseGriddingAlgorithm.
    */
   virtual ~BaseGriddingAlgorithm();

   /*!
    * @brief Construct the coarsest level in the hierarchy (i.e., level 0).
    *
    * Level 0 should occupy the full domain.  If level 0 already
    * exists, it may be re-balanced.
    *
    * @param level_time Simulation time when level is constructed
    *
    * @param override_mapped_box_level MappedBoxLevel representing a
    *                                  decomposition of level zero of the
    *                                  hierarchy.
    */
   virtual void
   makeCoarsestLevel(const double level_time) = 0;

   /*
    * @brief Same as makeCoarsestLevel(const double) but allows user to
    * specify the configuration of level 0.
    *
    * @param level_time Simulation time when level is constructed
    *
    * @param override_mapped_box_level MappedBoxLevel representing a
    * decomposition of level zero of the hierarchy.  This should cover
    * all of the domain.
    */
   virtual void
   makeCoarsestLevel(
      const double level_time,
      const hier::MappedBoxLevel& override_mapped_box_level) = 0;

   /*!
    * @brief Attempts to create a new level in the hierarchy finer
    * than the finest level currently residing in the hierarchy.
    *
    * It should select cells for refinement on the finest level and
    * construct a new finest level, if necessary.  If no cells are
    * selected for refinement, no new level should be added to the
    * hierarchy.
    *
    * The boolean argument initial_time indicates whether the routine
    * is called at the initial simulation time.  If true, this routine
    * is being used to build individual levels during the construction of
    * the AMR hierarchy at the initial simulation time.  If false, the
    * routine is being used to add new levels to the hierarchy at some
    * later point.  In either case, the time value is the current
    * simulation time.  Note that this routine cannot be used to
    * construct the coarsest level in the hierarchy (i.e., level 0).
    * The routine makeCoarsestLevel() above must be used for that
    * purpose.
    *
    * The tag buffer indicates the number of cells by which cells
    * selected for refinement should be buffered before new finer
    * level boxes are constructed.  All tagged cells should be refined
    * except where refinement would violate proper nesting.  The
    * buffer is meant to keep phenomena of interest on refined regions
    * of the mesh until adaptive regridding occurs next.  Callers of
    * this method should take into account how the simulation may
    * evolve before regridding occurs (e.g., number of timesteps
    * taken) when calculating the tag_buffer.
    *
    * @param level_time Simulation time when level is constructed
    * @param initial_time Must be true if level_time is the initial time
    *                     of the simulation, false otherwise
    * @param tag_buffer Size of buffer around tagged cells that will be
    *                   covered by the fine level.  Must be non-negative.
    * @param regrid_start_time The simulation time when the regridding
    *                          operation began (this parameter is ignored
    *                          except when using Richardson extrapolation)
    */
   virtual void
   makeFinerLevel(
      const double level_time,
      const bool initial_time,
      const int tag_buffer,
      const double regrid_start_time = 0.0) = 0;

   /*!
    * @brief Attempt to regrid each level in the PatchHierarchy
    * that is finer than the specified level.
    *
    * The given level number is that of the coarsest level on which
    * cells will be selected for refinement.  In other words, that
    * level is the finest level that will not be subject to a change
    * in its patch configuration during the regridding process.
    * Generally, this routine should be used to alter the pre-existing
    * AMR patch hierarchy based on the need to adapt the computational
    * mesh around some phenomenon of interest.  The routine
    * makeFinerLevel() above should be used to construct the initial
    * hierarchy configuration or to add more than one new level into
    * the hierarchy.  Also, this routine will not reconfigure the
    * patches on level 0 (i.e., the coarsest in the hierarchy).  The
    * routine makeCoarsestLevel() above is provided for that purpose.
    *
    * The tag buffer array indicates the number of cells by which
    * cells selected for refinement on a level will be buffered before
    * new finer level boxes are constructed.  The buffer is important
    * to keep phenomena of interest on refined regions of the mesh
    * until adaptive regridding occurs next.  Thus, the buffer size
    * should take into account how the simulation may evolve before
    * regridding occurs (e.g., number of timesteps taken on each
    * level).  Tag buffers must be non-negative.
    *
    * The boolean level_is_coarsest_to_sync is used for regridding in
    * time-dependent problems.  When true, it indicates that the
    * specified level is the coarsest level to synchronize at the
    * current regrid time before this regridding method is called.
    * This is a pretty idiosyncratic argument but allows some
    * flexibility in the way memory is managed during time-dependent
    * regridding operations.
    *
    * @param level_number Coarsest level on which cells will be tagged for
    *                     refinement
    * @param regrid_time Simulation time when regridding occurs
    * @param tag_buffer Size of buffer on each level around tagged cells
    *                   that will be covered by the next finer level
    * @param regrid_start_time The simulation time when the regridding
    *                          operation began on each level (this parameter is
    *                          ignored except when using Richardson
    *                          extrapolation)
    * @param level_is_coarsest_to_sync Level is the coarsest to sync
    */
   virtual void
   regridAllFinerLevels(
      const int level_number,
      const double regrid_time,
      const tbox::Array<int>& tag_buffer,
      tbox::Array<double> regrid_start_time = tbox::Array<double>(),
      const bool level_is_coarsest_to_sync = true) = 0;

   /*!
    * @brief Return true if error estimation process uses time integration;
    * otherwise, return false.
    *
    * QUESTION: Will this ever need to be anything other than
    * getTagAndInitializeStrategy()->usesTimeIntegration() ?  If not,
    * do we really need this interface?
    */
   virtual bool
   errorEstimationUsesTimeIntegration() const = 0;

   /*!
    * @brief Return pointer to level gridding strategy data member.
    */
   virtual
   tbox::Pointer<TagAndInitializeStrategy>
   getTagAndInitializeStrategy() const = 0;

   /*!
    * @brief Return efficiency tolerance for clustering tags on level.
    */
   virtual double
   getEfficiencyTolerance(
      const int level_number) const = 0;

   /*!
    * @brief Return combine efficiency for clustering tags on level.
    */
   virtual double
   getCombineEfficiency(
      const int level_number) const = 0;

   /*!
    * @brief Write object state out to the given database.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   virtual void
   putToDatabase(
      tbox::Pointer<tbox::Database> db) = 0;

private:

   // TODO: What is the purpose of this data?
   const tbox::Pointer<hier::PatchHierarchy> d_hierarchy;
};

}
}
#endif
