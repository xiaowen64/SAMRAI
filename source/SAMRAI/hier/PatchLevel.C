/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A collection of patches at one level of the AMR hierarchy
 *
 ************************************************************************/

#ifndef included_hier_PatchLevel_C
#define included_hier_PatchLevel_C

#include "SAMRAI/hier/PatchLevel.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#include <cstdio>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/PatchLevel.I"
#endif

namespace SAMRAI {
namespace hier {

const int PatchLevel::HIER_PATCH_LEVEL_VERSION = 3;

static tbox::Pointer<tbox::Timer> t_level_constructor;
static tbox::Pointer<tbox::Timer> t_constructor_setup;
static tbox::Pointer<tbox::Timer> t_constructor_phys_domain;
static tbox::Pointer<tbox::Timer> t_constructor_touch_boundaries;
static tbox::Pointer<tbox::Timer> t_constructor_set_geometry;
static tbox::Pointer<tbox::Timer> t_set_patch_touches;
static tbox::Pointer<tbox::Timer> t_constructor_compute_shifts;

tbox::StartupShutdownManager::Handler
PatchLevel::s_initialize_finalize_handler(
   PatchLevel::initializeCallback,
   0,
   0,
   PatchLevel::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *
 * Default patch level constructor sets default (non-usable) state.
 *
 *************************************************************************
 */

PatchLevel::PatchLevel(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_mapped_box_level(NULL),
   d_has_globalized_data(false),
   d_ratio_to_level_zero(hier::IntVector::getZero(dim)),
   d_physical_domain(0),
   d_ratio_to_coarser_level(hier::IntVector::getZero(dim))
{
   t_level_constructor->start();
   d_local_number_patches = 0;

   d_level_number = -1;
   d_next_coarser_level_number = -1;
   d_in_hierarchy = false;

   d_geometry.setNull();
   d_descriptor.setNull();

   d_factory = new hier::PatchFactory();

   t_level_constructor->stop();
}

/*
 *************************************************************************
 *
 * Create a new patch level using the specified boxes and processor
 * mapping.  Only those patches that are local to the processor are
 * allocated.  Allocate patches using the specified patch factory or
 * the standard patch factory if none is explicitly specified.
 *
 *************************************************************************
 */

PatchLevel::PatchLevel(
   const BoxLevel& mapped_box_level,
   const tbox::Pointer<GridGeometry> grid_geometry,
   const tbox::Pointer<PatchDescriptor> descriptor,
   tbox::Pointer<PatchFactory> factory,
   bool defer_boundary_box_creation):
   d_dim(grid_geometry->getDim()),
   d_mapped_box_level(new BoxLevel(mapped_box_level)),
   d_has_globalized_data(false),
   d_ratio_to_level_zero(d_mapped_box_level->getRefinementRatio()),
   d_physical_domain(grid_geometry->getNumberBlocks()),
   d_ratio_to_coarser_level(grid_geometry->getDim(), 0)

{
   d_number_blocks = grid_geometry->getNumberBlocks();

   TBOX_DIM_ASSERT_CHECK_ARGS2(mapped_box_level, *grid_geometry);

   t_level_constructor->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!grid_geometry.isNull());
   TBOX_ASSERT(!descriptor.isNull());
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(mapped_box_level.getRefinementRatio() !=
      hier::IntVector::getZero(getDim()));

   if (getDim().getValue() > 1) {
      for (int i = 0; i < getDim().getValue(); i++) {
         TBOX_ASSERT((mapped_box_level.getRefinementRatio() (i)
                      * mapped_box_level.getRefinementRatio() ((i
                                                                + 1)
                         % getDim().getValue()) > 0)
            || (mapped_box_level.getRefinementRatio() (i) == 1)
            || (mapped_box_level.getRefinementRatio() ((i + 1) % getDim().getValue()) ==
                1));
      }
   }
#endif

   t_constructor_setup->start();

   d_local_number_patches =
      static_cast<int>(d_mapped_box_level->getLocalNumberOfBoxes());
   d_descriptor = descriptor;

   d_geometry = grid_geometry;

   d_level_number = -1;
   d_next_coarser_level_number = -1;
   d_in_hierarchy = false;

   if (!factory.isNull()) {
      d_factory = factory;
   } else {
      d_factory = new hier::PatchFactory();
   }

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni) {
      const Box& mapped_box = *ni;
      const BoxId& ip = mapped_box.getId();
      tbox::Pointer<Patch>& patch = d_patches[ip];
      patch = d_factory->allocate(mapped_box, d_descriptor);
      patch->setPatchLevelNumber(d_level_number);
      patch->setPatchInHierarchy(d_in_hierarchy);
   }

   d_boundary_boxes_created = false;
   t_constructor_setup->stop();

   t_constructor_phys_domain->start();
   for (int nb = 0; nb < d_number_blocks; nb++) {
      grid_geometry->computePhysicalDomain(d_physical_domain[nb],
         d_ratio_to_level_zero, BlockId(nb));
   }
   t_constructor_phys_domain->stop();

   t_constructor_touch_boundaries->start();
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;
   grid_geometry->findPatchesTouchingBoundaries(
      touches_regular_bdry,
      touches_periodic_bdry,
      *this,
      grid_geometry->getPeriodicShift(d_ratio_to_level_zero),
      d_physical_domain);
   t_constructor_touch_boundaries->stop();

   t_constructor_set_geometry->start();
   grid_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);
   t_constructor_set_geometry->stop();

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   t_level_constructor->stop();
}

/*
 *************************************************************************
 *
 * Create a new patch level from information in the given database.
 *
 *************************************************************************
 */

PatchLevel::PatchLevel(
   tbox::Pointer<tbox::Database> level_database,
   tbox::Pointer<GridGeometry> grid_geometry,
   tbox::Pointer<PatchDescriptor> descriptor,
   tbox::Pointer<PatchFactory> factory,
   const ComponentSelector& component_selector,
   bool defer_boundary_box_creation):
   d_dim(grid_geometry->getDim()),
   d_has_globalized_data(false),
   d_ratio_to_level_zero(hier::IntVector(grid_geometry->getDim(),
                                         tbox::MathUtilities<int>::getMax())),
   d_physical_domain(grid_geometry->getNumberBlocks()),
   d_ratio_to_coarser_level(hier::IntVector(grid_geometry->getDim(),
                                            tbox::MathUtilities<int>::getMax()))
{
   d_number_blocks = grid_geometry->getNumberBlocks();

   TBOX_ASSERT(!level_database.isNull());
   TBOX_ASSERT(!grid_geometry.isNull());
   TBOX_ASSERT(!descriptor.isNull());

   t_level_constructor->start();

   d_geometry = grid_geometry;
   d_descriptor = descriptor;

   if (!factory.isNull()) {
      d_factory = factory;
   } else {
      d_factory = new PatchFactory();
   }

   getFromDatabase(level_database, component_selector);

   d_boundary_boxes_created = false;

   t_constructor_touch_boundaries->start();
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;
   grid_geometry->findPatchesTouchingBoundaries(
      touches_regular_bdry,
      touches_periodic_bdry,
      *this,
      grid_geometry->getPeriodicShift(d_ratio_to_level_zero),
      d_physical_domain);
   t_constructor_touch_boundaries->stop();

   t_constructor_set_geometry->start();
   grid_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);
   t_constructor_set_geometry->stop();

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   t_level_constructor->stop();
}

PatchLevel::~PatchLevel()
{
}

/*
 * ************************************************************************
 *
 * Allocate or deallocate data for single components or collections of
 * component on all patches on a patch level.
 *
 * ************************************************************************
 */

void PatchLevel::allocatePatchData(
   const int id,
   const double timestamp)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->allocatePatchData(id, timestamp);
   }
}

void PatchLevel::allocatePatchData(
   const ComponentSelector& components,
   const double timestamp)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->allocatePatchData(components, timestamp);
   }
}

bool PatchLevel::checkAllocated(
   const int id) const
{
   bool allocated = true;
   for (PatchContainer::const_iterator
        mi = d_patches.begin(); mi != d_patches.end(); ++mi) {
      allocated &= (*mi).second->checkAllocated(id);
   }
   return allocated;
}

void PatchLevel::deallocatePatchData(
   const int id)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->deallocatePatchData(id);
   }
}

void PatchLevel::deallocatePatchData(
   const ComponentSelector& components)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->deallocatePatchData(components);
   }
}

/*
 * ************************************************************************
 *
 * Set the simulation time for all patches in the patch level.
 *
 * ************************************************************************
 */

void PatchLevel::setTime(
   const double timestamp,
   const int id)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->setTime(timestamp, id);
   }
}

void PatchLevel::setTime(
   const double timestamp,
   const ComponentSelector& components)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->setTime(timestamp, components);
   }
}

void PatchLevel::setTime(
   const double timestamp)
{
   for (PatchLevel::Iterator ip(this); ip; ip++) {
      ip->setTime(timestamp);
   }
}

/*
 * ************************************************************************
 *
 * Set level numbers relating this level to "level", a level in
 * a hierarchy
 *
 * ************************************************************************
 */

void PatchLevel::setLevelNumber(
   const int level)
{
   d_level_number = level;

   for (PatchLevel::Iterator p(this); p; p++) {
      p->setPatchLevelNumber(d_level_number);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

void PatchLevel::setNextCoarserHierarchyLevelNumber(
   const int level)
{
   d_next_coarser_level_number = level;
}

/*
 * ************************************************************************
 *
 * Set whether this level resides in a hierarchy.
 *
 * ************************************************************************
 */

void PatchLevel::setLevelInHierarchy(
   bool in_hierarchy)
{
   d_in_hierarchy = in_hierarchy;

   for (PatchLevel::Iterator p(this); p; p++) {
      p->setPatchInHierarchy(d_in_hierarchy);
   }
}

/*
 * ************************************************************************
 *
 * Set data members of this patch level by refining information on
 * the argument level by the given ratio.
 *
 * ************************************************************************
 */

void PatchLevel::setRefinedPatchLevel(
   const tbox::Pointer<hier::PatchLevel> coarse_level,
   const hier::IntVector& refine_ratio,
   const tbox::Pointer<hier::GridGeometry> fine_grid_geometry,
   bool defer_boundary_box_creation)
{
   TBOX_ASSERT(!coarse_level.isNull());
   TBOX_ASSERT(refine_ratio > hier::IntVector::getZero(getDim()));
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, *coarse_level, refine_ratio);
   if (!fine_grid_geometry.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *fine_grid_geometry);
   }
#endif

   /*
    * The basic state of the new patch level is initialized from the state of
    * the given existing patch level.
    */

   // d_global_number_patches = coarse_level->d_global_number_patches;
   d_descriptor = coarse_level->d_descriptor;
   d_factory = coarse_level->d_factory;
   // d_mapping.setProcessorMapping( (coarse_level->d_mapping).getProcessorMapping() );

   /*
    * Compute the ratio to coarsest level (reference level in hierarchy --
    * usually level  zero) and set grid geometry for this (fine) level.  If
    * pointer to given fine grid geometry is null, then it is assumed that
    * this level is to use the same grid geometry as the given coarse level
    * and the ratio to level zero is set relative to the give coarse level.
    * Otherwise, use given grid geometry and copy ratio to level zero from
    * given coarse level.
    */

   if (fine_grid_geometry.isNull()) {

      d_geometry = coarse_level->d_geometry;

      const hier::IntVector& coarse_ratio = coarse_level->getRatioToLevelZero();
      for (int i = 0; i < getDim().getValue(); i++) {
         int coarse_rat = coarse_ratio(i);
         int refine_rat = refine_ratio(i);
         if (coarse_rat < 0) {
            if (tbox::MathUtilities<int>::Abs(coarse_rat) >= refine_rat) {
               d_ratio_to_level_zero(i) =
                  -(tbox::MathUtilities<int>::Abs(coarse_rat / refine_rat));
            } else {
               d_ratio_to_level_zero(i) =
                  tbox::MathUtilities<int>::Abs(refine_rat / coarse_rat);
            }
         } else {
            d_ratio_to_level_zero(i) = coarse_rat * refine_rat;
         }

      }

   } else {

      d_geometry = fine_grid_geometry;

      d_ratio_to_level_zero = coarse_level->d_ratio_to_level_zero;
   }

   /*
    * Set global box array and index space for level based on refining
    * coarse level information.
    */

   d_boxes = coarse_level->d_boxes;
   d_boxes.refine(refine_ratio);

   {
#if 0
      BoxContainer mapped_boxes(coarse_level->d_mapped_box_level->getBoxes());
      mapped_boxes.refine(refine_ratio); 
#endif
      d_mapped_box_level = new BoxLevel(
            d_ratio_to_level_zero,
            d_geometry,
            coarse_level->d_mapped_box_level->getMPI());
      coarse_level->d_mapped_box_level->refineBoxes(
         *d_mapped_box_level,
         refine_ratio,
         d_ratio_to_level_zero);
      d_mapped_box_level->finalize();
   }
   d_local_number_patches = coarse_level->getLocalNumberOfPatches();
   d_number_blocks = coarse_level->d_number_blocks;

   d_physical_domain.resizeArray(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      d_physical_domain[nb] = coarse_level->d_physical_domain[nb];
      d_physical_domain[nb].refine(refine_ratio);
   }

   /*
    * Allocate arrays of patches and patch information.  Then, allocate and
    * initialize patch objects.  Finally, set patch geometry and remaining
    * domain information.
    */

   // d_patch_touches_regular_boundary.resizeArray(d_global_number_patches);
   // d_patch_touches_periodic_boundary.resizeArray(d_global_number_patches);
   // d_shifts.resizeArray(d_global_number_patches);

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni) {
      const Box& mapped_box = *ni;
      const BoxId& mapped_box_id = mapped_box.getId();
      d_patches[mapped_box_id] = d_factory->allocate(mapped_box, d_descriptor);
      d_patches[mapped_box_id]->setPatchLevelNumber(d_level_number);
      d_patches[mapped_box_id]->setPatchInHierarchy(d_in_hierarchy);
   }

   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;

   for (PatchLevel::Iterator ip(coarse_level); ip; ip++) {
      tbox::Pointer<PatchGeometry> coarse_pgeom =
         (*ip)->getPatchGeometry();

      /* If map does not contain values create them */
      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_regular_bdry(
         touches_regular_bdry.find(ip->getBox().getId()));
      if (iter_touches_regular_bdry == touches_regular_bdry.end()) {
         iter_touches_regular_bdry = touches_regular_bdry.insert(
               iter_touches_regular_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_periodic_bdry(
         touches_periodic_bdry.find(ip->getBox().getId()));
      if (iter_touches_periodic_bdry == touches_periodic_bdry.end()) {
         iter_touches_periodic_bdry = touches_periodic_bdry.insert(
               iter_touches_periodic_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      PatchGeometry::TwoDimBool&
      touches_regular_bdry_ip((*iter_touches_regular_bdry).second);
      PatchGeometry::TwoDimBool&
      touches_periodic_bdry_ip((*iter_touches_periodic_bdry).second);

      for (int axis = 0; axis < getDim().getValue(); axis++) {
         for (int side = 0; side < 2; side++) {

            touches_regular_bdry_ip(axis, side) =
               coarse_pgeom->getTouchesRegularBoundary(axis, side);

            touches_periodic_bdry_ip(axis, side) =
               coarse_pgeom->getTouchesPeriodicBoundary(axis, side);
         }
      }
   }

   d_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

}

/*
 * ************************************************************************
 *
 * Set data members of this patch level by coarsening information on
 * the argument level by the given ratio.
 *
 * ************************************************************************
 */

void PatchLevel::setCoarsenedPatchLevel(
   const tbox::Pointer<hier::PatchLevel> fine_level,
   const hier::IntVector& coarsen_ratio,
   const tbox::Pointer<hier::GridGeometry> coarse_grid_geom,
   bool defer_boundary_box_creation)
{
   TBOX_ASSERT(!fine_level.isNull());
   TBOX_ASSERT(coarsen_ratio > hier::IntVector::getZero(getDim()));

#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, *fine_level, coarsen_ratio);
   if (!coarse_grid_geom.isNull()) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *coarse_grid_geom);
   }
#endif

   /*
    * The basic state of the new patch level is initialized from the state of
    * the given existing patch level.
    */

   // d_global_number_patches = fine_level->d_global_number_patches;
   d_descriptor = fine_level->d_descriptor;
   d_factory = fine_level->d_factory;
   // d_mapping.setProcessorMapping( (fine_level->d_mapping).getProcessorMapping() );

   /*
    * Compute the ratio to coarsest level (reference level in hierarchy --
    * usually level zero) and set grid geometry for this (coarse) level.  If
    * pointer to a given coarse grid geometry is null, then it is assumed
    * that this level is to use the same grid geometry as the given fine
    * level and the ratio to level zero is set relative to the given fine
    * level.  Otherwise, use given grid geometry and copy ratio to level zero
    * from given fine level.
    */

   if (coarse_grid_geom.isNull()) {

      d_geometry = fine_level->d_geometry;

      const hier::IntVector& fine_ratio =
         fine_level->d_ratio_to_level_zero;

      for (int i = 0; i < getDim().getValue(); i++) {
         int fine_rat = fine_ratio(i);
         int coarsen_rat = coarsen_ratio(i);
         if (fine_rat > 0) {
            if (fine_rat >= coarsen_rat) {
               d_ratio_to_level_zero(i) = fine_rat / coarsen_rat;
            } else {
               d_ratio_to_level_zero(i) =
                  -(tbox::MathUtilities<int>::Abs(coarsen_rat / fine_rat));
            }
         } else {
            d_ratio_to_level_zero(i) =
               -(tbox::MathUtilities<int>::Abs(fine_rat * coarsen_rat));
         }
      }

   } else {

      d_geometry = coarse_grid_geom;

      d_ratio_to_level_zero = fine_level->d_ratio_to_level_zero;
   }

   /*
    * Set global box array and index space for level based on coarsening
    * of fine level information.
    */

   d_boxes = fine_level->d_boxes;
   d_boxes.coarsen(coarsen_ratio);

   /*
    * Set coarse mapped_box_level to be the coarsened version of fine mapped_box_level.
    *
    * NOTE: Some parts of SAMRAI (CoarsenSchedule in particular)
    * assumes that the mapped_box identities are the same between the
    * fine and coarsened levels.
    */
   const BoxLevel& fine_mapped_box_level =
      *fine_level->d_mapped_box_level;
#if 0
   BoxContainer coarsened_mapped_boxes(fine_level->d_mapped_box_level->getBoxes());
   coarsened_mapped_boxes.coarsen(coarsen_ratio);
#endif
   d_mapped_box_level = new BoxLevel(
         d_ratio_to_level_zero,
         d_geometry,
         fine_mapped_box_level.getMPI());
   fine_level->d_mapped_box_level->coarsenBoxes(
      *d_mapped_box_level,
      coarsen_ratio,
      d_ratio_to_level_zero);
   d_mapped_box_level->finalize();
   d_local_number_patches = fine_level->getNumberOfPatches();
   d_number_blocks = fine_level->d_number_blocks;

   d_physical_domain.resizeArray(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      d_physical_domain[nb] = fine_level->d_physical_domain[nb];
      d_physical_domain[nb].coarsen(coarsen_ratio);
   }

   /*
    * Allocate arrays of patches and patch information.  Then, allocate and
    * initialize patch objects.  Finally, set patch geometry and remaining
    * domain information.
    */

   // d_patch_touches_regular_boundary.resizeArray(d_global_number_patches);
   // d_patch_touches_periodic_boundary.resizeArray(d_global_number_patches);
   // d_shifts.resizeArray(d_global_number_patches);

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni) {
      const Box& mapped_box = *ni;
      const BoxId& mapped_box_id = mapped_box.getId();
      d_patches[mapped_box_id] = d_factory->allocate(mapped_box, d_descriptor);
      d_patches[mapped_box_id]->setPatchLevelNumber(d_level_number);
      d_patches[mapped_box_id]->setPatchInHierarchy(d_in_hierarchy);
   }

   d_boundary_boxes_created = false;

   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;

   for (PatchLevel::Iterator ip(fine_level); ip; ip++) {
      tbox::Pointer<PatchGeometry> fine_pgeom =
         (*ip)->getPatchGeometry();

      /* If map does not contain values create them */
      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_regular_bdry(
         touches_regular_bdry.find(ip->getBox().getId()));
      if (iter_touches_regular_bdry == touches_regular_bdry.end()) {
         iter_touches_regular_bdry = touches_regular_bdry.insert(
               iter_touches_regular_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_periodic_bdry(
         touches_periodic_bdry.find(ip->getBox().getId()));
      if (iter_touches_periodic_bdry == touches_periodic_bdry.end()) {
         iter_touches_periodic_bdry = touches_periodic_bdry.insert(
               iter_touches_periodic_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      PatchGeometry::TwoDimBool&
      touches_regular_bdry_ip((*iter_touches_regular_bdry).second);
      PatchGeometry::TwoDimBool&
      touches_periodic_bdry_ip((*iter_touches_periodic_bdry).second);

      for (int axis = 0; axis < getDim().getValue(); axis++) {
         for (int side = 0; side < 2; side++) {
            touches_regular_bdry_ip(axis, side) =
               fine_pgeom->getTouchesRegularBoundary(axis, side);
            touches_periodic_bdry_ip(axis, side) =
               fine_pgeom->getTouchesPeriodicBoundary(axis, side);
         }
      }
   }

   d_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

}

/*
 * ************************************************************************
 *
 * Call the geometry routine to create and set boundary boxes, if they
 * have not already been created.
 *
 * ************************************************************************
 */

void PatchLevel::setBoundaryBoxes()
{
   if (!d_boundary_boxes_created) {
      d_geometry->setBoundaryBoxes(*this);
      d_boundary_boxes_created = true;
   }
}

/*
 * ************************************************************************
 *
 *  Check that class version and restart file number are the same.  If
 *  so, read in data from database and build patch level from data.
 *
 * ************************************************************************
 */

void PatchLevel::getFromDatabase(
   tbox::Pointer<tbox::Database> database,
   const ComponentSelector& component_selector)
{
   TBOX_ASSERT(!database.isNull());

   int ver = database->getInteger("HIER_PATCH_LEVEL_VERSION");
   if (ver != HIER_PATCH_LEVEL_VERSION) {
      TBOX_ERROR("PatchLevel::getFromDatabase() error...\n"
         << "   Restart file version different than class version.");
   }

   if (database->keyExists("d_boxes")) {
      d_boxes = database->getDatabaseBoxArray("d_boxes");
   }

   int* temp_ratio = &d_ratio_to_level_zero[0];
   database->getIntegerArray("d_ratio_to_level_zero", temp_ratio, getDim().getValue());

   d_number_blocks = database->getInteger("d_number_blocks");

   d_physical_domain.resizeArray(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      std::string domain_name = "d_physical_domain_"
         + tbox::Utilities::blockToString(nb);
      d_physical_domain[nb] = database->getDatabaseBoxArray(domain_name);
   }

   d_level_number = database->getInteger("d_level_number");
   d_next_coarser_level_number =
      database->getInteger("d_next_coarser_level_number");
   d_in_hierarchy = database->getBool("d_in_hierarchy");

   temp_ratio = &d_ratio_to_coarser_level[0];
   database->getIntegerArray("d_ratio_to_coarser_level", temp_ratio, getDim().getValue());

   /*
    * Put local patches in database.
    */

   tbox::Pointer<tbox::Database> mbl_database = database->getDatabase(
         "mapped_box_level");
   tbox::Pointer<BoxLevel> mapped_box_level(new BoxLevel(getDim()));
   tbox::ConstPointer<GridGeometry> grid_geometry(getGridGeometry());
   mapped_box_level->getFromDatabase(*mbl_database, grid_geometry);
   d_mapped_box_level = mapped_box_level;

   d_patches.clear();

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   tbox::Pointer<tbox::Database> patch_database;
   for (RealBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni) {
      const Box& mapped_box = *ni;
      const LocalId& local_id = mapped_box.getLocalId();
      const BoxId& mapped_box_id = mapped_box.getId();

      std::string patch_name = "level_" + tbox::Utilities::levelToString(
            d_level_number)
         + "-patch_" + tbox::Utilities::patchToString(local_id.getValue())
         + "-block_"
         + tbox::Utilities::blockToString(
            mapped_box_id.getBlockId().getBlockValue());
      if (!(database->isDatabase(patch_name))) {
         TBOX_ERROR("PatchLevel::getFromDatabase() error...\n"
            << "   patch name " << patch_name
            << " not found in database" << std::endl);
      }
      patch_database = database->getDatabase(patch_name);

      tbox::Pointer<Patch>& patch = d_patches[mapped_box_id];
      patch = d_factory->allocate(mapped_box, d_descriptor);
      patch->setPatchLevelNumber(d_level_number);
      patch->setPatchInHierarchy(d_in_hierarchy);
      patch->getFromDatabase(patch_database, component_selector);
   }

}

/*
 * ************************************************************************
 *
 *  Write out class version number and patch_level data members to the
 *  database, then has each patch on the local processor write itself
 *  to the database.   The following are written out to the database:
 *  d_physical_domain, d_ratio_to_level_zero, d_boxes, d_mapping,
 *  d_global_number_patches, d_level_number, d_next_coarser_level_number,
 *  d_in_hierarchy, d_patches[].
 *  The database key for all data members except for d_patches is
 *  the same as the variable name.  For the patches, the database keys
 *  are "level_Xpatch_Y" where X is the level number and Y is the index
 *  position of the patch in the patch in d_patches.
 *
 * ************************************************************************
 */
void PatchLevel::putToDatabase(
   tbox::Pointer<tbox::Database> database,
   const ComponentSelector& patchdata_write_table)
{
   TBOX_ASSERT(!database.isNull());

   database->putInteger("HIER_PATCH_LEVEL_VERSION", HIER_PATCH_LEVEL_VERSION);

   database->putBool("d_is_patch_level", true);

   tbox::Array<tbox::DatabaseBox> temp_boxes = d_boxes;
   if (temp_boxes.getSize() > 0) {
      database->putDatabaseBoxArray("d_boxes", temp_boxes);
   }

   // database->putInteger("d_global_number_patches",d_global_number_patches);

   // database->putIntegerArray("d_mapping", d_mapping.getProcessorMapping());

   int* temp_ratio_to_level_zero = &d_ratio_to_level_zero[0];
   database->putIntegerArray("d_ratio_to_level_zero",
      temp_ratio_to_level_zero, getDim().getValue());

   database->putInteger("d_number_blocks", d_number_blocks);

   for (int nb = 0; nb < d_number_blocks; nb++) {
      tbox::Array<tbox::DatabaseBox> temp_domain = d_physical_domain[nb];
      std::string domain_name = "d_physical_domain_"
         + tbox::Utilities::blockToString(nb);
      database->putDatabaseBoxArray(domain_name, temp_domain);
   }
   database->putInteger("d_level_number", d_level_number);
   database->putInteger("d_next_coarser_level_number",
      d_next_coarser_level_number);
   database->putBool("d_in_hierarchy", d_in_hierarchy);

   int* temp_ratio_to_coarser_level = &d_ratio_to_coarser_level[0];
   database->putIntegerArray("d_ratio_to_coarser_level",
      temp_ratio_to_coarser_level, getDim().getValue());

   /*
    * Put local patches in database.
    */

   tbox::Pointer<tbox::Database> mbl_database = database->putDatabase(
         "mapped_box_level");
   d_mapped_box_level->putToDatabase(*mbl_database);

   tbox::Pointer<tbox::Database> patch_database;
   for (PatchLevel::Iterator ip(this); ip; ip++) {

      std::string patch_name = "level_" + tbox::Utilities::levelToString(
            d_level_number)
         + "-patch_"
         + tbox::Utilities::patchToString(ip->getLocalId().getValue())
         + "-block_"
         + tbox::Utilities::blockToString(
            ip->getBox().getBlockId().getBlockValue());

      patch_database = database->putDatabase(patch_name);

      ip->putToDatabase(patch_database, patchdata_write_table);
   }

}

int PatchLevel::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int depth)
{
   int npatch = getGlobalNumberOfPatches();

// Disable Intel warnings on conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#pragma warning (disable:857)
#endif

   os << border << "Local/Global number of patches and cells = "
      << getLocalNumberOfPatches() << "/" << getGlobalNumberOfPatches() << "  "
      << getLocalNumberOfCells() << "/" << getGlobalNumberOfCells() << "\n";
   os << getBoxLevel()->format(border, 2) << std::endl;

   if (depth > 0) {
      for (Iterator pi(this); pi; pi++) {
         const tbox::Pointer<Patch> patch = *pi;
         os << border << "Patch " << patch->getLocalId() << '/' << npatch << "\n";
         patch->recursivePrint(os, border + "\t", depth - 1);

      }
   }
   return 0;
}

/*
 *************************************************************************
 * Private utility function to gather and store globalized data, if needed.
 *************************************************************************
 */
void PatchLevel::initializeGlobalizedBoxLevel() const
{
   if (!d_has_globalized_data) {

      const BoxLevel& globalized_mapped_box_level(
         d_mapped_box_level->getGlobalizedVersion());

      const int nboxes = globalized_mapped_box_level.getGlobalNumberOfBoxes();
      d_boxes.clear();
      d_mapping.setMappingSize(nboxes);

      /*
       * Backward compatibility with things requiring global sequential
       * indices (such as the VisIt writer) is provided by the implicit
       * ordering of the mapped_boxes in the nested loops below.
       *
       * Due to this necessary renumbering, the patch number obtained
       * by the PatchLevel::Iterator does not correspond to the
       * global sequential index.
       */
      int count = 0;
      const BoxContainer& mapped_boxes =
         globalized_mapped_box_level.getGlobalBoxes();
      for (hier::RealBoxConstIterator ni(mapped_boxes);
           ni.isValid();
           ++ni) {
         d_mapping.setProcessorAssignment(count, ni->getOwnerRank());
         d_boxes.pushBack(*ni);
         ++count;
      }

      d_has_globalized_data = true;
   }
}

/*
 * ************************************************************************
 * ************************************************************************
 */

void PatchLevel::initializeCallback()
{
   t_level_constructor = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::level_constructor");
   t_constructor_setup = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_setup");
   t_constructor_phys_domain = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_phys_domain");
   t_constructor_touch_boundaries = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_touch_boundaries");
   t_constructor_set_geometry = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::set_geometry");
   t_constructor_compute_shifts = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_compute_shifts");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */

void PatchLevel::finalizeCallback()
{
   t_level_constructor.setNull();
   t_constructor_setup.setNull();
   t_constructor_phys_domain.setNull();
   t_constructor_touch_boundaries.setNull();
   t_constructor_set_geometry.setNull();
   t_set_patch_touches.setNull();
   t_constructor_compute_shifts.setNull();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
