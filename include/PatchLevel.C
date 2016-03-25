//
// File:   $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/patches/PatchLevel.C $
// Package:   SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:   $LastChangedRevision: 1889 $
// Modified:   $LastChangedDate: 2008-01-22 16:46:52 -0800 (Tue, 22 Jan 2008) $
// Description:   A collection of patches at one level of the AMR hierarchy
//

#ifndef included_hier_PatchLevel_C
#define included_hier_PatchLevel_C

#include "PatchLevel.h"

#include <stdio.h>
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define HIER_PATCH_LEVEL_VERSION (3)

#define NAME_BUF_SIZE (32)

#ifdef DEBUG_NO_INLINE
#include "PatchLevel.I"
#endif

namespace SAMRAI {
    namespace hier {

/*
 *************************************************************************
 *                                                                       *
 * Default patch level constructor sets default (non-usable) state.      *
 *                                                                       *
 *************************************************************************
 */

template<int DIM>  PatchLevel<DIM>::PatchLevel()
{
   d_number_patches = 0;

   d_level_number = -1;
   d_next_coarser_level_number = -1;
   d_in_hierarchy = false;
   d_ratio_to_coarser_level = hier::IntVector<DIM>(0);

   d_geometry.setNull();
   d_descriptor.setNull();

   d_box_graph.setNull();
   d_box_top.setNull();
   d_box_tree.setNull();
   d_binary_tree.setNull();

   d_factory = new hier::PatchFactory<DIM>();

   d_patches.resizeArray(0);
   d_patch_touches_regular_boundary.resizeArray(0);
   d_patch_touches_periodic_boundary.resizeArray(0);
   d_shifts.resizeArray(0);
}

/*
 *************************************************************************
 *                                                                       *
 * Create a new patch level using the specified boxes and processor      *
 * mapping.  Only those patches that are local to the processor are      *
 * allocated.  Allocate patches using the specified patch factory or     *
 * the standard patch factory if none is explicitly specified.           *
 *                                                                       *
 *************************************************************************
 */

template<int DIM>  PatchLevel<DIM>::PatchLevel(
   const BoxArray<DIM>& boxes,
   const ProcessorMapping& mapping,
   const IntVector<DIM>& ratio_to_level_zero,
   const tbox::Pointer< GridGeometry<DIM> > grid_geometry,
   const tbox::Pointer< PatchDescriptor<DIM> > descriptor,
   tbox::Pointer< PatchFactory<DIM> > factory,
   bool defer_boundary_box_creation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(boxes.getNumberOfBoxes() == mapping.getSizeOfMappingArray());
   TBOX_ASSERT(!grid_geometry.isNull());
   TBOX_ASSERT(!descriptor.isNull());
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != hier::IntVector<DIM>(0));

   if (DIM > 1) {
      for (int i = 0; i < DIM; i++) {
	 TBOX_ASSERT( (ratio_to_level_zero(i)*ratio_to_level_zero((i+1)%DIM) > 0)
		 || (ratio_to_level_zero(i) == 1)
		 || (ratio_to_level_zero((i+1)%DIM) == 1) );
      }
   }
#endif

   d_number_patches = boxes.getNumberOfBoxes();
   d_boxes = boxes;
   d_mapping.setProcessorMapping(mapping.getProcessorMapping());
   d_descriptor = descriptor;

   d_geometry = grid_geometry;

   d_ratio_to_level_zero = ratio_to_level_zero;

   d_level_number = -1;
   d_next_coarser_level_number = -1;
   d_in_hierarchy = false;
   d_ratio_to_coarser_level = hier::IntVector<DIM>(0);

   d_box_graph.setNull();
   d_box_top.setNull();
   d_box_tree.setNull();
   d_binary_tree.setNull();

   if (!factory.isNull()) {
      d_factory = factory;
   } else {
      d_factory = new hier::PatchFactory<DIM>();
   }

   
   d_patches.resizeArray(d_number_patches);
   d_patch_touches_regular_boundary.resizeArray(d_number_patches);
   d_patch_touches_periodic_boundary.resizeArray(d_number_patches);
   d_shifts.resizeArray(d_number_patches);

   for (int ip = 0; ip < d_number_patches; ip++) {
      if (d_mapping.isMappingLocal(ip)) {
         d_patches[ip] = d_factory->allocate(d_boxes.getBox(ip), d_descriptor);
         d_patches[ip]->setPatchNumber(ip);
         d_patches[ip]->setPatchLevelNumber(d_level_number);
         d_patches[ip]->setPatchInHierarchy(d_in_hierarchy);
      }
   }

   d_boundary_boxes_created = false;

   grid_geometry->computePhysicalDomain(d_physical_domain,
                                        d_ratio_to_level_zero);

   tbox::Array< tbox::Array< tbox::Array<bool> > > touches_regular_bdry;
   tbox::Array< tbox::Array< tbox::Array<bool> > > touches_periodic_bdry;

   grid_geometry->findPatchesTouchingBoundaries(
      touches_regular_bdry,
      touches_periodic_bdry,
      *this,
      grid_geometry->getPeriodicShift(d_ratio_to_level_zero),
      d_physical_domain);

   grid_geometry->setGeometryOnPatches(*this,
                                       d_ratio_to_level_zero,
                                       touches_regular_bdry,
                                       touches_periodic_bdry,
                                       defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   setPatchTouchesBoundaryArrays();

   grid_geometry->computeShiftsForLevel(d_shifts, *this, d_physical_domain);

}

/*
 *************************************************************************
 *                                                                       *
 * Create a new patch level from information in the given database.      *
 *                                                                       *
 *************************************************************************
 */
template<int DIM>  PatchLevel<DIM>::PatchLevel(
   tbox::Pointer<tbox::Database> level_database,
   tbox::Pointer< GridGeometry<DIM> > grid_geometry,
   tbox::Pointer< PatchDescriptor<DIM> > descriptor,
   tbox::Pointer< PatchFactory<DIM> > factory,
   ComponentSelector component_selector,
   bool defer_boundary_box_creation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level_database.isNull());
   TBOX_ASSERT(!grid_geometry.isNull());
   TBOX_ASSERT(!descriptor.isNull());
#endif

   d_geometry = grid_geometry;
   d_descriptor = descriptor;

   if (!factory.isNull()) {
      d_factory = factory;
   } else {
      d_factory = new PatchFactory<DIM>();
   }

   d_box_graph.setNull();
   d_box_top.setNull();
   d_box_tree.setNull();
   d_binary_tree.setNull();

   getFromDatabase(level_database, component_selector);

   d_patch_touches_regular_boundary.resizeArray(d_number_patches);
   d_patch_touches_periodic_boundary.resizeArray(d_number_patches);
   d_shifts.resizeArray(d_number_patches);

   d_boundary_boxes_created = false;

   tbox::Array< tbox::Array< tbox::Array<bool> > > touches_regular_bdry;
   tbox::Array< tbox::Array< tbox::Array<bool> > > touches_periodic_bdry;

   grid_geometry->findPatchesTouchingBoundaries(
      touches_regular_bdry,
      touches_periodic_bdry,
      *this,
      grid_geometry->getPeriodicShift(d_ratio_to_level_zero),
      d_physical_domain);

   grid_geometry->setGeometryOnPatches(*this,
                                       d_ratio_to_level_zero,
                                       touches_regular_bdry,
                                       touches_periodic_bdry,
                                       defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   setPatchTouchesBoundaryArrays();

   grid_geometry->computeShiftsForLevel(d_shifts, *this, d_physical_domain);
}

template<int DIM>  PatchLevel<DIM>::~PatchLevel()
{
}


/*
 * ************************************************************************
 *                                                                       *
 * Allocate or deallocate data for single components or collections of   *
 * component on all patches on a patch level from specified memory pool. *
 *                                                                       *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::allocatePatchData(
   const int id,
   const double timestamp,
   tbox::Pointer<tbox::Arena> pool)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->allocatePatchData(id, timestamp, pool);
   }
}

template<int DIM> void PatchLevel<DIM>::allocatePatchData(
   const ComponentSelector& components,
   const double timestamp,
   tbox::Pointer<tbox::Arena> pool)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->allocatePatchData(components, timestamp, pool);
   }
}

template<int DIM> bool PatchLevel<DIM>::checkAllocated(const int id) const
{
   bool allocated = true;
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      allocated &= d_patches[ip()]->checkAllocated(id);
   }
   return(allocated);
}

template<int DIM> void PatchLevel<DIM>::deallocatePatchData(const int id)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->deallocatePatchData(id);
   }
}

template<int DIM> void PatchLevel<DIM>::deallocatePatchData(
   const ComponentSelector& components)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->deallocatePatchData(components);
   }
}

/*
 * ************************************************************************
 *                                                                        * 
 * Set the simulation time for all patches in the patch level.            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::setTime(const double timestamp, const int id)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->setTime(timestamp, id);
   }
}

template<int DIM> void PatchLevel<DIM>::setTime(const double timestamp,
                               const ComponentSelector& components)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->setTime(timestamp, components);
   }
}

template<int DIM> void PatchLevel<DIM>::setTime(const double timestamp)
{
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {
      d_patches[ip()]->setTime(timestamp);
   }
}

/*
 * ************************************************************************
 *                                                                       *
 * Set level numbers relating this level to la level in a hierarchy    *
 *                                                                       *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::setLevelNumber(const int l)
{
   d_level_number = l;

   for (int i = 0; i < d_number_patches; i++) {
      if (d_mapping.isMappingLocal(i)) {
         d_patches[i]->setPatchLevelNumber(d_level_number); 
      }
   }
}

template<int DIM> void PatchLevel<DIM>::setNextCoarserHierarchyLevelNumber(const int l)
{
   d_next_coarser_level_number = l;
}

/*
 * ************************************************************************
 *                                                                       *
 * Set whether this level resides in a hierarchy.                        *
 *                                                                       *
 * ************************************************************************
 */
 
template<int DIM> void PatchLevel<DIM>::setLevelInHierarchy(bool in_hierarchy)
{
   d_in_hierarchy = in_hierarchy;
 
   for (int i = 0; i < d_number_patches; i++) {
      if (d_mapping.isMappingLocal(i)) {
         d_patches[i]->setPatchInHierarchy(d_in_hierarchy);
      }
   }
}



/*
 * ************************************************************************
 *                                                                       *
 * Set data members of this patch level by refining information on       *
 * the argument level by the given ratio.                                *
 *                                                                       *
 * ************************************************************************
 */

template<int DIM> void hier::PatchLevel<DIM>::setRefinedPatchLevel(
   const tbox::Pointer<hier::PatchLevel<DIM> > coarse_level,
   const hier::IntVector<DIM> & refine_ratio,
   const tbox::Pointer<hier::GridGeometry<DIM> > fine_grid_geometry,
   bool defer_boundary_box_creation)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coarse_level.isNull());
   TBOX_ASSERT(refine_ratio > hier::IntVector<DIM>(0));
#endif

   /*
    * The basic state of the new patch level is initialized from the state of
    * the given existing patch level.
    */

   d_number_patches = coarse_level->d_number_patches;
   d_descriptor = coarse_level->d_descriptor;
   d_factory = coarse_level->d_factory;
   d_mapping.
      setProcessorMapping( (coarse_level->d_mapping).getProcessorMapping() );

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

      const hier::IntVector<DIM>& coarse_ratio = coarse_level->getRatio();
      for (int i = 0; i < DIM; i++) {
         int coarse_rat = coarse_ratio(i);
         int refine_rat = refine_ratio(i);
         if (coarse_rat < 0) {
            if ( tbox::MathUtilities<int>::Abs(coarse_rat) >= refine_rat ) {
               d_ratio_to_level_zero(i) =
                  - ( tbox::MathUtilities<int>::Abs(coarse_rat / refine_rat) );
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

   d_physical_domain = coarse_level->d_physical_domain;
   d_physical_domain.refine(refine_ratio);

   /*
    * Allocate arrays of patches and patch information.  Then, allocate and
    * initialize patch objects.  Finally, set patch geometry and remaining
    * domain information.
    */ 

   d_patches.resizeArray(d_number_patches);
   d_patch_touches_regular_boundary.resizeArray(d_number_patches);
   d_patch_touches_periodic_boundary.resizeArray(d_number_patches);
   d_shifts.resizeArray(d_number_patches);

   for (int p = 0; p < d_number_patches; p++) {
      if (d_mapping.isMappingLocal(p)) {
         d_patches[p] = d_factory->allocate(d_boxes.getBox(p), d_descriptor);
         d_patches[p]->setPatchNumber(p);
         d_patches[p]->setPatchLevelNumber(d_level_number);
         d_patches[p]->setPatchInHierarchy(d_in_hierarchy);
      }
   }

   tbox::Array< tbox::Array< tbox::Array<bool> > >
      touches_regular_bdry(d_number_patches);
   tbox::Array< tbox::Array< tbox::Array<bool> > >
      touches_periodic_bdry(d_number_patches);

   for (typename PatchLevel<DIM>::Iterator ip(coarse_level); ip; ip++) {
      tbox::Pointer< PatchGeometry<DIM> > coarse_pgeom =
         coarse_level->getPatch(ip())->getPatchGeometry();

      touches_regular_bdry[ip()].resizeArray(DIM);
      touches_periodic_bdry[ip()].resizeArray(DIM);

      for (int axis = 0; axis < DIM; axis++) {
         touches_regular_bdry[ip()][axis].resizeArray(2);
         touches_periodic_bdry[ip()][axis].resizeArray(2); 

         for (int side = 0; side < 2; side++) {
            touches_regular_bdry[ip()][axis][side] =
               coarse_pgeom->getTouchesRegularBoundary(axis, side);
            touches_periodic_bdry[ip()][axis][side] =
               coarse_pgeom->getTouchesPeriodicBoundary(axis, side);
         }
      }
   }

   d_geometry->setGeometryOnPatches(*this,
                                    d_ratio_to_level_zero,
                                    touches_regular_bdry,
                                    touches_periodic_bdry,
                                    defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   setPatchTouchesBoundaryArrays();

   d_geometry->computeShiftsForLevel(d_shifts, *this, d_physical_domain);

}

/*
 * ************************************************************************
 *                                                                       *
 * Set data members of this patch level by coarsening information on     *
 * the argument level by the given ratio.                                *
 *                                                                       *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::setCoarsenedPatchLevel(
   const tbox::Pointer<hier::PatchLevel<DIM> > fine_level,
   const hier::IntVector<DIM>& coarsen_ratio,
   const tbox::Pointer<hier::GridGeometry<DIM> > coarse_grid_geom,
   bool defer_boundary_box_creation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fine_level.isNull());
   TBOX_ASSERT(coarsen_ratio > hier::IntVector<DIM>(0));
#endif

   /*
    * The basic state of the new patch level is initialized from the state of
    * the given existing patch level.
    */

   d_number_patches = fine_level->d_number_patches;
   d_descriptor = fine_level->d_descriptor;
   d_factory = fine_level->d_factory;
   d_mapping.
      setProcessorMapping( (fine_level->d_mapping).getProcessorMapping() );

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
   
      const hier::IntVector<DIM>& fine_ratio =
         fine_level->d_ratio_to_level_zero;

      for (int i = 0; i < DIM; i++) {
         int fine_rat = fine_ratio(i);
         int coarsen_rat = coarsen_ratio(i);
         if (fine_rat > 0) {
            if (fine_rat >= coarsen_rat) {
               d_ratio_to_level_zero(i) = fine_rat / coarsen_rat;
            } else {
               d_ratio_to_level_zero(i) =
                  - ( tbox::MathUtilities<int>::Abs(coarsen_rat / fine_rat) );
            } 
         } else {
            d_ratio_to_level_zero(i) =
               - ( tbox::MathUtilities<int>::Abs(fine_rat * coarsen_rat) ); 
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

   d_physical_domain = fine_level->d_physical_domain;
   d_physical_domain.coarsen(coarsen_ratio);

   /*
    * Allocate arrays of patches and patch information.  Then, allocate and
    * initialize patch objects.  Finally, set patch geometry and remaining
    * domain information.
    */ 

   d_patches.resizeArray(d_number_patches);
   d_patch_touches_regular_boundary.resizeArray(d_number_patches);
   d_patch_touches_periodic_boundary.resizeArray(d_number_patches);
   d_shifts.resizeArray(d_number_patches);

   for (int p = 0; p < d_number_patches; p++) {
      if (d_mapping.isMappingLocal(p)) {
         d_patches[p] = d_factory->allocate(d_boxes.getBox(p), d_descriptor);
         d_patches[p]->setPatchNumber(p);
         d_patches[p]->setPatchLevelNumber(d_level_number);
         d_patches[p]->setPatchInHierarchy(d_in_hierarchy);
      }
   }

   d_boundary_boxes_created = false;

   tbox::Array< tbox::Array< tbox::Array<bool> > >
      touches_regular_bdry(d_number_patches);
   tbox::Array< tbox::Array< tbox::Array<bool> > >
      touches_periodic_bdry(d_number_patches);
                                                                                
   for (typename PatchLevel<DIM>::Iterator ip(fine_level); ip; ip++) {
      tbox::Pointer< PatchGeometry<DIM> > fine_pgeom =
         fine_level->getPatch(ip())->getPatchGeometry();
                                                                                
      touches_regular_bdry[ip()].resizeArray(DIM);
      touches_periodic_bdry[ip()].resizeArray(DIM);
                                                                                
      for (int axis = 0; axis < DIM; axis++) {
         touches_regular_bdry[ip()][axis].resizeArray(2);
         touches_periodic_bdry[ip()][axis].resizeArray(2);
                                                                                
         for (int side = 0; side < 2; side++) {
            touches_regular_bdry[ip()][axis][side] =
               fine_pgeom->getTouchesRegularBoundary(axis, side);
            touches_periodic_bdry[ip()][axis][side] =
               fine_pgeom->getTouchesPeriodicBoundary(axis, side);
         }
      }
   }

   d_geometry->setGeometryOnPatches(*this,
                                    d_ratio_to_level_zero,
                                    touches_regular_bdry,
                                    touches_periodic_bdry,
                                    defer_boundary_box_creation);

   setPatchTouchesBoundaryArrays();

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   d_geometry->computeShiftsForLevel(d_shifts, *this, d_physical_domain);

}

/*
 * ************************************************************************
 *                                                                        *
 * Call the geometry routine to create and set boundary boxes, if they    *
 * have not already been created.                                         *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::setBoundaryBoxes()
{
   if (!d_boundary_boxes_created) {
      d_geometry->setBoundaryBoxes(*this);
      d_boundary_boxes_created = true;
   }
}

/*
 * ************************************************************************
 *                                                                        *
 *  Check that class version and restart file number are the same.  If    *
 *  so, read in data from database and build patch level from data.       *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::getFromDatabase(
   tbox::Pointer<tbox::Database> database,
   ComponentSelector component_selector)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   int ver = database->getInteger("HIER_PATCH_LEVEL_VERSION");
   if (ver != HIER_PATCH_LEVEL_VERSION) {
      TBOX_ERROR("PatchLevel<DIM>::getFromDatabase() error...\n"
              << "   Restart file version different than class version.");
   }

   d_boxes = database->getDatabaseBoxArray("d_boxes");

   d_number_patches = database->getInteger("d_number_patches");

   if (!(d_number_patches > 0)) {
      TBOX_ERROR("PatchLevel<DIM>::getFromDatabase() error...\n"
                 << "   d_number_patches in database is not positive" << std::endl);
   }

   if (d_number_patches != d_boxes.getNumberOfBoxes()) {
      TBOX_ERROR("PatchLevel<DIM>::getFromDatabase() error...\n"
                 << "   d_number_patches not same as the number of boxes on "
                 << "this level" << std::endl);
   }

   d_mapping.setProcessorMapping(database->getIntegerArray("d_mapping"));

   if (d_number_patches != d_mapping.getSizeOfMappingArray()) {
      TBOX_ERROR("PatchLevel<DIM>::getFromDatabase() error...\n"
                 << "   d_number_patches not same as size of processor "
                 << "mapping array" << std::endl);
   }

   int* temp_ratio = d_ratio_to_level_zero;
   database->getIntegerArray("d_ratio_to_level_zero", temp_ratio, DIM);

   d_physical_domain = database->getDatabaseBoxArray("d_physical_domain"); 

   d_level_number = database->getInteger("d_level_number");
   d_next_coarser_level_number = 
      database->getInteger("d_next_coarser_level_number");
   d_in_hierarchy = database->getBool("d_in_hierarchy");

   temp_ratio = d_ratio_to_coarser_level;
   database->getIntegerArray("d_ratio_to_coarser_level", temp_ratio, DIM);

   /*
    * Create local processors from database
    */

   d_patches.resizeArray(d_number_patches);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( NAME_BUF_SIZE > (5 + 1 + 4 + 1 + 5 + 1 + 4 + 1) );
#endif
   char patch_name[NAME_BUF_SIZE];
   tbox::Pointer<tbox::Database> patch_database;

   for (int i = 0; i < d_number_patches; i++) {
      if (d_mapping.isMappingLocal(i)) {

         sprintf(patch_name, "level_%04d-patch_%04d", d_level_number, i);

         if (!(database->isDatabase(patch_name))) {
            TBOX_ERROR("PatchLevel<DIM>::getFromDatabase() error...\n"
                       << "   patch name " << patch_name 
                       << " not found in database" << std::endl);
         }
         patch_database = database->getDatabase(patch_name);

         d_patches[i] = d_factory->allocate(d_boxes(i), d_descriptor);
         d_patches[i]->setPatchNumber(i);
         d_patches[i]->setPatchLevelNumber(d_level_number);
         d_patches[i]->setPatchInHierarchy(d_in_hierarchy);
         d_patches[i]->getFromDatabase(patch_database, component_selector);
      }
   }

}

/*
 * ************************************************************************
 *                           *
 *  Write out class version number and patch_level data members to the    *
 *  database, then has each patch on the local processor write itself    *
 *  to the database.   The following are written out to the database:    *
 *  d_physical_domain, d_ratio_to_level_zero, d_boxes, d_mapping,        *
 *  d_number_patches, d_level_number, d_next_coarser_level_number,       *
 *  d_in_hierarchy, d_patches[].                                         *
 *  The database key for all data members except for d_patches is        *
 *  the same as the variable name.  For the patches, the database keys   *
 *  are "level_Xpatch_Y" where X is the level number and Y is the index  *
 *  position of the patch in the patch in d_patches.                     *
 *                           *
 * ************************************************************************
 */
template<int DIM> void PatchLevel<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> database,
   const ComponentSelector& patchdata_write_table)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("HIER_PATCH_LEVEL_VERSION", HIER_PATCH_LEVEL_VERSION);

   database->putBool("d_is_patch_level", true);

   tbox::Array<tbox::DatabaseBox> temp_boxes = d_boxes;
   database->putDatabaseBoxArray("d_boxes", temp_boxes);

   database->putInteger("d_number_patches",d_number_patches);

   database->putIntegerArray("d_mapping", d_mapping.getProcessorMapping());

   int* temp_ratio_to_level_zero = d_ratio_to_level_zero;
   database->putIntegerArray("d_ratio_to_level_zero",
                              temp_ratio_to_level_zero, DIM);

   tbox::Array<tbox::DatabaseBox> temp_domain = d_physical_domain;
   database->putDatabaseBoxArray("d_physical_domain", temp_domain);

   database->putInteger("d_level_number", d_level_number);
   database->putInteger("d_next_coarser_level_number",
                        d_next_coarser_level_number);
   database->putBool("d_in_hierarchy", d_in_hierarchy);

   int* temp_ratio_to_coarser_level = d_ratio_to_coarser_level;
   database->putIntegerArray("d_ratio_to_coarser_level", 
                             temp_ratio_to_coarser_level, DIM);

/*
 * Put local patches in database.
 */

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( NAME_BUF_SIZE > (5 + 1 + 4 + 1 + 5 + 1 + 4 + 1) );
#endif
   char patch_name[NAME_BUF_SIZE];
   tbox::Pointer<tbox::Database> patch_database;
   for (typename PatchLevel<DIM>::Iterator ip(this); ip; ip++) {

      sprintf(patch_name, "level_%04d-patch_%04d", d_level_number, ip());

      patch_database = database->putDatabase(patch_name);

      d_patches[ip()]->putToDatabase(patch_database, patchdata_write_table);
   }

}


template<int DIM> int PatchLevel<DIM>::recursivePrint( std::ostream &os ,
                                                       const std::string &border ,
                                                       unsigned short depth )
{

// Disable Intel warnings on conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#pragma warning (disable:857)
#endif

  int npatch = getNumberOfPatches();
  os << border << "Number of patches = " << npatch << "\n";
  if ( depth > 0 ) {
    for ( Iterator pi(this); pi; pi++ ) {
      int pn = *pi;
      os << border << "Patch " << pn << '/' << npatch << "\n";
      getPatch(pn)->recursivePrint( os, border+"\t", depth - 1);

    }
  }
  return 0;
}

/*
 * ************************************************************************
 * 
 * Returns a BoxGraph object that has been constructed using the
 * level's boxes.  The BoxGraph is constructed the first time
 * the method is called.
 * 
 * ************************************************************************
 */

template<int DIM> tbox::Pointer< BoxGraph<DIM> > PatchLevel<DIM>::getBoxGraph()
{
   if (d_box_graph.isNull()) {

      IntVector<DIM> dst_growth = d_descriptor->getMaxGhostWidth();
      int max_gcw = tbox::MathUtilities<int>::Max(dst_growth.max(),1);
      
      d_box_graph = new BoxGraph<DIM>(d_boxes,
                                       d_shifts,
                                       d_mapping,
                                       d_boxes,
                                       max_gcw);
   }
   return d_box_graph;
}

/*
 * ************************************************************************
 * 
 * Returns a BinaryTree object that has been constructed using the
 * level's boxes and processor mapping.  The BinaryTree is constructed
 * the first time the method is called.
 * 
 * ************************************************************************
 */

template<int DIM> tbox::Pointer< BinaryTree<DIM> > PatchLevel<DIM>::getBinaryTree()
{ 
   if (d_binary_tree.isNull()) {
      d_binary_tree = new hier::BinaryTree<DIM>(d_mapping, d_boxes);
   }
   
   return d_binary_tree;
}

/*
 * ************************************************************************
 * 
 * Returns a BoxTop object that has been constructed using the
 * level's boxes.  The BoxTop is constructed the first time
 * the method is called.
 * 
 * ************************************************************************
 */

template<int DIM> tbox::Pointer< BoxTop<DIM> > PatchLevel<DIM>::getBoxTop()
{
   if (d_box_top.isNull()) {
      d_box_top = new BoxTop<DIM>(d_boxes, d_shifts);
   }
   return d_box_top;
}

/*
 * ************************************************************************
 * 
 * Returns a BoxTree object that has been constructed using the
 * level's boxes.  The BoxTree is constructed the first time
 * the method is called.
 * 
 * ************************************************************************
 */

template<int DIM> tbox::Pointer< BoxTree<DIM> > PatchLevel<DIM>::getBoxTree()
{
   if (d_box_tree.isNull()) {
      d_box_tree = new BoxTree<DIM>(d_boxes, d_shifts, d_mapping);
   }
   return d_box_tree;
}

/*
 * ************************************************************************
 *
 * Private member function that sets arrays indicating whether patches
 * touch domain boundary.  NOTE: This function should only be called 
 * after the grid geometry has set the patch geometry on each patch.
 *
 * ************************************************************************
 */

template<int DIM> void PatchLevel<DIM>::setPatchTouchesBoundaryArrays()
{
   int npatches = getNumberOfPatches();

   d_patch_touches_regular_boundary.resizeArray(npatches);
   d_patch_touches_periodic_boundary.resizeArray(npatches);

   tbox::Array<int> tmp_reg_array(npatches);
   tbox::Array<int> tmp_per_array(npatches);

   for (int it = 0; it < npatches; it++) {
      tmp_reg_array[it] = 0;
      tmp_per_array[it] = 0;
   }

   for (typename PatchLevel<DIM>::Iterator p(this); p; p++) {
      if (getPatch(p())->getPatchGeometry()->getTouchesRegularBoundary()) {
         tmp_reg_array[p()] = 1;
      }
      if (getPatch(p())->getPatchGeometry()->getTouchesPeriodicBoundary()) {
         tmp_per_array[p()] = 1;
      }
   }

   tbox::SAMRAI_MPI::sumReduction(tmp_reg_array.getPointer(), npatches);
   tbox::SAMRAI_MPI::sumReduction(tmp_per_array.getPointer(), npatches);

   for (int ip = 0; ip < npatches; ip++) {
      d_patch_touches_regular_boundary[ip] = 
         ((tmp_reg_array[ip] == 1) ? true : false);
      d_patch_touches_periodic_boundary[ip] = 
         ((tmp_per_array[ip] == 1) ? true : false);
   }
   
}


/*
 * ************************************************************************
 * 
 * Create the patch level iterator and advance through the patches on   *
 * a level of the hierarchy.  The iterator will always point to a valid   *
 * patch or it will be at the end of the patch list.  The iterator will   *
 * only enumerate those patches that are local to the processor.      *
 *
 * ************************************************************************
 */

template<int DIM> PatchLevelIterator<DIM>::PatchLevelIterator(const PatchLevel<DIM>& pl)
:  d_patch(0)
{
   d_number_patches = pl.getProcessorMapping().getNumberOfLocalIndices();
   d_local_box_indices = &(pl.getProcessorMapping().getLocalIndices());
}

template<int DIM> PatchLevelIterator<DIM>::PatchLevelIterator(const PatchLevel<DIM>* pl)
:  d_patch(0)
{

   d_number_patches = pl->getProcessorMapping().getNumberOfLocalIndices();
   d_local_box_indices = &(pl->getProcessorMapping().getLocalIndices());
}

template<int DIM> void PatchLevelIterator<DIM>::initialize(const PatchLevel<DIM>& pl)
{
   d_patch = 0;
   d_number_patches = pl.getProcessorMapping().getNumberOfLocalIndices();
   d_local_box_indices = &(pl.getProcessorMapping().getLocalIndices());
}

template<int DIM> void PatchLevelIterator<DIM>::initialize(const PatchLevel<DIM>* pl)
{
   d_patch = 0;
   d_number_patches = pl->getProcessorMapping().getNumberOfLocalIndices();
   d_local_box_indices = &(pl->getProcessorMapping().getLocalIndices());
}

}
}

#endif
