/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management in AMR hierarchy
 *
 ************************************************************************/

#ifndef included_hier_GridGeometry_C
#define included_hier_GridGeometry_C

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/GridGeometry.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxContainerSetIterator.h"
#include "SAMRAI/hier/BoxContainerSetConstIterator.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/BoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <map>
#include <stdlib.h>

#define HIER_GRID_GEOMETRY_VERSION (3)

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/GridGeometry.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

tbox::StartupShutdownManager::Handler
GridGeometry::s_initialize_handler(
   GridGeometry::initializeCallback,
   0,
   0,
   GridGeometry::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

tbox::Pointer<tbox::Timer> GridGeometry::t_find_patches_touching_boundaries;
tbox::Pointer<tbox::Timer> GridGeometry::t_touching_boundaries_init;
tbox::Pointer<tbox::Timer> GridGeometry::t_touching_boundaries_loop;
tbox::Pointer<tbox::Timer> GridGeometry::t_set_geometry_on_patches;
tbox::Pointer<tbox::Timer> GridGeometry::t_set_boundary_boxes;
tbox::Pointer<tbox::Timer> GridGeometry::t_set_geometry_data_on_patches;
tbox::Pointer<tbox::Timer> GridGeometry::t_compute_boundary_boxes_on_level;
tbox::Pointer<tbox::Timer> GridGeometry::t_get_boundary_boxes;

/*
 *************************************************************************
 *
 * Constructors for GridGeometry.  Both set up operator
 * handlers.  However, one initializes data members based on arguments.
 * The other initializes the object based on input database information.
 *
 *************************************************************************
 */
GridGeometry::GridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   tbox::Pointer<TransferOperatorRegistry> op_reg,
   tbox::Pointer<tbox::Database> input_db,
   bool register_for_restart):
   d_dim(dim),
   d_object_name(object_name),
   d_physical_domain(1, BoxList(d_dim)),
   d_domain_tree(0),
   d_domain_with_images(1, BoxSet(d_dim)),
   d_periodic_shift(IntVector::getZero(d_dim)),
   d_max_data_ghost_width(IntVector(d_dim, -1)),
   d_has_enhanced_connectivity(false),
   d_transfer_operator_registry(op_reg)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart && d_registered_for_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);
}

/*
 *************************************************************************
 *
 * Constructors for GridGeometry.  Both set up operator
 * handlers.  However, one initializes data members based on arguments.
 * The other initializes the object based on input database information.
 *
 *************************************************************************
 */
GridGeometry::GridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   tbox::Pointer<TransferOperatorRegistry> op_reg):
   d_dim(dim),
   d_object_name(object_name),
   d_physical_domain(1, BoxList(d_dim)),
   d_domain_tree(0),
   d_domain_with_images(1, BoxSet(d_dim)),
   d_periodic_shift(IntVector::getZero(d_dim)),
   d_max_data_ghost_width(IntVector(d_dim, -1)),
   d_number_blocks(1),
   d_number_of_block_singularities(0),
   d_block_neighbors(1),
   d_singularity(1, BoxList(dim)),
   d_singularity_indices(1),
   d_reduced_connect(1),
   d_has_enhanced_connectivity(false),
   d_transfer_operator_registry(op_reg)
{
   TBOX_ASSERT(!object_name.empty());

   d_registered_for_restart = false;

}

GridGeometry::GridGeometry(
   const std::string& object_name,
   const tbox::Array<BoxList>& domain,
   tbox::Pointer<TransferOperatorRegistry> op_reg,
   bool register_for_restart):
   d_dim(domain[0].getDim()),
   d_object_name(object_name),
   d_physical_domain(domain),
   d_domain_tree(0),
   d_domain_with_images(1, BoxSet(d_dim)),
   d_periodic_shift(IntVector::getZero(d_dim)),
   d_max_data_ghost_width(IntVector(d_dim, -1)),
   d_number_of_block_singularities(0),
   d_has_enhanced_connectivity(false),
   d_transfer_operator_registry(op_reg)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(domain.size() > 0);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int b = 0; b < domain.size(); b++) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, domain[b]);
   }
#endif

   d_registered_for_restart = register_for_restart;
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   d_number_blocks = domain.size();
   d_reduced_connect.resizeArray(d_number_blocks, false);
   d_block_neighbors.resizeArray(d_number_blocks);

   setPhysicalDomain(domain);

}

/*
 *************************************************************************
 *
 * Empty destructor.
 *
 *************************************************************************
 */

GridGeometry::~GridGeometry()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(getObjectName());
   }
}

/*
 *************************************************************************
 *
 * Compute boundary boxes for all patches in patch level.  The domain
 * array describes the interior of the level index space.  Note that
 * boundaries is assumed to be an array of DIM * #patches Arrays of
 * BoundaryBoxes.
 *
 *************************************************************************
 */

void GridGeometry::computeBoundaryBoxesOnLevel(
   std::map<BoxId, PatchBoundaries>& boundaries,
   const PatchLevel& level,
   const IntVector& periodic_shift,
   const IntVector& ghost_width,
   const tbox::Array<BoxList>& domain,
   bool do_all_patches) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim,
      level,
      periodic_shift,
      ghost_width);

   t_compute_boundary_boxes_on_level->start();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_width >= IntVector::getZero(ghost_width.getDim()));

   int num_per_dirs = 0;
   for (int i = 0; i < d_dim.getValue(); i++) {
      if (periodic_shift(i)) {
         num_per_dirs++;
      }
   }
   if (num_per_dirs > 0) {
      TBOX_ASSERT(domain.size() == 1);
   }
#endif

   for (PatchLevel::Iterator ip(&level); ip; ip++) {
      tbox::Pointer<Patch> patch = *ip;
      const BoxId& patch_id = patch->getBox().getId();
      const int block_num = patch_id.getBlockId().getBlockValue();

      if (patch->getPatchGeometry()->getTouchesRegularBoundary() ||
          do_all_patches) {

         const Box& box(patch->getBox());

         /*
          * patch_boundaries is an array of DIM BoxLists for each patch.
          * patch_boundaries[DIM-1] will store boundary boxes of the
          * mapped_box type. If DIM > 1, patch_boundaries[DIM-2] will store
          * boundary boxes of the edge type, and if DIM > 2,
          * patch_boundaries[DIM-3] will store boundary boxes of the face
          * type.
          */

         /*
          * Create new map element if one does not exist.
          * Note can't use [] as this requires a default ctor which we do
          * not have for PatchBoundaries.
          */
         std::map<BoxId, PatchBoundaries>::iterator iter(
            boundaries.find(patch_id));
         if (iter == boundaries.end()) {
            std::pair<BoxId, PatchBoundaries> new_boundaries(patch_id,
                                                             PatchBoundaries(d_dim));
            iter = boundaries.insert(iter, new_boundaries);
         }
         getBoundaryBoxes((*iter).second, box, domain[block_num],
            ghost_width, periodic_shift);

#ifdef DEBUG_CHECK_ASSERTIONS
         for (int j = 0; j < d_dim.getValue(); j++) {
            iter = (boundaries.find(patch_id));
            TBOX_ASSERT(iter != boundaries.end());
            for (int k = 0; k < ((*iter).second)[j].getSize(); k++) {
               TBOX_ASSERT(checkBoundaryBox(((*iter).second)[j][k], *patch,
                     domain[block_num], num_per_dirs, ghost_width));
            }
         }
#endif
      }
   }
   t_compute_boundary_boxes_on_level->stop();
}

/*
 *************************************************************************
 *
 * For each patch in the level, use box intersection operation to
 * determine what kind of boundaries, if any the patch touches.  Call
 * Patch functions to set flags that store this information once it
 * is found.
 *
 *************************************************************************
 */

void GridGeometry::findPatchesTouchingBoundaries(
   std::map<BoxId, TwoDimBool>& touches_regular_bdry,
   std::map<BoxId, TwoDimBool>& touches_periodic_bdry,
   const PatchLevel& level,
   const IntVector& periodic_shift,
   const tbox::Array<BoxList>& domain) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(level, periodic_shift, domain[0]);

   tbox::Array<tbox::Pointer<BoxTree> > domain_tree(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      domain_tree[nb] = new BoxTree(d_dim, domain[nb], BlockId(nb));
   }
   findPatchesTouchingBoundaries(touches_regular_bdry,
      touches_periodic_bdry,
      level,
      periodic_shift,
      domain_tree);
}

void GridGeometry::findPatchesTouchingBoundaries(
   std::map<BoxId, TwoDimBool>& touches_regular_bdry,
   std::map<BoxId, TwoDimBool>& touches_periodic_bdry,
   const PatchLevel& level,
   const IntVector& periodic_shift,
   const tbox::Array<tbox::Pointer<BoxTree> >& domain_tree) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(level, periodic_shift);

   t_find_patches_touching_boundaries->start();

   t_touching_boundaries_init->start();
   touches_regular_bdry.clear();
   touches_periodic_bdry.clear();
   t_touching_boundaries_init->stop();

   t_touching_boundaries_loop->start();
   for (PatchLevel::Iterator ip(&level); ip; ip++) {
      tbox::Pointer<Patch> patch = *ip;
      const Box& box(patch->getBox());
      const int block_number = patch->getBox().getBlockId().getBlockValue();

      std::map<BoxId, TwoDimBool>::iterator iter_touches_regular_bdry(
         touches_regular_bdry.find(ip->getBox().getId()));
      if (iter_touches_regular_bdry == touches_regular_bdry.end()) {
         iter_touches_regular_bdry = touches_regular_bdry.insert(
               iter_touches_regular_bdry,
               std::pair<BoxId, TwoDimBool>(ip->getBox().getId(), TwoDimBool(d_dim)));
      }

      std::map<BoxId, TwoDimBool>::iterator iter_touches_periodic_bdry(
         touches_periodic_bdry.find(ip->getBox().getId()));
      if (iter_touches_periodic_bdry == touches_periodic_bdry.end()) {
         iter_touches_periodic_bdry = touches_periodic_bdry.insert(
               iter_touches_periodic_bdry,
               std::pair<BoxId, TwoDimBool>(ip->getBox().getId(), TwoDimBool(d_dim)));
      }

      computeBoxTouchingBoundaries((*iter_touches_regular_bdry).second,
         (*iter_touches_periodic_bdry).second,
         box,
         periodic_shift,
         *(domain_tree[block_number]));
   }
   t_touching_boundaries_loop->stop();
   t_find_patches_touching_boundaries->stop();
}

void GridGeometry::computeBoxTouchingBoundaries(
   TwoDimBool& touches_regular_bdry,
   TwoDimBool& touches_periodic_bdry,
   const Box& box,
   const IntVector& periodic_shift,
   const BoxTree& domain_tree) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, periodic_shift);

   /*
    * Create a list of boxes inside a mapped_box_level of one cell outside the
    * patch.  Remove the intersections with the domain's interior, so that only
    * boxes outside the physical domain (if any) remain in the list.
    */
   BoxList bdry_list(Box::grow(box, IntVector::getOne(d_dim)));
   bdry_list.removeIntersections(domain_tree);
   const bool touches_any_boundary = (bdry_list.size() > 0);

   if (!touches_any_boundary) {
      for (int d = 0; d < d_dim.getValue(); ++d) {
         touches_regular_bdry(d, 0) = touches_periodic_bdry(d, 0) =
               touches_regular_bdry(d, 1) = touches_periodic_bdry(d, 1) = false;
      }
   } else {
      bool bdry_located = false;
      for (int nd = 0; nd < d_dim.getValue(); nd++) {
         BoxList lower_list(bdry_list);
         BoxList upper_list(bdry_list);

         Box test_box(box);

         test_box.growLower(nd, 1);
         lower_list.intersectBoxes(test_box); // performance ok.  lower_list is short.

         test_box = box;
         test_box.growUpper(nd, 1);
         upper_list.intersectBoxes(test_box); // performance ok.  upper_list is short.

         if (lower_list.size()) {
            // Touches regular or periodic bdry on lower side.
            touches_periodic_bdry(nd, 0) = (periodic_shift(nd) != 0);
            touches_regular_bdry(nd, 0) = (periodic_shift(nd) == 0);
            bdry_located = true;
         }

         if (upper_list.size()) {
            // Touches regular or periodic bdry on upper side.
            touches_periodic_bdry(nd, 1) = (periodic_shift(nd) != 0);
            touches_regular_bdry(nd, 1) = (periodic_shift(nd) == 0);
            bdry_located = true;
         }
      }

      /*
       * By this point, bdry_located will have been set to true almost
       * every time whenever touches_any_boundary is true.  The only way
       * it will not be true is if the domain is not a parallelpiped, and
       * the patch touches the boundary only at a location such as the
       * concave corner of an L-shaped domain.
       */
      if (!bdry_located) {
         for (int nd = 0; nd < d_dim.getValue(); nd++) {
            touches_periodic_bdry(nd, 0) = touches_periodic_bdry(nd, 1) = false;

            bool lower_side = false;
            bool upper_side = false;
            for (BoxContainer::Iterator bl(bdry_list); bl; bl++) {
               if (bl().lower() (nd) < box.lower(nd)) {
                  lower_side = true;
               }
               if (bl().upper() (nd) > box.upper(nd)) {
                  upper_side = true;
               }
               if (lower_side && upper_side) {
                  break;
               }
            }
            touches_regular_bdry(nd, 0) = lower_side;
            touches_regular_bdry(nd, 1) = upper_side;
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Set geometry data for each patch on level.
 *
 *************************************************************************
 */

void GridGeometry::setGeometryOnPatches(
   PatchLevel& level,
   const IntVector& ratio_to_level_zero,
   std::map<BoxId, TwoDimBool>& touches_regular_bdry,
   std::map<BoxId, TwoDimBool>& touches_periodic_bdry,
   bool defer_boundary_box_creation)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, level, ratio_to_level_zero);

   t_set_geometry_on_patches->start();
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != IntVector::getZero(d_dim));
   if (d_dim.getValue() > 1) {
      for (int i = 0; i < d_dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % d_dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   t_set_geometry_data_on_patches->start();
   for (PatchLevel::Iterator ip(&level); ip; ip++) {
      tbox::Pointer<Patch> patch = *ip;
      setGeometryDataOnPatch(*patch, ratio_to_level_zero,
         (*touches_regular_bdry.find(ip->getBox().getId())).second,
         (*touches_periodic_bdry.find(ip->getBox().getId())).second);
   }
   t_set_geometry_data_on_patches->stop();

   if (!defer_boundary_box_creation) {
      setBoundaryBoxes(level);
   }
   t_set_geometry_on_patches->stop();
}

/*
 *************************************************************************
 *
 * Set boundary boxes for each patch on level.
 *
 *************************************************************************
 */

void GridGeometry::setBoundaryBoxes(
   PatchLevel& level)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, level);

   t_set_boundary_boxes->start();
   std::map<BoxId, PatchBoundaries> boundaries;

   tbox::Array<BoxList> domain(level.getPhysicalDomainArray());

   IntVector ghost_width(
      level.getPatchDescriptor()->getMaxGhostWidth(d_dim));

   if (d_max_data_ghost_width != -1 &&
       !(ghost_width <= d_max_data_ghost_width)) {

      TBOX_ERROR("Error in GridGeometry object with name = "
         << d_object_name << ": in computeMaxGhostWidth():  "
         << "Cannot add variables and increase maximum ghost "
         << "width after creating the GridGeometry!");
   }

   d_max_data_ghost_width = ghost_width;

   computeBoundaryBoxesOnLevel(
      boundaries,
      level,
      getPeriodicShift(IntVector::getOne(d_dim)),
      d_max_data_ghost_width,
      domain);

   for (std::map<BoxId, PatchBoundaries>::const_iterator mi = boundaries.begin();
        mi != boundaries.end(); ++mi) {
      tbox::Pointer<Patch> patch = level.getPatch((*mi).first);
      patch->getPatchGeometry()->setBoundaryBoxesOnPatch((*mi).second.getArrays());
   }

   t_set_boundary_boxes->stop();
}

/*
 *************************************************************************
 *
 * Create PatchGeometry geometry object, initializing its
 * boundary and assigning it to the given patch.
 *
 *************************************************************************
 */

void GridGeometry::setGeometryDataOnPatch(
   Patch& patch,
   const IntVector& ratio_to_level_zero,
   const PatchGeometry::TwoDimBool& touches_regular_bdry,
   const PatchGeometry::TwoDimBool& touches_periodic_bdry)
const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, patch, ratio_to_level_zero,
      touches_regular_bdry, touches_periodic_bdry);

   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < dim.getValue(); i++) {
      TBOX_ASSERT(ratio_to_level_zero(i) != 0);
   }
   if (dim > tbox::Dimension(1)) {
      for (i = 0; i < dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % dim.getValue()) == 1));
      }
   }
#endif

   tbox::Pointer<PatchGeometry>
   geometry(new PatchGeometry(ratio_to_level_zero,
               touches_regular_bdry,
               touches_periodic_bdry));

   patch.setPatchGeometry(geometry);

}

/*
 *************************************************************************
 *
 * Create and return pointer to coarsened version of this
 * grid geometry object coarsened by the given ratio.
 *
 *************************************************************************
 */

tbox::Pointer<GridGeometry>
GridGeometry::makeCoarsenedGridGeometry(
   const std::string& coarse_geom_name,
   const IntVector& coarsen_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!coarse_geom_name.empty());
   TBOX_ASSERT(coarse_geom_name != getObjectName());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, coarsen_ratio);
   TBOX_ASSERT(coarsen_ratio > IntVector::getZero(dim));

   tbox::Array<BoxList> coarse_domain(d_number_blocks, BoxList(dim));

   for (int b = 0; b < d_number_blocks; b++) {
      BlockId block_id(b);
      coarse_domain[b] = getPhysicalDomain(block_id);
      coarse_domain[b].coarsen(coarsen_ratio);

      /*
       * Need to check that domain can be coarsened by given ratio.
       */
      const BoxList& fine_domain = getPhysicalDomain(block_id);
      const int nboxes = fine_domain.size();
      BoxContainer::ConstIterator coarse_domain_itr(coarse_domain[b]);
      BoxContainer::ConstIterator fine_domain_itr(fine_domain);
      for (int ib = 0; ib < nboxes; ib++, coarse_domain_itr++, fine_domain_itr++) {
         Box testbox = Box::refine(*coarse_domain_itr, coarsen_ratio);
         if (!testbox.isSpatiallyEqual(*fine_domain_itr)) {
#ifdef DEBUG_CHECK_ASSERTIONS
            tbox::plog
            << "GridGeometry::makeCoarsenedGridGeometry : Box # "
            << ib << std::endl;
            tbox::plog << "   fine box = " << *fine_domain_itr << std::endl;
            tbox::plog << "f   coarse box = " << *coarse_domain_itr << std::endl;
            tbox::plog << "   refined coarse box = " << testbox << std::endl;
#endif
            TBOX_ERROR(
               "GridGeometry::makeCoarsenedGridGeometry() error...\n"
               << "    geometry object with name = " << getObjectName()
               << "\n    Cannot be coarsened by ratio " << coarsen_ratio
               << std::endl);
         }
      }
   }

   GridGeometry* coarse_geometry =
      new GridGeometry(coarse_geom_name,
         coarse_domain,
         d_transfer_operator_registry,
         register_for_restart);

   coarse_geometry->initializePeriodicShift(getPeriodicShift(
         IntVector::getOne(dim)));

   return tbox::Pointer<GridGeometry>(coarse_geometry);
}

/*
 *************************************************************************
 *
 * Create and return pointer to refined version of this
 * grid geometry object refined by the given ratio.
 *
 *************************************************************************
 */

tbox::Pointer<GridGeometry>
GridGeometry::makeRefinedGridGeometry(
   const std::string& fine_geom_name,
   const IntVector& refine_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!fine_geom_name.empty());
   TBOX_ASSERT(fine_geom_name != getObjectName());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, refine_ratio);
   TBOX_ASSERT(refine_ratio > IntVector::getZero(dim));

   tbox::Array<BoxList> fine_domain(d_number_blocks, BoxList(dim));

   for (int b = 0; b < d_number_blocks; b++) {
      fine_domain[b] = getPhysicalDomain(BlockId(b));
      fine_domain[b].refine(refine_ratio);
   }

   GridGeometry* fine_geometry =
      new GridGeometry(fine_geom_name,
         fine_domain,
         d_transfer_operator_registry,
         register_for_restart);

   fine_geometry->initializePeriodicShift(getPeriodicShift(
         IntVector::getOne(dim)));

   return tbox::Pointer<GridGeometry>(fine_geometry);
}

/*
 *************************************************************************
 * Checks to see if the version number for the class is the same as
 * as the version number of the restart file.
 * If they are equal, then the data from the database are read to local
 * variables and the setPhysicalDomain() method is called.
 *
 *************************************************************************
 */
void GridGeometry::getFromRestart()
{
   const tbox::Dimension dim(getDim());

   tbox::Pointer<tbox::Database> restart_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;

   if (restart_db->isDatabase(getObjectName())) {
      db = restart_db->getDatabase(getObjectName());
   } else {
      TBOX_ERROR("Restart database corresponding to "
         << getObjectName() << " not found in the restart file.");
   }

   int ver = db->getInteger("HIER_GRID_GEOMETRY_VERSION");
   if (ver != HIER_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR(
         getObjectName() << ":  "
                         << "Restart file version is different than class version.");
   }

   d_number_blocks = db->getInteger("d_number_blocks");

   std::string domain_name;
   tbox::Array<hier::BoxList> domain(d_number_blocks, hier::BoxList(dim));

   for (int b = 0; b < d_number_blocks; b++) {
      domain_name = "d_physical_domain_" + tbox::Utilities::intToString(b);
      if (db->keyExists(domain_name)) {
         domain[b] = db->getDatabaseBoxArray(domain_name);
      } else {
         TBOX_ERROR(
            getObjectName() << ":  "
                            << "No '" << domain_name << "' restart data found for "
                            << "Block " << b << " physical domain. ");
      }
   }
   setPhysicalDomain(domain);

   hier::IntVector periodic_shift(dim);
   int* temp_shift = &periodic_shift[0];
   db->getIntegerArray("d_periodic_shift", temp_shift, dim.getValue());
   initializePeriodicShift(periodic_shift);

   d_has_enhanced_connectivity = db->getBool("d_has_enhanced_connectivity");
}

/*
 *************************************************************************
 *
 * Data is read from input only if the simulation is not from restart.
 * Otherwise, all values specifed in the input database are ignored.
 * In this method data from the database are read to local
 * variables and the setPhysicalDomain() method is called.
 *
 *************************************************************************
 */

void GridGeometry::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{

   TBOX_ASSERT(!db.isNull());

   const tbox::Dimension dim(getDim());

   if (!is_from_restart) {

      d_number_blocks = db->getIntegerWithDefault("num_blocks", 1);

      std::string domain_name;
      tbox::Array<hier::BoxList> domain(d_number_blocks, hier::BoxList(dim));

      for (int b = 0; b < d_number_blocks; b++) {

         domain_name = "domain_boxes_" + tbox::Utilities::intToString(b);

         if (db->keyExists(domain_name)) {
            domain[b] = db->getDatabaseBoxArray(domain_name);
            if (domain[b].size() == 0) {
               TBOX_ERROR(
                  getObjectName() << ":  "
                                  << "No boxes for " << domain_name
                                  << " array found in input.");
            }
         } else {
            TBOX_ERROR(
               getObjectName() << ":  "
                               << "Key data '" << domain_name << "' not found in input.");
         }
      }

      int pbc[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      hier::IntVector per_bc(dim, 0);
      if (db->keyExists("periodic_dimension")) {
         db->getIntegerArray("periodic_dimension", pbc, dim.getValue());
         for (int i = 0; i < dim.getValue(); i++) {
            per_bc(i) = ((pbc[i] == 0) ? 0 : 1);
         }
      }

      setPhysicalDomain(domain);

      initializePeriodicShift(per_bc);
   }

   readBlockDataFromInput(db);
}

/*
 *************************************************************************
 *
 * Writes out version number and data members for the class.
 *
 *************************************************************************
 */

void GridGeometry::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
   TBOX_ASSERT(!db.isNull());

   const tbox::Dimension dim(getDim());

   db->putInteger("HIER_GRID_GEOMETRY_VERSION",
      HIER_GRID_GEOMETRY_VERSION);

   db->putInteger("d_number_blocks", d_number_blocks);

   std::string domain_name;

   for (int b = 0; b < d_number_blocks; b++) {

      domain_name = "d_physical_domain_" + tbox::Utilities::intToString(b);

      tbox::Array<tbox::DatabaseBox> temp_box_array =
         getPhysicalDomain(BlockId(b));
      db->putDatabaseBoxArray(domain_name, temp_box_array);
   }

   hier::IntVector level0_shift(getPeriodicShift(hier::IntVector::getOne(dim)));
   int* temp_shift = &level0_shift[0];
   db->putIntegerArray("d_periodic_shift", temp_shift, dim.getValue());

   db->putBool("d_has_enhanced_connectivity", d_has_enhanced_connectivity);
}

/*
 * ************************************************************************
 *
 * Compute the valid periodic shifts for the given box.
 *
 * ************************************************************************
 */

void GridGeometry::computeShiftsForBox(
   std::vector<IntVector>& shifts,
   const Box& box,
   const BoxTree& domain_search_tree,
   const IntVector& periodic_shift) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, box, periodic_shift);

   shifts.clear();

   int num_periodic_dirs = 0;

   for (int i = 0; i < d_dim.getValue(); i++) {
      if (periodic_shift(i) != 0) {
         num_periodic_dirs++;
      }
   }

   if (num_periodic_dirs > 0) {

      const PeriodicShiftCatalog* periodic_shift_catalog =
         PeriodicShiftCatalog::getCatalog(d_dim);

      shifts.reserve(periodic_shift_catalog->getNumberOfShifts());

      BoundaryLookupTable* blut =
         BoundaryLookupTable::getLookupTable(d_dim);

      const tbox::Array<int>& location_index_max =
         blut->getMaxLocationIndices();

      for (int d = 0; d < num_periodic_dirs; d++) {

         const int codim = d + 1;

         for (int loc = 0; loc < location_index_max[d]; loc++) {

            const tbox::Array<int>& dirs = blut->getDirections(loc, codim);

            bool need_to_test = true;
            for (int k = 0; k < dirs.size(); k++) {
               if (periodic_shift(dirs[k]) == 0) {
                  need_to_test = false;
                  break;
               }
            }

            if (need_to_test) {

               Box border(box);
               IntVector border_shift(d_dim, 0);

               tbox::Array<bool> is_upper(codim);
               for (int j = 0; j < codim; j++) {
                  if (blut->isUpper(loc, codim, j)) {
                     border.lower(dirs[j]) = box.upper(dirs[j]);
                     border.upper(dirs[j]) = box.upper(dirs[j]);
                     border_shift(dirs[j]) = 1;
                     is_upper[j] = true;
                  } else {
                     border.lower(dirs[j]) = box.lower(dirs[j]);
                     border.upper(dirs[j]) = box.lower(dirs[j]);
                     border_shift(dirs[j]) = -1;
                     is_upper[j] = false;
                  }
               }

               border.shift(border_shift);
               BoxList border_list(border);

               border_list.removeIntersections(domain_search_tree);

               if (border_list.size() > 0) {

                  const Box& domain_bound_box =
                     domain_search_tree.getBoundingBox();

                  if (codim == 1) {

                     IntVector new_shift(d_dim, 0);
                     if (is_upper[0]) {
                        new_shift(dirs[0]) =
                           -domain_bound_box.numberCells(dirs[0]);
                     } else {
                        new_shift(dirs[0]) =
                           domain_bound_box.numberCells(dirs[0]);
                     }
                     // shifts.addItem(new_shift);
                     shifts.insert(shifts.end(), new_shift);

                  } else {

                     bool shift_to_add = true;
                     for (int c = 0; c < codim; c++) {

                        if (is_upper[c]) {
                           if (border.upper(dirs[c]) <=
                               domain_bound_box.upper(dirs[c])) {
                              shift_to_add = false;
                              break;
                           }
                        } else {
                           if (border.lower(dirs[c]) >=
                               domain_bound_box.lower(dirs[c])) {
                              shift_to_add = false;
                              break;
                           }
                        }

                     }

                     if (shift_to_add) {
                        IntVector new_shift(d_dim, 0);
                        for (int b = 0; b < codim; b++) {
                           if (is_upper[b]) {
                              new_shift(dirs[b]) =
                                 -domain_bound_box.numberCells(dirs[b]);
                           } else {
                              new_shift(dirs[b]) =
                                 domain_bound_box.numberCells(dirs[b]);
                           }
                        }
                        // shifts.addItem(new_shift);
                        shifts.insert(shifts.end(), new_shift);
                     }
                  }
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Decompose patch boundary region into pieces depending on spatial dim.
 *
 *************************************************************************
 */

void GridGeometry::getBoundaryBoxes(
   PatchBoundaries& patch_boundaries,
   const Box& box,
   const BoxList& domain_boxes,
   const IntVector& ghosts,
   const IntVector& periodic_shift) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS5(*this, box, domain_boxes, ghosts, periodic_shift);

   t_get_boundary_boxes->start();

   const Index ifirst = box.lower();
   const Index ilast = box.upper();

   int num_per_dirs = 0;
   for (int d = 0; d < d_dim.getValue(); d++) {
      num_per_dirs += (periodic_shift(d) ? 1 : 0);
   }

   if (num_per_dirs == d_dim.getValue()) {
      for (int k = 0; k < d_dim.getValue(); k++) {
         patch_boundaries[k].resizeArray(0, BoundaryBox(d_dim));
      }

   } else {
      BoundaryLookupTable* blut =
         BoundaryLookupTable::getLookupTable(d_dim);

      const tbox::Array<int>& location_index_max =
         blut->getMaxLocationIndices();
      tbox::Array<BoxList> codim_boxlist(d_dim.getValue(), BoxList(d_dim));

      for (int d = 0; d < d_dim.getValue() - num_per_dirs; d++) {

         int codim = d + 1;

         patch_boundaries[d].resizeArray(location_index_max[d],
            BoundaryBox(d_dim));
         int bdry_array_size = location_index_max[d];
         int num_bboxes = 0;

         for (int loc = 0; loc < location_index_max[d]; loc++) {
            const tbox::Array<int>& dirs = blut->getDirections(loc, codim);

            tbox::Array<bool> periodic_dir(codim);
            bool all_is_per = true;
            for (int p = 0; p < codim; p++) {
               if (periodic_shift(dirs[p]) == 0) {
                  periodic_dir[p] = false;
                  all_is_per = false;
               } else {
                  periodic_dir[p] = true;
               }
            }

            if (!all_is_per) {
               Box border(box);
               IntVector border_shift(d_dim, 0);

               for (int i = 0; i < codim; i++) {
                  if (blut->isUpper(loc, codim, i)) {
                     border.lower(dirs[i]) = box.upper(dirs[i]);
                     border.upper(dirs[i]) = box.upper(dirs[i]);
                     border_shift(dirs[i]) = 1;
                  } else {
                     border.lower(dirs[i]) = box.lower(dirs[i]);
                     border.upper(dirs[i]) = box.lower(dirs[i]);
                     border_shift(dirs[i]) = -1;
                  }
               }

               // grow in non-dirs directions
               for (int j = 0; j < d_dim.getValue(); j++) {
                  bool dir_used = false;
                  for (int du = 0; du < codim; du++) {
                     if (dirs[du] == j) {
                        dir_used = true;
                        break;
                     }
                  }
                  if (!dir_used) {
                     border.upper(j) = ilast(j) + ghosts(j);
                     border.lower(j) = ifirst(j) - ghosts(j);
                  }
               }

               BoxList domain_list(domain_boxes);
               if (num_per_dirs != 0) {
                  domain_list.grow(periodic_shift);
               }

               /*
                * Intersect border_list with domain, then shift so that
                * true boundary boxes are outside domain.  Then remove
                * intersections with the domain.
                */

               BoxList border_list(border);
               border_list.intersectBoxes(domain_list); // Should use BoxTree here if possible.

               border_list.shift(border_shift);

               border_list.removeIntersections(domain_list);

               if (border_list.size() > 0) {
                  for (int bd = 0; bd < d; bd++) {
                     border_list.removeIntersections(codim_boxlist[bd]);

                     if (border_list.size() == 0) {
                        break;
                     }
                  }
               }

               if (border_list.size() > 0) {
                  border_list.coalesce();
                  for (BoxContainer::Iterator bl(border_list);
                       bl; bl++) {
                     if (num_bboxes == bdry_array_size) {
                        patch_boundaries[d].resizeArray(
                           bdry_array_size + location_index_max[d],
                           BoundaryBox(d_dim));
                        bdry_array_size = patch_boundaries[d].size();
                     }

                     BoundaryBox boundary_box(bl(), codim, loc);

                     patch_boundaries[d][num_bboxes] = boundary_box;

                     num_bboxes++;
                  }

                  codim_boxlist[d].spliceFront(border_list);
               }
            }

            if (loc + 1 == location_index_max[d]) {
               patch_boundaries[d].resizeArray(num_bboxes,
                  BoundaryBox(d_dim));
            }
         }
      }
   }
   t_get_boundary_boxes->stop();
}

/*
 *************************************************************************
 *
 * Compute physical domain for index space related to reference domain
 * by specified ratio.  If any entry of ratio is negative, the reference
 * domain will be coarsened.  Otherwise, it will be refined.
 *
 *************************************************************************
 */

void GridGeometry::computePhysicalDomain(
   BoxList& domain,
   const IntVector& ratio_to_level_zero,
   const BlockId& block_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, domain, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < d_dim.getValue(); i++) {
      TBOX_ASSERT(ratio_to_level_zero(i) != 0);
   }
   if (d_dim.getValue() > 1) {
      for (i = 0; i < d_dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % d_dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   domain = d_physical_domain[block_id.getBlockValue()];

   bool coarsen = false;
   IntVector tmp_rat = ratio_to_level_zero;
   for (int id = 0; id < d_dim.getValue(); id++) {
      if (ratio_to_level_zero(id) < 0) coarsen = true;
      tmp_rat(id) = abs(ratio_to_level_zero(id));
   }

   if (coarsen) {
      domain.coarsen(tmp_rat);
   } else {
      domain.refine(tmp_rat);
   }

}

/*
 *************************************************************************
 *
 * Compute physical domain for index space related to reference domain
 * by specified ratio.  If any entry of ratio is negative, the reference
 * domain will be coarsened.  Otherwise, it will be refined.
 *
 *************************************************************************
 */
#if 0
void GridGeometry::computePhysicalDomain(
   BoxSet& domain_mapped_boxes,
   const IntVector& ratio_to_level_zero,
   const BlockId& block_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < d_dim.getValue(); i++) {
      TBOX_ASSERT(ratio_to_level_zero(i) != 0);
   }
   if (d_dim.getValue() > 1) {
      for (i = 0; i < d_dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % d_dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   domain_mapped_boxes = d_domain_with_images[block_id.getBlockValue()];

   if (ratio_to_level_zero != IntVector::getOne(d_dim)) {
      bool coarsen = false;
      IntVector tmp_rat = ratio_to_level_zero;
      for (int id = 0; id < d_dim.getValue(); id++) {
         if (ratio_to_level_zero(id) < 0) coarsen = true;
         tmp_rat(id) = abs(ratio_to_level_zero(id));
      }
      BoxSet tmp_mapped_boxes;
      if (coarsen) {
         domain_mapped_boxes.coarsen(tmp_mapped_boxes, tmp_rat);
      } else {
         domain_mapped_boxes.refine(tmp_mapped_boxes, tmp_rat);
      }
      domain_mapped_boxes.swap(tmp_mapped_boxes);
   }
}
#endif

void GridGeometry::computePhysicalDomain(
   BoxSet& domain_mapped_boxes,
   const IntVector& ratio_to_level_zero) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < d_dim.getValue(); i++) {
      TBOX_ASSERT(ratio_to_level_zero(i) != 0);
   }
   if (d_dim.getValue() > 1) {
      for (i = 0; i < d_dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % d_dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   for (int nb = 0; nb < d_number_blocks; nb++) {
      BoxSet block_domain_boxes = d_domain_with_images[nb];

      if (ratio_to_level_zero != IntVector::getOne(d_dim)) {
         bool coarsen = false;
         IntVector tmp_rat = ratio_to_level_zero;
         for (int id = 0; id < d_dim.getValue(); id++) {
            if (ratio_to_level_zero(id) < 0) coarsen = true;
            tmp_rat(id) = abs(ratio_to_level_zero(id));
         }
//         BoxSet tmp_mapped_boxes(d_dim);
         if (coarsen) {
            block_domain_boxes.coarsen(tmp_rat);
//            block_domain_boxes.coarsen(tmp_mapped_boxes, tmp_rat);
         } else {
            block_domain_boxes.refine(tmp_rat);
//            block_domain_boxes.refine(tmp_mapped_boxes, tmp_rat);
         }
//         block_domain_boxes.swap(tmp_mapped_boxes);
      }

      for (BoxSet::SetConstIterator bi = block_domain_boxes.setBegin();
           bi != block_domain_boxes.setEnd(); ++bi) {

         domain_mapped_boxes.insert(*bi);

      }
   }
}

/*
 *************************************************************************
 *
 * Set physical domain data member from input box array and determine
 * whether domain is a single box.
 *
 *************************************************************************
 */

void GridGeometry::setPhysicalDomain(
   const tbox::Array<BoxList>& domain)
{
   TBOX_ASSERT(domain.size() > 0);

   int number_blocks = domain.size();

   d_domain_is_single_box.resizeArray(number_blocks);
   d_physical_domain.resizeArray(number_blocks, BoxList(d_dim));
   d_domain_tree.resizeArray(number_blocks);
   d_domain_with_images.resizeArray(number_blocks, BoxSet(d_dim));
   d_number_blocks = number_blocks;

   for (int b = 0; b < number_blocks; b++) {

      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, domain[b]);

      BlockId block_id(b);
      BoxList domain_boxes(domain[b]);
      BoxList bounding_box(domain_boxes.getBoundingBox());

      bounding_box.removeIntersections(domain_boxes);
      if (bounding_box.size() == 0) {
         d_domain_is_single_box[b] = true;
         d_physical_domain[b] = BoxList(domain_boxes.getBoundingBox());
      } else {
         d_domain_is_single_box[b] = false;
         d_physical_domain[b] = domain[b];
      }

      resetDomainBoxSet(block_id);

      d_domain_tree[b] = new BoxTree(d_dim,
            d_physical_domain[b],
            BlockId(b));
   }
}

/*
 *************************************************************************
 *
 * Reset the domain BoxSet based on current definition of
 * physical domain and periodic shift.
 *
 *************************************************************************
 */

void GridGeometry::resetDomainBoxSet(const BlockId& block_id)
{
   const int& block_number = block_id.getBlockValue();
   BoxContainer domain_set(d_dim);
   BoxContainer::Iterator itr(d_physical_domain[block_number]);
   for (LocalId i(0); i < d_physical_domain[block_number].size(); ++i, itr++) {
      Box mapped_box(*itr, i, 0, BlockId(block_number));
      d_domain_with_images[block_number].insert(mapped_box);
   }
   d_physical_domain[block_number].clear();
   d_physical_domain[block_number] = d_domain_with_images[block_number];

   const IntVector periodic_shift_distances = getPeriodicShift(
         IntVector::getOne(d_dim));
   bool is_periodic = false;
   for (int i = 0; i < d_dim.getValue(); i++) {
      if (periodic_shift_distances(i) != 0) {
         is_periodic = true;
         break;
      }
   }
   if (is_periodic) {

      PeriodicShiftCatalog::initializeShiftsByIndexDirections(d_periodic_shift);
      const PeriodicShiftCatalog* periodic_shift_catalog =
         PeriodicShiftCatalog::getCatalog(d_dim);

      // Add periodic images of the domain d_domain_with_images.
      std::vector<IntVector> shifts;
      const BoxSet& mapped_boxes =
         d_domain_with_images[block_number];
      for (RealBoxConstIterator ni(mapped_boxes); ni.isValid(); ++ni) {
         const Box& mapped_box = *ni;
         computeShiftsForBox(shifts,
            mapped_box,
            *d_domain_tree[block_number],
            periodic_shift_distances);
         for (std::vector<IntVector>::const_iterator si = shifts.begin();
              si != shifts.end(); ++si) {
            PeriodicId shift_number =
               periodic_shift_catalog->shiftDistanceToShiftNumber(*si);
            Box image_mapped_box(
               mapped_box,
               shift_number,
               IntVector::getOne(d_dim));
            d_domain_with_images[block_number].insert(image_mapped_box);
         }
      }
   }
}

/*
 *************************************************************************
 *
 * The argument is an IntVector of length DIM.  It is set to 1
 * for periodic directions and 0 for all other directions.  In the
 * periodic directions, the coarse-level shift is calculated and stored
 * in the IntVector d_periodic_shift. The shift is the number of cells
 * in each periodic direction and is zero in all other directions.
 *
 *************************************************************************
 */

void GridGeometry::initializePeriodicShift(
   const IntVector& directions)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, directions);

   d_periodic_shift = directions;

   if (d_physical_domain.size() == 1) {

      int id;
      /*
       * Check incoming array and reset values if necessary.
       */
      for (id = 0; id < d_dim.getValue(); id++) {
         d_periodic_shift(id) = ((d_periodic_shift(id) == 0) ? 0 : 1);
      }

      /*
       * Check if the physical domain is valid for the specified
       * periodic conditions.  If so, compute the shift in each
       * dimension based on the the number of cells.
       */
      if (checkPeriodicValidity(d_physical_domain[0])) {

         BoxList domain_box_list(d_physical_domain[0]);
         Box bounding_box = domain_box_list.getBoundingBox();

         for (id = 0; id < d_dim.getValue(); id++) {
            d_periodic_shift(id) *= bounding_box.numberCells(id);
         }

      } else {
         TBOX_ERROR("Error in GridGeometry object with name = "
            << d_object_name << ": in intializePeriodicShift():  "
            << "Domain is not periodic for one (or more) of the dimensions "
            << "specified in the geometry input file!");
      }
   } else {
      TBOX_ASSERT(directions == IntVector::getZero(d_dim));
   }

   for (int b = 0; b < d_physical_domain.size(); b++) {
      resetDomainBoxSet(BlockId(b));
   }
}

/*
 *************************************************************************
 *
 * This returns an IntVector of length d_dim that is set to the width of
 * the domain in periodic directions and 0 in all other directions.
 * the argument contains the refinement ratio relative to the coarsest
 * level, which is multiplied by d_periodic_shift to get the return
 * vector.
 *
 *************************************************************************
 */

IntVector GridGeometry::getPeriodicShift(
   const IntVector& ratio_to_level_zero) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio vector must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   int k;
   for (k = 0; k < d_dim.getValue(); k++) {
      TBOX_ASSERT(ratio_to_level_zero(k) != 0);
   }
   if (d_dim.getValue() > 1) {
      for (k = 0; k < d_dim.getValue(); k++) {
         TBOX_ASSERT((ratio_to_level_zero(k)
                      * ratio_to_level_zero((k + 1) % d_dim.getValue()) > 0)
            || (ratio_to_level_zero(k) == 1)
            || (ratio_to_level_zero((k + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   IntVector periodic_shift(d_dim);
   for (int i = 0; i < d_dim.getValue(); i++) {
      if (ratio_to_level_zero(i) > 0) {
         periodic_shift(i) = d_periodic_shift(i) * ratio_to_level_zero(i);
      } else {
         int abs_ratio = abs(ratio_to_level_zero(i));
         periodic_shift(i) = d_periodic_shift(i) / abs_ratio;
      }
   }
   return periodic_shift;
}

/*
 *************************************************************************
 *
 * This checks if the periodic directions given to the constructor are
 * valid for the domain.  Periodic directions are valid if the domain
 * has exactly two physical boundaries normal to the periodic direction.
 *
 *************************************************************************
 */

bool GridGeometry::checkPeriodicValidity(
   const BoxList& domain)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, domain);

   bool is_valid = true;

   IntVector valid_direction(d_dim, 1);
   IntVector grow_direction(d_dim, 1);

   /*
    * Compute the bounding box of a "duplicate" domain + 1
    * cell and set the min and max indices of this grown box.
    */
   BoxList dup_domain(domain);

   Box domain_box = dup_domain.getBoundingBox();
   domain_box.grow(grow_direction);
   int i;
   Index min_index(d_dim, 0), max_index(d_dim, 0);
   for (i = 0; i < d_dim.getValue(); i++) {
      //set min/max of the bounding box
      min_index(i) = domain_box.lower(i);
      max_index(i) = domain_box.upper(i);
   }

   /*
    * Next, for each dimension, grow another "duplicate" domain
    * by 1.  Remove the intersections with the original domain,
    * and loop through the remaining box list, checking if the
    * upper index of the box matches the bounding box max or the
    * lower index of the box matches the bounding box min.  If
    * not, this dimension is not a valid periodic dimension.
    */
   for (i = 0; i < d_dim.getValue(); i++) {
      BoxList dup_domain2(domain);
      IntVector grow_one(d_dim, 0);
      grow_one(i) = 1;
      dup_domain2.grow(grow_one);
      dup_domain2.removeIntersections(BoxList(domain));

      BoxContainer::Iterator n(dup_domain2);
      for ( ; n; n++) {
         Box this_box = n();
         Index box_lower = this_box.lower();
         Index box_upper = this_box.upper();
         if (d_periodic_shift(i) != 0) {
            if (!((box_lower(i) == min_index(i)) ||
                  (box_upper(i) == max_index(i)))) {
               valid_direction(i) = 0;
            }
         }
      }
   }

   for (i = 0; i < d_dim.getValue(); i++) {
      if ((valid_direction(i) == 0) &&
          (d_periodic_shift(i) != 0)) {
         is_valid = false;
      }
   }

   return is_valid;
}

/*
 *************************************************************************
 *
 * Perform an error check on a recently-constructed boundary box to
 * make sure that it is the proper size, is adjacent to a patch, and is
 * outside the physical domain.
 *
 *************************************************************************
 */

bool GridGeometry::checkBoundaryBox(
   const BoundaryBox& boundary_box,
   const Patch& patch,
   const BoxList& domain,
   const int num_per_dirs,
   const IntVector& max_data_ghost_width) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS5(*this,
      boundary_box,
      patch,
      domain,
      max_data_ghost_width);

   bool return_val = true;

   const Box& bbox = boundary_box.getBox();

   /*
    * Test to see that the box is of size 1 in at least 1 direction.
    */
   IntVector box_size(d_dim);

   for (int i = 0; i < d_dim.getValue(); i++) {
      box_size(i) = bbox.numberCells(i);
   }

   if (box_size.min() != 1) {
      return_val = false;
   }

   /*
    * Quick and dirty test to see that boundary box is adjacent to patch
    * boundary, or a patch boundary extended through the ghost region.
    */
   Box patch_box = patch.getBox();

   Box grow_patch_box(patch_box);

   grow_patch_box.grow(IntVector::getOne(d_dim));

   if (!grow_patch_box.isSpatiallyEqual((grow_patch_box + bbox))) {
      bool valid_box = false;
      grow_patch_box = patch_box;
      for (int j = 0; j < d_dim.getValue(); j++) {
         if (num_per_dirs == 0) {

            for (int k = 1; k < d_dim.getValue(); k++) {

               grow_patch_box.grow((j + k) % d_dim.getValue(),
                  max_data_ghost_width((j + k) % d_dim.getValue()));

            }

         } else {

            for (int k = 1; k < d_dim.getValue(); k++) {

               grow_patch_box.grow((j + k) % d_dim.getValue(),
                  2 * max_data_ghost_width((j + k) % d_dim.getValue()));

            }

         }
         grow_patch_box.grow(j, 1);
         if (grow_patch_box.isSpatiallyEqual((grow_patch_box + bbox))) {
            valid_box = true;
         }
         grow_patch_box = patch_box;
      }
      if (!valid_box) {
         return_val = false;
      }
   }

   /*
    * check that the boundary box is outside the physical domain.
    */
   BoxList domain_list(domain);
   BoxList bbox_list(bbox);

   domain_list.intersectBoxes(bbox_list);

   if (domain_list.size()) {
      return_val = false;
   }

   return return_val;
}

/*
 ***************************************************************************
 * Read multiblock metadata from input database
 ***************************************************************************
 */
void GridGeometry::readBlockDataFromInput(
   const tbox::Pointer<tbox::Database>& input_db)
{
   TBOX_ASSERT(!input_db.isNull());

   d_singularity.resizeArray(d_number_blocks, BoxList(d_dim));
   d_singularity_indices.resizeArray(d_number_blocks);
   d_reduced_connect.resizeArray(d_number_blocks);
   d_block_neighbors.resizeArray(d_number_blocks);

   std::string sing_name;
   std::string neighbor_name;

   for (int i = 0; i < d_number_blocks; i++) {

      d_reduced_connect[i] = false;

   }

   for (d_number_of_block_singularities = 0; true; ++d_number_of_block_singularities) {

      sing_name = "Singularity" + tbox::Utilities::intToString(d_number_of_block_singularities);

      if (!input_db->keyExists(sing_name)) {
         break;
      }

      tbox::Pointer<tbox::Database> sing_db =
         input_db->getDatabase(sing_name);

      tbox::Array<int> blocks = sing_db->getIntegerArray("blocks");

      for (int i = 0; i < blocks.size(); i++) {

         const int block_number = blocks[i];

         std::string block_box_name = "sing_box_"
            + tbox::Utilities::intToString(block_number);

         Box sing_box(sing_db->getDatabaseBox(block_box_name));

         d_singularity[block_number].pushFront(sing_box);

         d_singularity_indices[block_number].push_back(d_number_of_block_singularities);
      }
   }

   for (int bn = 0; true; bn++) {
      neighbor_name = "BlockNeighbors" + tbox::Utilities::intToString(bn);

      if (!input_db->keyExists(neighbor_name)) {
         break;
      }
      tbox::Pointer<tbox::Database> pair_db =
         input_db->getDatabase(neighbor_name);

      BlockId block_a(pair_db->getInteger("block_a"));
      BlockId block_b(pair_db->getInteger("block_b"));
      Transformation::RotationIdentifier rotation_b_to_a;

      IntVector shift(d_dim, 0);
      if (d_dim.getValue() == 1) {
         rotation_b_to_a = Transformation::NO_ROTATE;
      } else {
         tbox::Array<std::string> rstr =
            pair_db->getStringArray("rotation_b_to_a");
         rotation_b_to_a = Transformation::getRotationIdentifier(rstr, d_dim);

         tbox::Array<int> b_array =
            pair_db->getIntegerArray("point_in_b_space");
         tbox::Array<int> a_array =
            pair_db->getIntegerArray("point_in_a_space");

         Index b_index(d_dim);
         Index a_index(d_dim);

         for (int p = 0; p < d_dim.getValue(); p++) {
            b_index(p) = b_array[p];
            a_index(p) = a_array[p];
         }

         Box b_box(b_index, b_index);
         Box a_box(a_index, a_index);

         b_box.rotate(rotation_b_to_a);
         Index b_rotated_point(b_box.lower());
         Index a_point = (a_box.lower());

         shift = a_point - b_rotated_point;
      }

      bool is_singularity =
         pair_db->getBoolWithDefault("is_singularity", false);

      registerNeighbors(block_a, block_b,
         rotation_b_to_a, shift, is_singularity);

   }

   if (d_number_blocks > 1) {
      for (int b = 0; b < d_number_blocks; b++) {
         BoxList pseudo_domain(d_dim);
         getDomainOutsideBlock(pseudo_domain, BlockId(b));

         hier::BoxList physical_domain(d_physical_domain[b]);
         pseudo_domain.spliceFront(physical_domain);

         for (BoxContainer::Iterator
              si(d_singularity[b]); si; si++) {
            BoxList test_domain(pseudo_domain);
            test_domain.intersectBoxes(si());
            if (test_domain.size() == 0) {
               d_reduced_connect[b] = true;
               break;
            }
         }
      }
   }
}

/*
 * ************************************************************************
 *
 * Get a BoxList representing all of the domain outside the given block.
 *
 * ************************************************************************
 */

void
GridGeometry::getDomainOutsideBlock(
   BoxList& domain_outside_block,
   const BlockId& block_id) const
{
   for (tbox::List<Neighbor>::Iterator
        nei(d_block_neighbors[block_id.getBlockValue()]); nei; nei++) {
      hier::BoxList transformed_domain(nei().getTransformedDomain()); 
      domain_outside_block.spliceFront(transformed_domain);
   }
}

/*
 * ************************************************************************
 *
 * Register a neighbor relationship between two blocks.
 *
 * ************************************************************************
 */

void GridGeometry::registerNeighbors(
   const BlockId& block_a,
   const BlockId& block_b,
   const Transformation::RotationIdentifier rotation,
   const IntVector& shift,
   const int is_singularity)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, shift);

   const int& a = block_a.getBlockValue();
   const int& b = block_b.getBlockValue();
   BoxList b_domain_in_a_space(d_physical_domain[b]);
   BoxList a_domain_in_b_space(d_physical_domain[a]);

   Transformation::RotationIdentifier back_rotation =
      Transformation::getReverseRotationIdentifier(rotation, d_dim);
   IntVector back_shift(d_dim);

   if (d_dim.getValue() == 2 || d_dim.getValue() == 3) {
      Transformation::calculateReverseShift(back_shift, shift, rotation);
   } else {
      TBOX_ERROR("GridGeometry::registerNeighbors error...\n"
         << "  object name = " << d_object_name
         << " Multiblock only works for 2D and 3D" << std::endl);
   }

   bool rotation_needed;
   if (rotation != 0) {
      rotation_needed = true;
   } else {
      rotation_needed = false;
   }

   if (rotation_needed) {
      b_domain_in_a_space.rotate(rotation);
      a_domain_in_b_space.rotate(back_rotation);
   }
   b_domain_in_a_space.shift(shift);
   a_domain_in_b_space.shift(back_shift);

   Transformation transformation(rotation, shift);
   Transformation back_transformation(back_rotation, back_shift);
   Neighbor neighbor_of_b(block_a, a_domain_in_b_space,
                          back_transformation,
                          is_singularity);
   Neighbor neighbor_of_a(block_b, b_domain_in_a_space,
                          transformation,
                          is_singularity);

   d_block_neighbors[a].addItem(neighbor_of_a);
   d_block_neighbors[b].addItem(neighbor_of_b);

   if (is_singularity) {
      d_has_enhanced_connectivity = true;
   }

}

/*
 *************************************************************************
 * Rotate and shift a box according to the rotation and shift that is
 * used to transform the index space of input_block into the
 * index space of base_block.
 *************************************************************************
 */

bool
GridGeometry::transformBox(
   Box& box,
   const IntVector& ratio,
   const BlockId& output_block,
   const BlockId& input_block) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, box, ratio);

   for (tbox::List<Neighbor>::Iterator
        ni(d_block_neighbors[output_block.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == input_block) {
         IntVector refined_shift = (ni().getShift()) * (ratio);
         box.rotate(ni().getRotationIdentifier());
         box.shift(refined_shift);
         return true;
      }
   }
   return false;
}

/*
 * ************************************************************************
 *
 * Rotate and shift the boxes in the given array according to the
 * rotation and shift that is used to transformed the index space of
 * input_block into the index space of output_block.
 *
 * ************************************************************************
 */

bool
GridGeometry::transformBoxList(
   BoxList& boxes,
   const IntVector& ratio,
   const BlockId& output_block,
   const BlockId& input_block) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, boxes, ratio);

   for (tbox::List<Neighbor>::Iterator
        ni(d_block_neighbors[output_block.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == input_block) {
         IntVector refined_shift = (ni().getShift()) * (ratio);
         boxes.rotate(ni().getRotationIdentifier());
         boxes.shift(refined_shift);
         return true;
      }
   }
   return false;
}

/*
 * ************************************************************************
 *
 * Set block to be the domain of transformed_block in the index space of
 * base_block.
 *
 * ************************************************************************
 */

void
GridGeometry::getTransformedBlock(
   BoxList& block,
   const BlockId& base_block,
   const BlockId& transformed_block)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, block);

   for (tbox::List<Neighbor>::Iterator
        ni(d_block_neighbors[base_block.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == transformed_block) {
         block = ni().getTransformedDomain();
         break;
      }
   }
}

/*
 * ************************************************************************
 *
 * Adjust all of the boundary boxes on the level so that they are
 * multiblock-aware.
 *
 * ************************************************************************
 */

void GridGeometry::adjustMultiblockPatchLevelBoundaries(
   PatchLevel& patch_level)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch_level);
   TBOX_ASSERT(patch_level.getGridGeometry()->getNumberBlocks() == d_number_blocks);

   if (d_number_blocks > 1) {

      const BoxSet& d_mapped_boxes =
         patch_level.getBoxLevel()->getBoxes();

      IntVector gcw(patch_level.getPatchDescriptor()->getMaxGhostWidth(d_dim));

      for (int nb = 0; nb < d_number_blocks; nb++) {

         const BlockId block_id(nb);

         BoxList singularity(d_singularity[nb]);
         singularity.refine(patch_level.getRatioToLevelZero());

         BoxList pseudo_domain(d_dim);

         for (tbox::List<GridGeometry::Neighbor>::Iterator
              nei(d_block_neighbors[nb]); nei; nei++) {
            hier::BoxList transformed_domain(nei().getTransformedDomain());
            pseudo_domain.spliceFront(transformed_domain);
         }

         pseudo_domain.refine(patch_level.getRatioToLevelZero());

         hier::BoxList physical_domain(patch_level.getPhysicalDomain(block_id));
         hier::BoxList sing_boxes(singularity); 
         pseudo_domain.spliceFront(physical_domain);
         pseudo_domain.spliceFront(sing_boxes);
         pseudo_domain.coalesce();

         BoxSetSingleBlockIterator mbi(d_mapped_boxes, block_id);

         for ( ; mbi.isValid(); mbi++) {
            const BoxId& mapped_box_id = (*mbi).getId();
            tbox::Pointer<Patch> patch(patch_level.getPatch(mapped_box_id));

            adjustBoundaryBoxesOnPatch(*patch,
               pseudo_domain,
               gcw,
               singularity);
         }
      }
   }
}

/*
 * ************************************************************************
 *
 * Adjust all of the boundary boxes on the patch so that they are
 * multiblock-aware.
 *
 * ************************************************************************
 */

void GridGeometry::adjustBoundaryBoxesOnPatch(
   const Patch& patch,
   const BoxList& pseudo_domain,
   const IntVector& gcw,
   const BoxList& singularity)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, patch, gcw);

   /*
    * Avoid adjusting boundary boxes for the case where we just use
    * a single block, since this is equivalent to not using multiblocks
    * at all.
    */
   if (d_number_blocks > 1) {
      PatchBoundaries boundaries(d_dim);

      getBoundaryBoxes(boundaries,
         patch.getBox(),
         BoxList(pseudo_domain),
         gcw,
         IntVector::getZero(d_dim));

      tbox::Array<BoundaryBox> codim_boundaries[tbox::Dimension::
                                                MAXIMUM_DIMENSION_VALUE];
      tbox::List<int> boundaries_in_sing[tbox::Dimension::
                                         MAXIMUM_DIMENSION_VALUE];
      for (int codim = 2; codim <= d_dim.getValue(); codim++) {

         codim_boundaries[codim - 1] =
            patch.getPatchGeometry()->getCodimensionBoundaries(codim);

         int num_boxes = codim_boundaries[codim - 1].size();

         for (int n = 0; n < num_boxes; n++) {
            Box border_box(codim_boundaries[codim - 1][n].getBox());
            BoxList sing_test_list(singularity);
            sing_test_list.intersectBoxes(border_box);
            if (sing_test_list.size() != 0) {
               boundaries_in_sing[codim - 1].addItem(n);
            }
         }
      }

      for (int i = 0; i < d_dim.getValue(); i++) {
         if (boundaries_in_sing[i].size() != 0) {
            int old_size = boundaries[i].size();
            boundaries[i].resizeArray(old_size + boundaries_in_sing[i].size(),
               BoundaryBox(d_dim));
            int nb = 0;
            for (tbox::List<int>::Iterator b(boundaries_in_sing[i]); b; b++) {
               boundaries[i][old_size + nb] = codim_boundaries[i][b()];
               boundaries[i][old_size + nb].setIsMultiblockSingularity(true);
               nb++;
            }
         }
         patch.getPatchGeometry()->setCodimensionBoundaries(boundaries[i],
            i + 1);
      }

   }

}

/*
 *************************************************************************
 *
 *************************************************************************
 */
Transformation::RotationIdentifier
GridGeometry::getRotationIdentifier(const BlockId& dst,
                                    const BlockId& src) const
{
   TBOX_ASSERT(areNeighbors(dst, src));

   Transformation::RotationIdentifier rotate = Transformation::NO_ROTATE;
   for (tbox::List<GridGeometry::Neighbor>::Iterator
        ni(d_block_neighbors[dst.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == src.getBlockValue()) {
         rotate = ni().getTransformation().getRotation();
         break;
      }
   }

   return rotate;
}

/*
 *************************************************************************
 *
 *************************************************************************
 */
const IntVector&
GridGeometry::getOffset(const BlockId& dst,
                        const BlockId& src) const
{
   TBOX_ASSERT(areNeighbors(dst, src));

   for (tbox::List<GridGeometry::Neighbor>::Iterator
        ni(d_block_neighbors[dst.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == src.getBlockValue()) {
         return ni().getTransformation().getOffset();
      }
   }

   return IntVector::getOne(d_dim);
}

/*
 *************************************************************************
 *
 *************************************************************************
 */
bool GridGeometry::areNeighbors(const BlockId& block_a,
                                const BlockId& block_b) const
{
   bool are_neighbors = false;

   for (tbox::List<GridGeometry::Neighbor>::Iterator
        ni(d_block_neighbors[block_a.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == block_b.getBlockValue()) {
         are_neighbors = true;
         break;
      }
   }

   return are_neighbors;
}

/*
 *************************************************************************
 *
 *************************************************************************
 */
bool GridGeometry::areSingularityNeighbors(const BlockId& block_a,
                                           const BlockId& block_b) const
{
   bool are_sing_neighbors = false;

   for (tbox::List<GridGeometry::Neighbor>::Iterator
        ni(d_block_neighbors[block_a.getBlockValue()]); ni; ni++) {
      if (ni().getBlockId() == block_b.getBlockValue()) {
         if (ni().isSingularity()) {
            are_sing_neighbors = true;
            break;
         }
      }
   }

   return are_sing_neighbors;
}

/*
 *************************************************************************
 *
 * Add operator to appropriate lookup list.
 *
 *************************************************************************
 */

void GridGeometry::addCoarsenOperator(
   tbox::Pointer<CoarsenOperator> coarsen_op)
{
   d_transfer_operator_registry->addCoarsenOperator(coarsen_op);
}

void GridGeometry::addRefineOperator(
   tbox::Pointer<RefineOperator> refine_op)
{
   d_transfer_operator_registry->addRefineOperator(refine_op);
}

void GridGeometry::addTimeInterpolateOperator(
   tbox::Pointer<TimeInterpolateOperator> time_op)
{
   d_transfer_operator_registry->addTimeInterpolateOperator(time_op);
}

/*
 *************************************************************************
 *
 * Search operator lists for operator matching request.
 *
 *************************************************************************
 */

tbox::Pointer<CoarsenOperator>
GridGeometry::lookupCoarsenOperator(
   const tbox::Pointer<Variable>& var,
   const std::string& op_name)
{
   return d_transfer_operator_registry->lookupCoarsenOperator(
      var, op_name);
}

tbox::Pointer<RefineOperator>
GridGeometry::lookupRefineOperator(
   const tbox::Pointer<Variable>& var,
   const std::string& op_name)
{
   return d_transfer_operator_registry->lookupRefineOperator(
      var, op_name);
}

tbox::Pointer<TimeInterpolateOperator>
GridGeometry::lookupTimeInterpolateOperator(
   const tbox::Pointer<Variable>& var,
   const std::string& op_name)
{
   return d_transfer_operator_registry->lookupTimeInterpolateOperator(
      var, op_name);
}

/*
 *************************************************************************
 * Compute the max operator stencil width from all constructed
 * refine and coarsen operators and the user-specified minimum value.
 *************************************************************************
 */
IntVector
GridGeometry::getMaxTransferOpStencilWidth()
{
   return d_transfer_operator_registry->getMaxTransferOpStencilWidth();
}

/*
 *************************************************************************
 * Set the mininum value to be returned by
 * getMaxTransferOpStencilWidth().
 *************************************************************************
 */

void
GridGeometry::setMinTransferOpStencilWidth(
   const IntVector& min_width)
{
   d_transfer_operator_registry->setMinTransferOpStencilWidth(min_width);
}

/*
 *************************************************************************
 *
 * Print object data to the specified output stream.
 *
 *************************************************************************
 */

void GridGeometry::printClassData(
   std::ostream& stream) const
{

   stream << "\nGridGeometry::printClassData..." << std::endl;
   stream << "GridGeometry: this = "
          << (GridGeometry *)this << std::endl;
   stream << "d_object_name = " << d_object_name << std::endl;

   for (int b = 0; b < d_physical_domain.size(); b++) {
      const int n = d_physical_domain[b].size();
      stream << "Block number " << b << std::endl;
      stream << "Number of boxes describing physical domain = " << n << std::endl;
      stream << "Boxes describing physical domain..." << std::endl;
      d_physical_domain[b].print(stream);
   }

   stream << "\nd_periodic_shift = " << d_periodic_shift << std::endl;

   stream << "d_max_data_ghost_width = " << d_max_data_ghost_width << std::endl;

   stream << "Block neighbor data:\n";

   for (int bn = 0; bn < d_number_blocks; ++bn) {

      stream << "   Block " << bn << '\n';

      const BlockId block_id(bn);
      const tbox::List<Neighbor>& block_neighbors(getNeighbors(block_id));

      const BoxList& singularity_boxlist(getSingularityBoxList(block_id));

      for (tbox::List<Neighbor>::Iterator li(block_neighbors); li; li++) {
         const Neighbor& neighbor(*li);
         stream << "      neighbor block " << neighbor.getBlockId() << ':';
         stream << " singularity = " << neighbor.isSingularity() << '\n';
      }

      stream << "      singularity Boxes (" << singularity_boxlist.size() << ")\n";
      for (BoxContainer::ConstIterator bi(singularity_boxlist); bi; bi++) {
         stream << "         " << *bi << '\n';
      }

   }

}

/*
 * ************************************************************************
 * ************************************************************************
 */

void GridGeometry::initializeCallback()
{
   t_find_patches_touching_boundaries = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::findPatchesTouchingBoundaries()");
   TBOX_ASSERT(!t_find_patches_touching_boundaries.isNull());
   t_touching_boundaries_init = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::...TouchingBoundaries()_init");
   TBOX_ASSERT(!t_touching_boundaries_init.isNull());
   t_touching_boundaries_loop = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::...TouchingBoundaries()_loop");
   TBOX_ASSERT(!t_touching_boundaries_loop.isNull());
   t_set_geometry_on_patches = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::setGeometryOnPatches()");
   TBOX_ASSERT(!t_set_geometry_on_patches.isNull());
   t_set_boundary_boxes = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::setBoundaryBoxes()");
   TBOX_ASSERT(!t_set_boundary_boxes.isNull());
   t_set_geometry_data_on_patches = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::set_geometry_data_on_patches");
   TBOX_ASSERT(!t_set_geometry_data_on_patches.isNull());
   t_compute_boundary_boxes_on_level = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::computeBoundaryBoxesOnLevel()");
   TBOX_ASSERT(!t_compute_boundary_boxes_on_level.isNull());
   t_get_boundary_boxes = tbox::TimerManager::getManager()->
      getTimer("hier::GridGeometry::getBoundaryBoxes()");
   TBOX_ASSERT(!t_get_boundary_boxes.isNull());
}

/*
 *************************************************************************
 *************************************************************************
 */

void GridGeometry::finalizeCallback()
{
   t_find_patches_touching_boundaries.setNull();
   t_touching_boundaries_init.setNull();
   t_touching_boundaries_loop.setNull();
   t_set_geometry_on_patches.setNull();
   t_set_boundary_boxes.setNull();
   t_set_geometry_data_on_patches.setNull();
   t_compute_boundary_boxes_on_level.setNull();
   t_get_boundary_boxes.setNull();
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
