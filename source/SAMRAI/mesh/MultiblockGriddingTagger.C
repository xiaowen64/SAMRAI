/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_mesh_MultiblockGriddingTagger_C
#define included_mesh_MultiblockGriddingTagger_C

#include "SAMRAI/mesh/MultiblockGriddingTagger.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *
 * The default constructor and virtual destructor do nothing
 * particularly interesting.
 *
 *************************************************************************
 */

MultiblockGriddingTagger::MultiblockGriddingTagger(
   const tbox::Dimension& dim):
   xfer::RefinePatchStrategy(dim),
   d_dim(dim)
{
}

MultiblockGriddingTagger::~MultiblockGriddingTagger()
{
}

void MultiblockGriddingTagger::setScratchTagPatchDataIndex(
   int buf_tag_indx)
{

   tbox::Pointer<hier::Variable> check_var;
   bool indx_maps_to_variable =
      hier::VariableDatabase::getDatabase()->mapIndexToVariable(buf_tag_indx,
         check_var);
   if (!indx_maps_to_variable || check_var.isNull()) {
      TBOX_ERROR(
         "MultiblockGriddingTagger::setScratchTagPatchDataIndex error...\n"
         << "Given patch data index = " << buf_tag_indx
         << " is not in VariableDatabase."
         << std::endl);
   } else {
      tbox::Pointer<pdat::CellVariable<int> > t_check_var = check_var;
      if (t_check_var.isNull()) {
         TBOX_ERROR(
            "MultiblockGriddingTagger::setScratchTagPatchDataIndex error...\n"
            << "Given patch data index = " << buf_tag_indx
            << " does not map to cell-centered"
            << "\ninteger data in VariableDatabase." << std::endl);
      }
   }

   d_buf_tag_indx = buf_tag_indx;
}

void MultiblockGriddingTagger::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, patch);

   NULL_USE(fill_time);

   const tbox::Pointer<pdat::CellData<int> > tag_data =
      patch.getPatchData(d_buf_tag_indx);

   hier::IntVector gcw =
      hier::IntVector::min(ghost_width_to_fill,
         tag_data->getGhostCellWidth());

   tbox::Pointer<hier::PatchGeometry> pgeom = patch.getPatchGeometry();

   for (int d = 0; d < d_dim.getValue(); d++) {

      tbox::Array<hier::BoundaryBox> bbox =
         pgeom->getCodimensionBoundaries(d + 1);

      for (int b = 0; b < bbox.size(); b++) {
         if (!bbox[b].getIsMultiblockSingularity()) {
            hier::Box fill_box = pgeom->getBoundaryFillBox(bbox[b],
                  patch.getBox(),
                  gcw);

            tag_data->fillAll(0, fill_box);
         }
      }
   }
}

void MultiblockGriddingTagger::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   const hier::Connector& dst_to_encon,
   const double fill_time,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box,
   const tbox::Pointer<hier::GridGeometry>& grid_geometry)
{
   NULL_USE(fill_time);
   NULL_USE(boundary_box);
   NULL_USE(grid_geometry);

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim, patch, fill_box, boundary_box);

   const tbox::Dimension& dim = fill_box.getDim();

   const hier::BoxId& dst_mb_id = patch.getBox().getId();

   const hier::BlockId& patch_blk_id = dst_mb_id.getBlockId();

   const tbox::Pointer<pdat::CellData<int> > tag_data =
      patch.getPatchData(d_buf_tag_indx);

   hier::Box sing_fill_box(tag_data->getGhostBox() * fill_box);
   tag_data->fillAll(0, sing_fill_box);

   if (grid_geometry->hasEnhancedConnectivity()) {

      const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
         grid_geometry->getNeighbors(patch_blk_id);

      const hier::NeighborhoodSet& dst_to_encon_nbrhood_set =
         dst_to_encon.getNeighborhoodSets();

      hier::NeighborhoodSet::const_iterator ni =
         dst_to_encon_nbrhood_set.find(dst_mb_id);

      if (ni != dst_to_encon_nbrhood_set.end()) {

         const hier::BoxSet& encon_nbrs = ni->second;

         for (hier::BoxSet::OrderedConstIterator ei = encon_nbrs.orderedBegin();
              ei != encon_nbrs.orderedEnd(); ++ei) {

            tbox::Pointer<hier::Patch> encon_patch(
               encon_level.getPatch(ei->getId()));

            const hier::BlockId& encon_blk_id = ei->getBlockId();

            hier::Transformation::RotationIdentifier rotation =
               hier::Transformation::NO_ROTATE;
            hier::IntVector offset(dim);

            for (tbox::List<hier::GridGeometry::Neighbor>::Iterator
                 ni(neighbors); ni; ni++) {

               if (ni().getBlockId() == encon_blk_id) {
                  rotation = ni().getRotationIdentifier();
                  offset = ni().getShift();
                  break;
               }
            }

            offset *= patch.getPatchGeometry()->getRatio();

            hier::Transformation transformation(rotation, offset);
            hier::Box encon_patch_box(encon_patch->getBox());
            transformation.transform(encon_patch_box);

            hier::Box encon_fill_box(encon_patch_box * sing_fill_box);
            if (!encon_fill_box.empty()) {

               const hier::Transformation::RotationIdentifier back_rotate =
                  hier::Transformation::getReverseRotationIdentifier(
                     rotation, dim);

               hier::IntVector back_shift(dim);

               hier::Transformation::calculateReverseShift(
                  back_shift, offset, rotation);

               hier::Transformation back_trans(back_rotate, back_shift);

               tbox::Pointer<pdat::CellData<int> > sing_data(
                  encon_patch->getPatchData(d_buf_tag_indx));

               for (pdat::CellIterator ci(encon_fill_box); ci; ci++) {
                  pdat::CellIndex src_index(ci());
                  pdat::CellGeometry::transform(src_index, back_trans);

                  (*tag_data)(ci()) += (*sing_data)(src_index);
               }

            }
         }
      }
   }
}

void MultiblockGriddingTagger::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(d_dim, fine, coarse, fine_box, ratio);

   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
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
