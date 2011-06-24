/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data. 
 *
 ************************************************************************/

#ifndef included_mesh_MultiblockGriddingTagger_C
#define included_mesh_MultiblockGriddingTagger_C

#include "SAMRAI/mesh/MultiblockGriddingTagger.h"

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
 *									*
 * The default constructor and virtual destructor do nothing             *
 * particularly interesting.		                                *
 *									*
 *************************************************************************
 */

MultiblockGriddingTagger::MultiblockGriddingTagger(
   const tbox::Dimension& dim):
   xfer::MultiblockRefinePatchStrategy(dim),
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
   tbox::List<tbox::Pointer<hier::Patch> >& singularity_patches,
   const double fill_time,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim, patch, fill_box, boundary_box);

   NULL_USE(boundary_box);
   NULL_USE(fill_time);

   const tbox::Dimension dim(fill_box.getDim());

   const tbox::Pointer<pdat::CellData<int> > tag_data =
      patch.getPatchData(d_buf_tag_indx);

   int num_sing_patches = singularity_patches.getNumberOfItems();
   tbox::Pointer<pdat::CellData<int> >* sing_tag_data;
   sing_tag_data =
      new tbox::Pointer<pdat::CellData<int> >[num_sing_patches];

   int sn = 0;

   for (tbox::List<tbox::Pointer<hier::Patch> >::Iterator
        sp(singularity_patches); sp; sp++) {
      sing_tag_data[sn] =
         sp()->getPatchData(d_buf_tag_indx);
      sn++;
   }

   hier::Box sing_box(dim);
   if (num_sing_patches == 0) {
      sing_box = fill_box;
   } else {
      sing_box = sing_tag_data[0]->getBox();
      for (int i = 1; i < num_sing_patches; i++) {
         sing_box = sing_box * sing_tag_data[i]->getBox();
      }
   }

   const hier::Box tag_fill_box =
      fill_box * sing_box * tag_data->getGhostBox();

   for (pdat::CellData<int>::Iterator cdi(tag_fill_box);
        cdi; cdi++) {
      for (int n = 0; n < num_sing_patches; n++) {
         int sing_tag_val = (*sing_tag_data[n])(cdi());
         if (sing_tag_val != 0) {
            (*tag_data)(cdi()) = sing_tag_val;
            break;
         }
      }
   }

   if (num_sing_patches > 1) {
      for (int st = 0; st < num_sing_patches; st++) {
         if ((sing_tag_data[st]->getBox() + sing_box) != sing_box) {
            hier::Box sing_tag_box(sing_tag_data[st]->getBox());

            hier::Box new_fill_box = fill_box * sing_tag_box
               * tag_data->getGhostBox();
            for (pdat::CellData<int>::Iterator cdi(new_fill_box);
                 cdi; cdi++) {
               if (!sing_box.contains(cdi())) {
                  int sing_tag_val = (*sing_tag_data[st])(cdi());
                  if (sing_tag_val != 0) {
                     (*tag_data)(cdi()) = sing_tag_val;
                  }
               }
            }
         }
      }
   }

   delete[] sing_tag_data;
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
