//
// File:	MultiblockGriddingTagger.C
// Package:	SAMRAI multiblock package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 455 $
// Modified:	$Date: 2005-06-16 12:45:23 -0700 (Thu, 16 Jun 2005) $
// Description:	Strategy interface to user routines for refining AMR data.
//

#ifndef included_mblk_MultiblockGriddingTagger_C
#define included_mblk_MultiblockGriddingTagger_C 

#include "MultiblockGriddingTagger.h"

#include "CellData.h"
#include "CellVariable.h"
#include "tbox/Utilities.h"

namespace SAMRAI {
    namespace mblk {

/*
*************************************************************************
*									*
* The default constructor and virtual destructor do nothing             *
* particularly interesting.		                                *
*									*
*************************************************************************
*/

template<int DIM>
MultiblockGriddingTagger<DIM>::MultiblockGriddingTagger()
{
}

template<int DIM>
MultiblockGriddingTagger<DIM>::~MultiblockGriddingTagger()
{
}

template<int DIM>
void MultiblockGriddingTagger<DIM>::setScratchTagPatchDataIndex(int buf_tag_indx)
{

   tbox::Pointer< hier::Variable<DIM> > check_var;
   bool indx_maps_to_variable = 
      hier::VariableDatabase<DIM>::getDatabase()->mapIndexToVariable(buf_tag_indx, check_var);
   if (!indx_maps_to_variable || check_var.isNull()) {
      TBOX_ERROR("MultiblockGriddingTagger<DIM>::setScratchTagPatchDataIndex error...\n"
                 << "Given patch data index = " << buf_tag_indx << " is not in VariableDatabase." 
                 << endl);
   } else {
      tbox::Pointer< pdat::CellVariable<DIM,int> > t_check_var = check_var;
      if (t_check_var.isNull()) {
         TBOX_ERROR("MultiblockGriddingTagger<DIM>::setScratchTagPatchDataIndex error...\n"
                    << "Given patch data index = " << buf_tag_indx << " does not map to cell-centered"
                    << "\ninteger data in VariableDatabase." << endl);
      }
   }

   d_buf_tag_indx = buf_tag_indx;
}

template<int DIM>
void MultiblockGriddingTagger<DIM>::setPhysicalBoundaryConditions(
   hier::Patch<DIM>& patch,
   const double fill_time,
   const hier::IntVector<DIM>& ghost_width_to_fill)
{
   (void) fill_time;

   const tbox::Pointer< pdat::CellData<DIM,int> > tag_data =
      patch.getPatchData(d_buf_tag_indx);

   hier::IntVector<DIM> gcw =
      hier::IntVector<DIM>::min(ghost_width_to_fill,
                                tag_data->getGhostCellWidth());

   tbox::Pointer< hier::PatchGeometry<DIM> > pgeom = patch.getPatchGeometry();

   for (int d = 0; d < DIM; d++) {

      tbox::Array< hier::BoundaryBox<DIM> > bbox =
         pgeom->getCodimensionBoundary(d+1);

      for (int b = 0; b < bbox.size(); b++) {
         if (!bbox[b].getIsMultiblockSingularity()) {
            hier::Box<DIM> fill_box = pgeom->getBoundaryFillBox(bbox[b],
                                                           patch.getBox(),
                                                           gcw);

            tag_data->fillAll(0, fill_box);
         }
      }
   }
}

template<int DIM>
void MultiblockGriddingTagger<DIM>::fillSingularityBoundaryConditions(
   hier::Patch<DIM>& patch,
   tbox::List<typename MultiblockRefineSchedule<DIM>::SingularityPatch>&
      singularity_patches,
   const double fill_time,
   const hier::Box<DIM>& fill_box,
   const hier::BoundaryBox<DIM>& boundary_box)
{
   (void) boundary_box;

   const tbox::Pointer< pdat::CellData<DIM,int> > tag_data =
      patch.getPatchData(d_buf_tag_indx);

   int num_sing_patches = singularity_patches.getNumberItems();
      tbox::Pointer< pdat::CellData<DIM,int> > *sing_tag_data;
   sing_tag_data =
      new tbox::Pointer< pdat::CellData<DIM,int> >[num_sing_patches];

   int sn = 0;

   for (typename tbox::List
        <typename MultiblockRefineSchedule<DIM>::SingularityPatch>::Iterator
        sp(singularity_patches); sp; sp++) {
      sing_tag_data[sn] =
         sp().d_patch->getPatchData(d_buf_tag_indx);
         sn++;
   }

   hier::Box<DIM> sing_box;
   if (num_sing_patches == 0) {
      sing_box = fill_box;
   } else {
      sing_box = sing_tag_data[0]->getBox();
      for (int i = 1; i < num_sing_patches; i++) {
         sing_box = sing_box * sing_tag_data[i]->getBox();
      }
   } 

   const hier::Box<DIM> tag_fill_box =
      fill_box * sing_box * tag_data->getGhostBox();

   for (typename pdat::CellData<DIM,int>::Iterator cdi(tag_fill_box);
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
            hier::Box<DIM> sing_tag_box(sing_tag_data[st]->getBox());

            hier::Box<DIM> new_fill_box = fill_box * sing_tag_box *
                                     tag_data->getGhostBox();
            for (typename pdat::CellData<DIM,int>::Iterator cdi(new_fill_box);
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

template<int DIM>
void MultiblockGriddingTagger<DIM>::postprocessRefine(
   hier::Patch<DIM>& fine,
   const hier::Patch<DIM>& coarse,
   const hier::Box<DIM>& fine_box,
   const hier::IntVector<DIM>& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}



}
}
#endif
