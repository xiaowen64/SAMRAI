//
// File:        PatchMultiblockTestStrategy.C
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 1.6 $
// Modified:    $Date: 2004/02/11 23:46:08 $
// Description: Base class for patch data test operations.
//

#include "PatchMultiblockTestStrategy.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CellData.h"
#include "PatchLevel.h"
#include "tbox/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

using namespace SAMRAI;

// These are used in the cell tagging routine.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/*
*************************************************************************
*									*
* The constructor and destructor.                                       *
*									*
*************************************************************************
*/

PatchMultiblockTestStrategy::PatchMultiblockTestStrategy()
{
   d_variable_src_name.resizeArray(0);
   d_variable_dst_name.resizeArray(0);
   d_variable_depth.resizeArray(0);
   d_variable_src_ghosts.resizeArray(0);
   d_variable_dst_ghosts.resizeArray(0);
   d_variable_coarsen_op.resizeArray(0);
   d_variable_refine_op.resizeArray(0);
}

PatchMultiblockTestStrategy::~PatchMultiblockTestStrategy()
{
   if (d_grid_geometry) delete[] d_grid_geometry;
}

/*
*************************************************************************
*                                                                       *
* Routines for reading variable and refinement data from input.         *
*                                                                       *
*************************************************************************
*/

void PatchMultiblockTestStrategy::readVariableInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   tbox::Array<string> var_keys = db->getAllKeys();
   int nkeys = var_keys.getSize();

   d_variable_src_name.resizeArray(nkeys);
   d_variable_dst_name.resizeArray(nkeys);
   d_variable_depth.resizeArray(nkeys);
   d_variable_src_ghosts.resizeArray(nkeys);
   d_variable_dst_ghosts.resizeArray(nkeys);
   d_variable_coarsen_op.resizeArray(nkeys);
   d_variable_refine_op.resizeArray(nkeys);

   for (int i = 0; i < nkeys; i++) {

      tbox::Pointer<tbox::Database> var_db = db->getDatabase(var_keys[i]);

      if (var_db->keyExists("src_name")) {
         d_variable_src_name[i] = var_db->getString("src_name");
      } else {
         TBOX_ERROR("Variable input error: No `src_name' string found for "
                    << "key = " << var_keys[i] << endl);
      }

      if (var_db->keyExists("dst_name")) {
         d_variable_dst_name[i] = var_db->getString("dst_name");
      } else {
         TBOX_ERROR("Variable input error: No `dst_name' string found for "
                    << "key = " << var_keys[i] << endl);
      }

      if (var_db->keyExists("depth")) {
         d_variable_depth[i] = var_db->getInteger("depth");
      } else {
         d_variable_depth[i] = 1; 
      }

      if (var_db->keyExists("src_ghosts")) {
         int* tmp_ghosts = d_variable_src_ghosts[i]; 
         var_db->getIntegerArray("src_ghosts", tmp_ghosts, NDIM);
      } else {
         d_variable_src_ghosts[i] = hier::IntVector<NDIM>(0); 
      }

      if (var_db->keyExists("dst_ghosts")) {
         int* tmp_ghosts = d_variable_dst_ghosts[i]; 
         var_db->getIntegerArray("dst_ghosts", tmp_ghosts, NDIM);
      } else {
         d_variable_dst_ghosts[i] = hier::IntVector<NDIM>(0); 
      }

      if (var_db->keyExists("coarsen_operator")) {
         d_variable_coarsen_op[i] = var_db->getString("coarsen_operator");
      } else {
         d_variable_coarsen_op[i] = "NO_COARSEN"; 
      }

      if (var_db->keyExists("refine_operator")) {
         d_variable_refine_op[i] = var_db->getString("refine_operator");
      } else {
         d_variable_refine_op[i] = "NO_REFINE"; 
      }

   }

}

void PatchMultiblockTestStrategy::readRefinementInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   tbox::Array<string> box_keys = db->getAllKeys();
   int nkeys = box_keys.getSize();

   d_refine_level_boxes.resizeArray(nkeys);
   for (int i = 0; i < nkeys; i++) {
      d_refine_level_boxes[i] = db->getDatabaseBoxArray(box_keys[i]); 
   }

}

/*
*************************************************************************
*                                                                       *
* Tag cells on level specified in input box array for refinement.       *
*                                                                       *
*************************************************************************
*/

void PatchMultiblockTestStrategy::tagCellsInInputBoxes(
   hier::Patch<NDIM>& patch, 
   int level_number, 
   int tag_index)
{

   if (level_number < d_refine_level_boxes.getSize()) { 

      tbox::Pointer< pdat::CellData<NDIM,int> > tags =
         patch.getPatchData(tag_index);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!tags.isNull());
#endif
      tags->fillAll(0);

      const hier::Box<NDIM> pbox = patch.getBox();

      hier::BoxArray<NDIM> boxes = d_refine_level_boxes[level_number];
      int nboxes = boxes.getNumberOfBoxes();
      for (int k = 0; k < nboxes; k++) {
         tags->fill(1, boxes.getBox(k)*pbox, 0);
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Blank physical boundary and pre/postprocess coarsen/refine operations *
* so tester isn't required to implement them when not needed.           *
*                                                                       *
*************************************************************************
*/

void PatchMultiblockTestStrategy::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   (void) patch;
   (void) time;
   (void) gcw;
}

void PatchMultiblockTestStrategy::preprocessRefine(
   hier::Patch<NDIM>& fine,
   const hier::Patch<NDIM>& coarse,
   const tbox::Pointer<hier::VariableContext>& context,
   const hier::Box<NDIM>& fine_box,
   const hier::IntVector<NDIM>& ratio) const
{
   (void) fine;
   (void) coarse;
   (void) context;
   (void) fine_box;
   (void) ratio;
}

void PatchMultiblockTestStrategy::postprocessRefine(hier::Patch<NDIM>& fine,
   const hier::Patch<NDIM>& coarse,
   const tbox::Pointer<hier::VariableContext>& context,
   const hier::Box<NDIM>& fine_box,
   const hier::IntVector<NDIM>& ratio) const
{
   (void) fine;
   (void) coarse;
   (void) context;
   (void) fine_box;
   (void) ratio;
}

void PatchMultiblockTestStrategy::preprocessCoarsen(hier::Patch<NDIM>& coarse,
   const hier::Patch<NDIM>& fine,
   const tbox::Pointer<hier::VariableContext>& context,
   const hier::Box<NDIM>& coarse_box,
   const hier::IntVector<NDIM>& ratio) const
{
   (void) coarse;
   (void) fine;
   (void) context;
   (void) coarse_box;
   (void) ratio;
}

void PatchMultiblockTestStrategy::postprocessCoarsen(hier::Patch<NDIM>& coarse,
   const hier::Patch<NDIM>& fine,
   const tbox::Pointer<hier::VariableContext>& context,
   const hier::Box<NDIM>& coarse_box,
   const hier::IntVector<NDIM>& ratio) const
{
   (void) coarse;
   (void) fine;
   (void) context;
   (void) coarse_box;
   (void) ratio;
}

