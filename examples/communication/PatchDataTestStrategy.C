//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/communication/PatchDataTestStrategy.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Base class for patch data test operations.
//

#include "PatchDataTestStrategy.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CellData.h"
#include "PatchLevel.h"
#include "tbox/Utilities.h"

namespace SAMRAI {

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

PatchDataTestStrategy::PatchDataTestStrategy()
{
   d_variable_src_name.resizeArray(0);
   d_variable_dst_name.resizeArray(0);
   d_variable_depth.resizeArray(0);
   d_variable_src_ghosts.resizeArray(0);
   d_variable_dst_ghosts.resizeArray(0);
   d_variable_coarsen_op.resizeArray(0);
   d_variable_refine_op.resizeArray(0);
}

PatchDataTestStrategy::~PatchDataTestStrategy()
{
}

/*
*************************************************************************
*                                                                       *
* Routines for reading variable and refinement data from input.         *
*                                                                       *
*************************************************************************
*/

void PatchDataTestStrategy::readVariableInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
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

void PatchDataTestStrategy::readRefinementInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
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
* Blank physical boundary and pre/postprocess coarsen/refine operations *
* so tester isn't required to implement them when not needed.           *
*                                                                       *
*************************************************************************
*/

void PatchDataTestStrategy::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   (void) patch;
   (void) time;
   (void) gcw;
}

void PatchDataTestStrategy::preprocessRefine(hier::Patch<NDIM>& fine,
                                             const hier::Patch<NDIM>& coarse,
                                             const hier::Box<NDIM>& fine_box,
                                             const hier::IntVector<NDIM>& ratio) const
{
   (void) fine;
   (void) coarse;
   (void) fine_box;
   (void) ratio;
}

void PatchDataTestStrategy::postprocessRefine(hier::Patch<NDIM>& fine,
                                              const hier::Patch<NDIM>& coarse,
                                              const hier::Box<NDIM>& fine_box,
                                              const hier::IntVector<NDIM>& ratio) const
{
   (void) fine;
   (void) coarse;
   (void) fine_box;
   (void) ratio;
}

void PatchDataTestStrategy::preprocessCoarsen(hier::Patch<NDIM>& coarse,
                                              const hier::Patch<NDIM>& fine,
                                              const hier::Box<NDIM>& coarse_box,
                                              const hier::IntVector<NDIM>& ratio) const
{
   (void) coarse;
   (void) fine;
   (void) coarse_box;
   (void) ratio;
}

void PatchDataTestStrategy::postprocessCoarsen(hier::Patch<NDIM>& coarse,
                                               const hier::Patch<NDIM>& fine,
                                               const hier::Box<NDIM>& coarse_box,
                                               const hier::IntVector<NDIM>& ratio) const
{
   (void) coarse;
   (void) fine;
   (void) coarse_box;
   (void) ratio;
}

}
