//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/multiblock/MBUtilities.C $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2009 $
// Modified:	$LastChangedDate: 2008-02-26 15:38:52 -0800 (Tue, 26 Feb 2008) $
// Description:	utility functions for multiblock
//

#ifndef included_hier_MBUtilities_C
#define included_hier_MBUtilities_C

#include "MBUtilities.h"

#include "MultiblockDataTranslator.h"
#include "tbox/ShutdownRegistry.h"




namespace SAMRAI {
    namespace hier {

#ifndef NULL
#defined NULL (0)
#endif

/*
*************************************************************************
*                                                                       *
* Constructor and destructor do nothing, as all member functions in     *
* this class are static.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM> MBUtilities<DIM>::MBUtilities()
{
}

template<int DIM> MBUtilities<DIM>::~MBUtilities()
{
}

/*
*************************************************************************
*                                                                       *
* Determines the patch data type and calls the appropriate routine      *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MBUtilities<DIM>::translateAndCopyData(
   hier::Patch<DIM>& dst_patch,
   const int dst_id,
   const hier::Patch<DIM>& src_patch,
   const int src_id,
   const hier::IntVector<DIM>& shift,
   const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   tbox::Pointer< PatchDataFactory<DIM> > dst_pdf = 
      VariableDatabase<DIM>::getDatabase()->getPatchDescriptor()->
         getPatchDataFactory(dst_id);

   MultiblockDataTranslator<DIM>* mb_trans =
      dst_pdf->getMultiblockDataTranslator();

   mb_trans->translateAndCopyData(dst_patch,
                                  dst_id,
                                  src_patch,
                                  src_id,
                                  shift,
                                  rotate);
}

template<int DIM>
void MBUtilities<DIM>::translateAndFillData(
   hier::Patch<DIM>& dst_patch,
   const int dst_id,
   const hier::Patch<DIM>& src_patch,
   const int src_id,
   const hier::IntVector<DIM>& shift,
   const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   tbox::Pointer< PatchDataFactory<DIM> > dst_pdf =
      VariableDatabase<DIM>::getDatabase()->getPatchDescriptor()->
         getPatchDataFactory(dst_id);
                                                                                
   MultiblockDataTranslator<DIM>* mb_trans =
      dst_pdf->getMultiblockDataTranslator();
                                                                                
   mb_trans->translateAndFillData(dst_patch,
                                  dst_id,
                                  src_patch,
                                  src_id,
                                  shift,
                                  rotate);
}

/*
*************************************************************************
*                                                                       *
* rotate an index around the origin.                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MBUtilities<DIM>::rotateIndex(
   int* index,
   const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotation)
{
   if (DIM == 2) {
      int num_rotations = (int) rotation;

      for (int j = 0; j < num_rotations; j++) {
         int tmp_in[DIM];
         tmp_in[0] = index[0];
         tmp_in[1] = index[1];

         index[0] = tmp_in[1];
         index[1] = -tmp_in[0]-1;
      }
   } else if (DIM == 3) {
      if (rotation == MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP) {
         return;
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP) {
         rotateAboutAxis(index,0,3);
         rotateAboutAxis(index,2,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP) {
         rotateAboutAxis(index,1,1);
         rotateAboutAxis(index,2,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP) {
         rotateAboutAxis(index,1,2);
         rotateAboutAxis(index,0,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN) {
         rotateAboutAxis(index,1,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP) {
         rotateAboutAxis(index,2,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP) {
         rotateAboutAxis(index,1,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP) {
         rotateAboutAxis(index,0,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN) {
         rotateAboutAxis(index,0,2);
         rotateAboutAxis(index,2,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP) {
         rotateAboutAxis(index,0,3);
         rotateAboutAxis(index,2,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN) {
         rotateAboutAxis(index,1,2);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN) {
         rotateAboutAxis(index,0,3);
         rotateAboutAxis(index,1,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP) {
         rotateAboutAxis(index,2,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN) {
         rotateAboutAxis(index,0,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP) {
         rotateAboutAxis(index,0,2);
         rotateAboutAxis(index,1,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN) {
         rotateAboutAxis(index,0,1);
         rotateAboutAxis(index,1,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP) {
         rotateAboutAxis(index,0,2);
         rotateAboutAxis(index,1,2);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN) {
         rotateAboutAxis(index,0,1);
         rotateAboutAxis(index,2,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP) {
         rotateAboutAxis(index,0,3);
         rotateAboutAxis(index,1,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN) {
         rotateAboutAxis(index,0,1);
         rotateAboutAxis(index,2,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN) {
         rotateAboutAxis(index,0,2);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN) {
         rotateAboutAxis(index,0,2);
         rotateAboutAxis(index,2,1);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN) {
         rotateAboutAxis(index,0,2);
         rotateAboutAxis(index,1,3);
      } else if (rotation == MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN) {
         rotateAboutAxis(index,1,2);
         rotateAboutAxis(index,0,1);
      }
   } else {
      TBOX_ERROR("MBUtilities<DIM>::rotateIndex : DIM = 1 or > 3 not implemented");
   }

}

/*
*************************************************************************
*                                                                       *
* Private routine to rotate an index about an axis.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MBUtilities<DIM>::rotateAboutAxis(int* index,
                                       const int axis,
                                       const int num_rotations)
{
   if (DIM == 3) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(axis < DIM);
#endif

      const int a = (axis+1)%DIM;
      const int b = (axis+2)%DIM;

      for (int j = 0; j < num_rotations; j++) {
         int tmp_in[3] = {index[0], index[1], index[2]};
         index[a] = tmp_in[b];
         index[b] = -tmp_in[a]-1;
      }
   }
}
   
}
}
#endif
