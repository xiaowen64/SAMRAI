//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/multiblock/MultiblockNodeDataTranslator.C $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1168 $
// Modified:	$LastChangedDate: 2006-07-11 16:29:55 -0700 (Tue, 11 Jul 2006) $
// Description:	Templated operations for copying patch data.
//

#ifndef included_pdat_MultiblockNodeDataTranslator_C
#define included_pdat_MultiblockNodeDataTranslator_C

#include "MultiblockNodeDataTranslator.h"

#include "NodeData.h"


namespace SAMRAI {
    namespace pdat {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor do nothing, as all member functions in     *
* this class are static.                                                *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
MultiblockNodeDataTranslator<DIM,TYPE>::MultiblockNodeDataTranslator()
{
}

template <int DIM, class TYPE>
MultiblockNodeDataTranslator<DIM,TYPE>::~MultiblockNodeDataTranslator()
{
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for node data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MultiblockNodeDataTranslator<DIM,TYPE>::translateAndCopyData(
   hier::Patch<DIM>& dst_patch,
   const int dst_id,
   const hier::Patch<DIM>& src_patch,
   const int src_id,
   const hier::IntVector<DIM>& shift,
   const typename
      hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{

   tbox::Pointer< NodeData<DIM,TYPE> > dst = dst_patch.getPatchData(dst_id);
   tbox::Pointer< NodeData<DIM,TYPE> > src = src_patch.getPatchData(src_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(dst.isNull()));
   TBOX_ASSERT(!(src.isNull()));
#endif

   bool no_rotate;
   if (rotate != 0) {
      no_rotate = false;
   } else {
      no_rotate = true;
   }
     
   if (no_rotate) {
      dst->getArrayData().copy(src->getArrayData(), 
                              dst->getArrayData().getBox(), shift);
   } else if (DIM == 2) {
      hier::Box<DIM> rotatebox(src->getArrayData().getBox());
      int num_rotations = rotate;

      rotatebox.rotate(num_rotations);

      const hier::Box<DIM> copybox = dst->getArrayData().getBox() *
                                (pdat::NodeGeometry<DIM>::toNodeBox(
                                   hier::Box<DIM>::shift(rotatebox, shift)));

      if (!copybox.empty()) {
         TYPE *const dst_ptr = dst->getArrayData().getPointer();
         const TYPE *const src_ptr = src->getArrayData().getPointer();

         const int depth = (dst->getArrayData().getDepth() < 
                            src->getArrayData().getDepth() ?
                            dst->getArrayData().getDepth() : 
                            src->getArrayData().getDepth());

         const int box_w0 = copybox.numberCells(0);

         const int dst_w0 = dst->getArrayData().getBox().numberCells(0);
         int src_w0 = src->getArrayData().getBox().numberCells(0);
         if (num_rotations == 3) {
            src_w0 = -src_w0;
         }
         const int box_w1 = copybox.numberCells(1);

         const int dst_offset = dst->getArrayData().getOffset();
         const int src_offset = src->getArrayData().getOffset();

         int dst_bd = dst->getArrayData().getBox().offset(copybox.lower());
         hier::Index<DIM> src_index(copybox.lower()-shift);

         // rotate src_index 4-num_rotations;
         for (int r = 0; r < 4-num_rotations; r++) {
            hier::Index<DIM> tmp_index(src_index);
            src_index(0) = tmp_index(1);
            src_index(1) = -tmp_index(0);
         }

         int src_bd_orig = src->getArrayData().getBox().offset(src_index);

         for (int d = 0; d < depth; d++) {
            int src_bd = src_bd_orig;
            int dst_b2 = dst_bd;
            int src_b2 = src_bd;
            int dst_b1 = dst_b2;
            int src_b1 = src_b2;

            for (int i1 = 0; i1 < box_w1; i1++) {

               for (int i0 = 0; i0 < box_w0; i0++) {
                  if (i0) {
                     if (num_rotations%2) {
                        src_b1 += src_w0;
                     } else {
                        src_b1--;
                     }
                  }
                  dst_ptr[dst_b1+i0] = src_ptr[src_b1];
               }
               dst_b1 += dst_w0;
               if (num_rotations == 1) {
                  src_b1 = --src_bd;
               } else if (num_rotations == 2) {
                  src_b1 = src_bd - src_w0;
                  src_bd = src_b1;
               } else if (num_rotations == 3) {
                  src_b1 = ++src_bd;
               }
            }

            dst_bd += dst_offset;
            src_bd_orig += src_offset;
         }
      }
   } else if (DIM == 3) {
      for (pdat::NodeIterator<DIM> ni(dst->getBox()); ni; ni++) {
         pdat::NodeIndex<DIM> dst_index(ni());

         typename
         hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier back_rotate =
         hier::MultiblockPatchHierarchy<DIM>::getReverseRotationIdentifier(rotate);

         hier::Box<DIM> src_box(dst_index, dst_index);

         src_box.rotate(back_rotate);

         hier::IntVector<DIM> back_shift;
         hier::MultiblockPatchHierarchy<DIM>::calculateReverseShift(
            back_shift, shift, back_rotate);

         src_box.shift(back_shift);

         pdat::NodeIndex<DIM> src_index;
         for (int n = 0; n < DIM; n++) {
            src_index(n) = src_box.lower()(n);
         }

         switch (rotate) {

            case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
            case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
            case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
               break;

            case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
            case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
               src_index(0)++;
               break;

            case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
            case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
            case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
               src_index(2)++;
               break;

            case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
            case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
               src_index(0)++;
               src_index(2)++;
               break;

            case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
            case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
               src_index(1)++;
               break;
 
            case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
            case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
               src_index(0)++;
               src_index(1)++;
               break;
 
            case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
            case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
               src_index(1)++;
               src_index(2)++;
               break;
 
            case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:
            case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:
               src_index(0)++;
               src_index(1)++;
               src_index(2)++;
               break;

            default:
               TBOX_ERROR(" ");
               break;
         }

         //back rotate src_index into src index space
         //back shift src_index to needed location
         //copy data from src_index to dst_index

         for (int d = 0; d < dst->getDepth(); d++) {
            (*dst)(dst_index, d) = (*src)(src_index, d);
         }
      }
   }
}

}
}
#endif
