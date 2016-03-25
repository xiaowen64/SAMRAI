//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/multiblock/MultiblockSideDataTranslator.C $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1168 $
// Modified:	$LastChangedDate: 2006-07-11 16:29:55 -0700 (Tue, 11 Jul 2006) $
// Description:	Templated operations for copying patch data.
//

#ifndef included_pdat_MultiblockSideDataTranslator_C
#define included_pdat_MultiblockSideDataTranslator_C

#include "MultiblockSideDataTranslator.h"

#include "SideData.h"


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
MultiblockSideDataTranslator<DIM,TYPE>::MultiblockSideDataTranslator()
{
}

template <int DIM, class TYPE>
MultiblockSideDataTranslator<DIM,TYPE>::~MultiblockSideDataTranslator()
{
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for side data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MultiblockSideDataTranslator<DIM,TYPE>::translateAndCopyData(
   hier::Patch<DIM>& dst_patch,
   const int dst_id,
   const hier::Patch<DIM>& src_patch,
   const int src_id,
   const hier::IntVector<DIM>& shift,
   const typename
      hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
 
   tbox::Pointer< SideData<DIM,TYPE> > dst = dst_patch.getPatchData(dst_id);
   tbox::Pointer< SideData<DIM,TYPE> > src = src_patch.getPatchData(src_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(dst.isNull()));
   TBOX_ASSERT(!(src.isNull()));
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
#endif

   hier::IntVector<DIM> dir_vector(dst->getDirectionVector());
   if (rotate == 0) {
      for (int axis = 0; axis < DIM; axis++) {
         if (dir_vector(axis)) {
            translateAndCopyArrayData(dst->getArrayData(axis),
                                      src->getArrayData(axis),
                                      shift,
                                      rotate);
         }
      }
   } else if (DIM == 2) {
      for (int axis = 0; axis < DIM; axis++) {
         if (dir_vector(axis)) {
            for (pdat::SideIterator<DIM> fi(dst->getBox(),axis); fi; fi++) {
               pdat::SideIndex<DIM> dst_index(fi());
               hier::Index<DIM> dst_xyz_index(dst_index);

               hier::Index<DIM> src_xyz_index(dst_xyz_index);
   
               int num_rotations = (4-rotate)%4;
               hier::IntVector<DIM> copy_shift(shift);

               int src_axis;
               if (num_rotations%2) {
                  src_axis = (axis+1)%DIM;
               } else {
                  src_axis = axis;
               }

               for (int r = 0; r < num_rotations; r++) {
                  hier::Index<DIM> tmp_index(src_xyz_index);
                  src_xyz_index(0) = tmp_index(1);
                  src_xyz_index(1) = -tmp_index(0)-1;
                  hier::IntVector<DIM> tmp_shift(copy_shift);
                  copy_shift(0) = tmp_shift(1);
                  copy_shift(1) = -tmp_shift(0);
               }

               for (int i=0; i<DIM; i++) {
                  src_xyz_index(i) -= copy_shift(i);
               }

               pdat::SideIndex<DIM> src_index;
               if (src_axis == 1) {
                  src_index(0) = src_xyz_index(0);
                  src_index(1) = src_xyz_index(1);
                  if (num_rotations == 1 || num_rotations == 2) {
                     src_index(1)++;
                  }
               } else {
                  src_index(0) = src_xyz_index(0);
                  if ((num_rotations == 3) || (num_rotations == 2)) {
                     src_index(0)++;
                  }
                  src_index(1) = src_xyz_index(1);
               }
               src_index.setAxis(src_axis);
   
               for (int d = 0; d < dst->getDepth(); d++) {
                  (*dst)(dst_index, d) = (*src)(src_index, d);
               }
            }
         }
      }
   } else if (DIM == 3) {
      for (int axis = 0; axis < DIM; axis++) {
         if (dir_vector(axis)) {
            int src_axis;
            if (axis == 0) {
               switch (rotate) {

                  case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:

                     src_axis = 0;
                     break;

                  case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:

                     src_axis = 1;
                     break;

                  default:

                     src_axis = 2;
                     break;

               }
            } else if (axis == 1) {

               switch (rotate) {
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:

                     src_axis = 0;
                     break;

                  case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:

                     src_axis = 1;
                     break;

                  default:

                     src_axis = 2;
                     break;

               }

            } else {

               switch (rotate) {
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:
                     src_axis = 0;
                     break;

                  case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
                  case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
                  case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:

                     src_axis = 1;
                     break;

                  default:

                     src_axis = 2;
                     break;

               }
            }

            for (pdat::SideIterator<DIM> si(dst->getBox(),axis); si; si++) {
               pdat::SideIndex<DIM> dst_index(si());
               hier::Index<DIM> dst_xyz_index(dst_index);

               typename
               hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier back_rotate =
                  hier::MultiblockPatchHierarchy<DIM>::
                     getReverseRotationIdentifier(rotate);

               hier::Box<DIM> src_box(dst_xyz_index, dst_xyz_index);

               src_box.rotate(back_rotate);

               hier::IntVector<DIM> back_shift;
               hier::MultiblockPatchHierarchy<DIM>::calculateReverseShift(
                  back_shift, shift, back_rotate);

               src_box.shift(back_shift);

               hier::Index<DIM> src_xyz_index;
               for (int i = 0; i < DIM; i++) {
                  src_xyz_index(i) = src_box.lower()(i);
               }

               pdat::SideIndex<DIM> src_index;
               src_index(0) = src_xyz_index(0);
               src_index(1) = src_xyz_index(1);
               src_index(2) = src_xyz_index(2);

               if (axis == 0) {
                  switch (rotate) {

                     case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
                        break;

                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:
                        src_index(src_axis)++;
                        break;

                     default:
                        TBOX_ERROR(" ");
                        break;
                  }
               } else if (axis == 1) {
                  switch (rotate) {

                     case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
                        break;
                                                                                
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:
                        src_index(src_axis)++;
                        break;
                                                                                
                     default:
                        TBOX_ERROR(" ");
                        break;
                  }
               } else {
                  switch (rotate) {

                     case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
                        break;

                     case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:
                     case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:
                        src_index(src_axis)++;
                        break;

                     default:
                        TBOX_ERROR(" ");
                        break;
                  }
               }

               src_index.setAxis(src_axis);

               for (int d = 0; d < dst->getDepth(); d++) {
                  (*dst)(dst_index, d) = (*src)(src_index, d);
               }
            }
         }
      }
   }
}


/*
*************************************************************************
*                                                                       *
* Translation and copy for array data                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void MultiblockSideDataTranslator<DIM,TYPE>::translateAndCopyArrayData(
   pdat::ArrayData<DIM,TYPE>& dst,
   const pdat::ArrayData<DIM,TYPE>& src,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   bool no_rotate = true;
   if (rotate != 0) {
      no_rotate = false;
   }

   if (no_rotate) {
      dst.copy(src, dst.getBox(), shift);
   } else if (DIM < 3) {
      hier::Box<DIM> rotatebox(src.getBox());
      int num_rotations = rotate;

      rotatebox.rotate(num_rotations);

      const hier::Box<DIM> copybox = dst.getBox() *
                                hier::Box<DIM>::shift(rotatebox, shift);

      if (!copybox.empty()) {
         TYPE *const dst_ptr = dst.getPointer();
         const TYPE *const src_ptr = src.getPointer();

         const int depth = (dst.getDepth() < src.getDepth() ?
                            dst.getDepth() : src.getDepth());

         const int box_w0 = copybox.numberCells(0);

         const int dst_w0 = dst.getBox().numberCells(0);
         int src_w0 = src.getBox().numberCells(0);
         if (num_rotations == 3) {
            src_w0 = -src_w0;
         }
         const int box_w1 = copybox.numberCells(1);

         const int dst_offset = dst.getOffset();
         const int src_offset = src.getOffset();

         int dst_bd = dst.getBox().offset(copybox.lower());
         hier::Index<DIM> src_index(copybox.lower()-shift);

         // rotate src_index 4-num_rotations;
         for (int r = 0; r < 4-num_rotations; r++) {
            hier::Index<DIM> tmp_index(src_index);
            src_index(0) = tmp_index(1);
            src_index(1) = -tmp_index(0)-1;
         }

         int src_bd_orig = src.getBox().offset(src_index);

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
   } else {
      TBOX_ERROR("MultiblockSideDataTranslator<DIM,TYPE>::translateAndCopyArrayData : DIM = 1 or > 3 not implemented");

   }
}

}
}
#endif
