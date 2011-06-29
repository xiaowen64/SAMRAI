/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for copying patch data. 
 *
 ************************************************************************/

#ifndef included_pdat_MultiblockFaceDataTranslator_C
#define included_pdat_MultiblockFaceDataTranslator_C

#include "SAMRAI/pdat/MultiblockFaceDataTranslator.h"

#include "SAMRAI/pdat/FaceData.h"

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

template<class TYPE>
MultiblockFaceDataTranslator<TYPE>::MultiblockFaceDataTranslator()
{
}

template<class TYPE>
MultiblockFaceDataTranslator<TYPE>::~MultiblockFaceDataTranslator()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Translation and copy for face data                                    *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MultiblockFaceDataTranslator<TYPE>::translateAndCopyData(
   hier::Patch& dst_patch,
   const int dst_id,
   const hier::Patch& src_patch,
   const int src_id,
   const hier::IntVector& shift,
   const hier::Transformation::RotationIdentifier rotate)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(dst_patch, src_patch, shift);

   const tbox::Dimension& dim(shift.getDim());

   tbox::Pointer<FaceData<TYPE> > dst = dst_patch.getPatchData(dst_id);
   tbox::Pointer<FaceData<TYPE> > src = src_patch.getPatchData(src_id);

   TBOX_ASSERT(!(dst.isNull()));
   TBOX_ASSERT(!(src.isNull()));

   if (rotate == 0) {
      for (int axis = 0; axis < dim.getValue(); axis++) {
         hier::Index face_shift(dim);
         for (int d = 0; d < dim.getValue(); d++) {
            face_shift(d) = shift((axis + d) % dim.getValue());
         }
         translateAndCopyArrayData(dst->getArrayData(axis),
            src->getArrayData(axis),
            face_shift,
            rotate);
      }
   } else if (dim == tbox::Dimension(2)) {
      for (int axis = 0; axis < dim.getValue(); axis++) {
         for (pdat::FaceIterator fi(dst->getBox(), axis); fi; fi++) {
            pdat::FaceIndex dst_index(fi());
            hier::Index dst_xyz_index(dst_index);
            if (axis == 1) {
               dst_xyz_index(0) = dst_index(1);
               dst_xyz_index(1) = dst_index(0);
            }

            hier::Index src_xyz_index(dst_xyz_index);

            int num_rotations = (4 - rotate) % 4;
            hier::IntVector copy_shift(shift);

            int src_axis;
            if (num_rotations % 2) {
               src_axis = (axis + 1) % dim.getValue();
            } else {
               src_axis = axis;
            }

            for (int r = 0; r < num_rotations; r++) {
               hier::Index tmp_index(src_xyz_index);
               src_xyz_index(0) = tmp_index(1);
               src_xyz_index(1) = -tmp_index(0) - 1;
               hier::IntVector tmp_shift(copy_shift);
               copy_shift(0) = tmp_shift(1);
               copy_shift(1) = -tmp_shift(0);
            }

            for (int i = 0; i < dim.getValue(); i++) {
               src_xyz_index(i) -= copy_shift(i);
            }

            pdat::FaceIndex src_index(dim);
            if (src_axis == 1) {
               src_index(0) = src_xyz_index(1);
               if (num_rotations == 1 || num_rotations == 2) {
                  src_index(0)++;
               }
               src_index(1) = src_xyz_index(0);
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
   } else if (dim == tbox::Dimension(3)) {
      for (int axis = 0; axis < dim.getValue(); axis++) {
         int src_axis;
         if (axis == 0) {
            switch (rotate) {

               case hier::Transformation::IUP_JUP_KUP:
               case hier::Transformation::IDOWN_KUP_JUP:
               case hier::Transformation::IUP_KDOWN_JUP:
               case hier::Transformation::IDOWN_JUP_KDOWN:
               case hier::Transformation::IUP_KUP_JDOWN:
               case hier::Transformation::IDOWN_JDOWN_KUP:
               case hier::Transformation::IUP_JDOWN_KDOWN:
               case hier::Transformation::IDOWN_KDOWN_JDOWN:

                  src_axis = 0;
                  break;

               case hier::Transformation::JUP_KUP_IUP:
               case hier::Transformation::JUP_IDOWN_KUP:
               case hier::Transformation::JUP_IUP_KDOWN:
               case hier::Transformation::JUP_KDOWN_IDOWN:
               case hier::Transformation::JDOWN_IUP_KUP:
               case hier::Transformation::JDOWN_KUP_IDOWN:
               case hier::Transformation::JDOWN_KDOWN_IUP:
               case hier::Transformation::JDOWN_IDOWN_KDOWN:

                  src_axis = 1;
                  break;

               default:

                  src_axis = 2;
                  break;

            }
         } else if (axis == 1) {
            switch (rotate) {
               case hier::Transformation::KUP_IUP_JUP:
               case hier::Transformation::JUP_IDOWN_KUP:
               case hier::Transformation::JUP_IUP_KDOWN:
               case hier::Transformation::KDOWN_IDOWN_JUP:
               case hier::Transformation::JDOWN_IUP_KUP:
               case hier::Transformation::KUP_IDOWN_JDOWN:
               case hier::Transformation::KDOWN_IUP_JDOWN:
               case hier::Transformation::JDOWN_IDOWN_KDOWN:

                  src_axis = 0;
                  break;

               case hier::Transformation::IUP_JUP_KUP:
               case hier::Transformation::KUP_JUP_IDOWN:
               case hier::Transformation::KDOWN_JUP_IUP:
               case hier::Transformation::IDOWN_JUP_KDOWN:
               case hier::Transformation::KUP_JDOWN_IUP:
               case hier::Transformation::IDOWN_JDOWN_KUP:
               case hier::Transformation::IUP_JDOWN_KDOWN:
               case hier::Transformation::KDOWN_JDOWN_IDOWN:

                  src_axis = 1;
                  break;

               default:

                  src_axis = 2;
                  break;

            }

         } else {

            switch (rotate) {
               case hier::Transformation::JUP_KUP_IUP:
               case hier::Transformation::KUP_JUP_IDOWN:
               case hier::Transformation::KDOWN_JUP_IUP:
               case hier::Transformation::JUP_KDOWN_IDOWN:
               case hier::Transformation::KUP_JDOWN_IUP:
               case hier::Transformation::JDOWN_KUP_IDOWN:
               case hier::Transformation::JDOWN_KDOWN_IUP:
               case hier::Transformation::KDOWN_JDOWN_IDOWN:
                  src_axis = 0;
                  break;

               case hier::Transformation::KUP_IUP_JUP:
               case hier::Transformation::IDOWN_KUP_JUP:
               case hier::Transformation::IUP_KDOWN_JUP:
               case hier::Transformation::KDOWN_IDOWN_JUP:
               case hier::Transformation::IUP_KUP_JDOWN:
               case hier::Transformation::KUP_IDOWN_JDOWN:
               case hier::Transformation::KDOWN_IUP_JDOWN:
               case hier::Transformation::IDOWN_KDOWN_JDOWN:

                  src_axis = 1;
                  break;

               default:

                  src_axis = 2;
                  break;

            }
         }
         for (pdat::FaceIterator fi(dst->getBox(), axis); fi; fi++) {
            pdat::FaceIndex dst_index(fi());
            hier::Index dst_xyz_index(dst_index);
            if (axis == 1) {
               dst_xyz_index(0) = dst_index(2);
               dst_xyz_index(1) = dst_index(0);
               dst_xyz_index(2) = dst_index(1);
            } else if (axis == 2) {
               dst_xyz_index(0) = dst_index(1);
               dst_xyz_index(1) = dst_index(2);
               dst_xyz_index(2) = dst_index(0);
            }

            hier::Transformation::RotationIdentifier back_rotate =
               hier::Transformation::
               getReverseRotationIdentifier(rotate, dim);

            hier::Box src_box(dst_xyz_index, dst_xyz_index);

            src_box.rotate(back_rotate);

            hier::IntVector back_shift(dim);
            hier::Transformation::calculateReverseShift(
               back_shift, shift, rotate);

            src_box.shift(back_shift);

            hier::Index src_xyz_index(dim);
            for (int i = 0; i < dim.getValue(); i++) {
               src_xyz_index(i) = src_box.lower() (i);
            }

            pdat::FaceIndex src_index(dim);
            if (src_axis == 0) {
               src_index(0) = src_xyz_index(0);
               src_index(1) = src_xyz_index(1);
               src_index(2) = src_xyz_index(2);
            } else if (src_axis == 1) {
               src_index(0) = src_xyz_index(1);
               src_index(1) = src_xyz_index(2);
               src_index(2) = src_xyz_index(0);
            } else {
               src_index(0) = src_xyz_index(2);
               src_index(1) = src_xyz_index(0);
               src_index(2) = src_xyz_index(1);
            }

            if (axis == 0) {
               switch (rotate) {

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                     break;

                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:
                     src_index(0)++;
                     break;

                  default:
                     TBOX_ERROR(" ");
                     break;
               }
            } else if (axis == 1) {
               switch (rotate) {

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                     break;

                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:
                     src_index(0)++;
                     break;

                  default:
                     TBOX_ERROR(" ");
                     break;
               }
            } else {
               switch (rotate) {

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                     break;

                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:
                     src_index(0)++;
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

/*
 *************************************************************************
 *                                                                       *
 * Translation and copy for array data                                   *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MultiblockFaceDataTranslator<TYPE>::translateAndCopyArrayData(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const hier::IntVector& shift,
   const hier::Transformation::RotationIdentifier rotate)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(dst, src, shift);

   const tbox::Dimension& dim(dst.getDim());

   bool no_rotate = true;
   if (rotate != 0) {
      no_rotate = false;
   }

   if (no_rotate) {
      dst.copy(src, dst.getBox(), shift);
   } else if (dim < tbox::Dimension(3)) {
      hier::Box rotatebox(src.getBox());
      int num_rotations = rotate;

      rotatebox.rotate(rotate);

      const hier::Box copybox = dst.getBox()
         * hier::Box::shift(rotatebox, shift);

      if (!copybox.empty()) {
         TYPE * const dst_ptr = dst.getPointer();
         const TYPE * const src_ptr = src.getPointer();

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
         hier::Index src_index(copybox.lower() - shift);

         // rotate src_index 4-num_rotations;
         for (int r = 0; r < 4 - num_rotations; r++) {
            hier::Index tmp_index(src_index);
            src_index(0) = tmp_index(1);
            src_index(1) = -tmp_index(0) - 1;
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
                     if (num_rotations % 2) {
                        src_b1 += src_w0;
                     } else {
                        src_b1--;
                     }
                  }
                  dst_ptr[dst_b1 + i0] = src_ptr[src_b1];
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
      TBOX_ERROR(
         "MultiblockFaceDataTranslator<TYPE>::translateAndCopyData : dim = 1 or > 3 not implemented");
   }
}

}
}
#endif
