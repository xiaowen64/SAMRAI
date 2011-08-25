/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for copying patch data.
 *
 ************************************************************************/

#ifndef included_pdat_MultiblockNodeDataTranslator_C
#define included_pdat_MultiblockNodeDataTranslator_C

#include "SAMRAI/pdat/MultiblockNodeDataTranslator.h"

#include "SAMRAI/pdat/NodeData.h"

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
MultiblockNodeDataTranslator<TYPE>::MultiblockNodeDataTranslator()
{
}

template<class TYPE>
MultiblockNodeDataTranslator<TYPE>::~MultiblockNodeDataTranslator()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Translation and copy for node data                                    *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MultiblockNodeDataTranslator<TYPE>::translateAndCopyData(
   hier::Patch& dst_patch,
   const int dst_id,
   const hier::Patch& src_patch,
   const int src_id,
   const hier::IntVector& shift,
   const hier::Transformation::RotationIdentifier rotate)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(dst_patch, src_patch, shift);

   const tbox::Dimension& dim(dst_patch.getDim());

   tbox::Pointer<NodeData<TYPE> > dst = dst_patch.getPatchData(dst_id);
   tbox::Pointer<NodeData<TYPE> > src = src_patch.getPatchData(src_id);

   TBOX_ASSERT(!(dst.isNull()));
   TBOX_ASSERT(!(src.isNull()));

   bool no_rotate;
   if (rotate != 0) {
      no_rotate = false;
   } else {
      no_rotate = true;
   }

   if (no_rotate) {
      dst->getArrayData().copy(src->getArrayData(),
         dst->getArrayData().getBox(), shift);
   } else if (dim == tbox::Dimension(2)) {
      hier::Box rotatebox(src->getArrayData().getBox());
      int num_rotations = rotate;

      rotatebox.rotate(rotate);

      const hier::Box copybox = dst->getArrayData().getBox()
         * (pdat::NodeGeometry::toNodeBox(
               hier::Box::shift(rotatebox, shift)));

      if (!copybox.empty()) {
         TYPE * const dst_ptr = dst->getArrayData().getPointer();
         const TYPE * const src_ptr = src->getArrayData().getPointer();

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
         hier::Index src_index(copybox.lower() - shift);

         // rotate src_index 4-num_rotations;
         for (int r = 0; r < 4 - num_rotations; r++) {
            hier::Index tmp_index(src_index);
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
   } else if (dim == tbox::Dimension(3)) {
      for (pdat::NodeIterator ni(dst->getBox()); ni; ni++) {
         pdat::NodeIndex dst_index(ni());

         hier::Transformation::RotationIdentifier back_rotate =
            hier::Transformation::getReverseRotationIdentifier(rotate, dim);

         hier::Box src_box(dst_index, dst_index);

         src_box.rotate(back_rotate);

         hier::IntVector back_shift(dim);
         hier::Transformation::calculateReverseShift(
            back_shift, shift, rotate);

         src_box.shift(back_shift);

         pdat::NodeIndex src_index(dim);
         for (int n = 0; n < dim.getValue(); n++) {
            src_index(n) = src_box.lower() (n);
         }

         switch (rotate) {

            case hier::Transformation::IUP_JUP_KUP:
            case hier::Transformation::KUP_IUP_JUP:
            case hier::Transformation::JUP_KUP_IUP:
               break;

            case hier::Transformation::IDOWN_KUP_JUP:
            case hier::Transformation::KUP_JUP_IDOWN:
            case hier::Transformation::JUP_IDOWN_KUP:
               src_index(0)++;
               break;

            case hier::Transformation::KDOWN_JUP_IUP:
            case hier::Transformation::IUP_KDOWN_JUP:
            case hier::Transformation::JUP_IUP_KDOWN:
               src_index(2)++;
               break;

            case hier::Transformation::KDOWN_IDOWN_JUP:
            case hier::Transformation::IDOWN_JUP_KDOWN:
            case hier::Transformation::JUP_KDOWN_IDOWN:
               src_index(0)++;
               src_index(2)++;
               break;

            case hier::Transformation::JDOWN_IUP_KUP:
            case hier::Transformation::IUP_KUP_JDOWN:
            case hier::Transformation::KUP_JDOWN_IUP:
               src_index(1)++;
               break;

            case hier::Transformation::JDOWN_KUP_IDOWN:
            case hier::Transformation::IDOWN_JDOWN_KUP:
            case hier::Transformation::KUP_IDOWN_JDOWN:
               src_index(0)++;
               src_index(1)++;
               break;

            case hier::Transformation::JDOWN_KDOWN_IUP:
            case hier::Transformation::KDOWN_IUP_JDOWN:
            case hier::Transformation::IUP_JDOWN_KDOWN:
               src_index(1)++;
               src_index(2)++;
               break;

            case hier::Transformation::JDOWN_IDOWN_KDOWN:
            case hier::Transformation::KDOWN_JDOWN_IDOWN:
            case hier::Transformation::IDOWN_KDOWN_JDOWN:
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
   } else {
      TBOX_ERROR(
         "MultiblockNodeDataTranslator<TYPE>::translateAndCopyData : dim = 1 or > 3 not implemented");
   }
}

}
}
#endif
