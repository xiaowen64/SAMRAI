/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for copying patch data.
 *
 ************************************************************************/

#ifndef included_pdat_MultiblockCellDataTranslator_C
#define included_pdat_MultiblockCellDataTranslator_C

#include "SAMRAI/pdat/MultiblockCellDataTranslator.h"

#include "SAMRAI/pdat/CellData.h"

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
MultiblockCellDataTranslator<TYPE>::MultiblockCellDataTranslator()
{
}

template<class TYPE>
MultiblockCellDataTranslator<TYPE>::~MultiblockCellDataTranslator()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Translation and copy for cell data                                    *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void MultiblockCellDataTranslator<TYPE>::translateAndCopyData(
   hier::Patch& dst_patch,
   const int dst_id,
   const hier::Patch& src_patch,
   const int src_id,
   const hier::IntVector& shift,
   const hier::Transformation::RotationIdentifier rotate)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(dst_patch, src_patch, shift);

   const tbox::Dimension& dim(shift.getDim());

   tbox::Pointer<CellData<TYPE> > dst = dst_patch.getPatchData(dst_id);
   tbox::Pointer<CellData<TYPE> > src = src_patch.getPatchData(src_id);

   TBOX_ASSERT(!(dst.isNull()));
   TBOX_ASSERT(!(src.isNull()));

   if (dim == tbox::Dimension(2)) {
      translateAndCopyArrayData(dst->getArrayData(),
         src->getArrayData(),
         shift,
         rotate);
   } else if (dim == tbox::Dimension(3)) {
      if (rotate == 0) {
         translateAndCopyArrayData(dst->getArrayData(),
            src->getArrayData(),
            shift,
            rotate);
      } else {

         const hier::Transformation::RotationIdentifier back_rotate =
            hier::Transformation::getReverseRotationIdentifier(
               rotate, dim);

         hier::IntVector back_shift(dim);

         hier::Transformation::calculateReverseShift(
            back_shift, shift, rotate);

         TYPE * const dst_ptr = dst->getArrayData().getPointer();
         const TYPE * const src_ptr = src->getArrayData().getPointer();

         const int depth = ((dst->getDepth() < src->getDepth()) ?
                            dst->getDepth() : src->getDepth());

         const hier::Box dst_box = dst->getBox();
         const hier::Box src_box = src->getBox();
         const hier::Box dst_ghost_box = dst->getGhostBox();
         const hier::Box src_ghost_box = src->getGhostBox();

         int dst_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         int dst_w[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         int src_glo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         int dst_gw[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         int src_gw[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

         int bshift[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

         for (int nd = 0; nd < dim.getValue(); nd++) {
            dst_lo[nd] = dst_box.lower(nd);
            dst_w[nd] = dst_box.numberCells(nd);
            src_glo[nd] = src_ghost_box.lower(nd);
            dst_gw[nd] = dst_ghost_box.numberCells(nd);
            src_gw[nd] = src_ghost_box.numberCells(nd);

            bshift[nd] = back_shift(nd);
         }

         const int dst_offset = dst_ghost_box.size();
         const int src_offset = src_ghost_box.size();

         int dst_bd = dst_ghost_box.offset(dst_box.lower());
         int src_bd = 0;

         int dst_ba;
         int src_ba;

         int src_in[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

         /*
          * two_id and one_id used to avoid array-bounds warnings when
          * compiling with dim < 3.  Since this is within a ((dim == tbox::Dimension(3)))
          * conditional, two_id will always be 2 and one_id will always
          * be 1. zero_id is used to keep the style consistent.
          */
         const int two_id = (dim > tbox::Dimension(2)) ? 2 : (dim.getValue() - 1);
         const int one_id = (dim > tbox::Dimension(1)) ? 1 : 0;
         const int zero_id = 0;
         for (int d = 0; d < depth; d++) {

            for (int i2 = 0; i2 < dst_w[two_id]; i2++) {
               for (int i1 = 0; i1 < dst_w[one_id]; i1++) {
                  for (int i0 = 0; i0 < dst_w[zero_id]; i0++) {

                     dst_ba = (i0
                               + i1 * dst_gw[zero_id]
                               + i2 * dst_gw[one_id] * dst_gw[zero_id]
                               ) + dst_bd;

                     src_in[zero_id] = dst_lo[zero_id] + i0;
                     src_in[one_id] = dst_lo[one_id] + i1;
                     src_in[two_id] = dst_lo[two_id] + i2;

                     hier::Transformation::rotateIndex(
                        src_in, dim, back_rotate);
                     for (int sb = 0; sb < dim.getValue(); sb++) {
                        src_in[sb] += bshift[sb];
                     }

                     src_ba = ((src_in[zero_id] - src_glo[zero_id])
                               + (src_in[one_id] - src_glo[one_id])
                               * src_gw[zero_id]
                               + (src_in[two_id] - src_glo[two_id])
                               * src_gw[one_id] * src_gw[zero_id]
                               ) + src_bd;

                     dst_ptr[dst_ba] = src_ptr[src_ba];
                  }
               }
            }

            dst_bd += dst_offset;
            src_bd += src_offset;

         }
      }
   } else {
      TBOX_ERROR(
         "MultiblockCellDataTranslator<TYPE>::translateAndCopyData : dim = 1 or > 3 not implemented");
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
void MultiblockCellDataTranslator<TYPE>::translateAndCopyArrayData(
   ArrayData<TYPE>& dst,
   const ArrayData<TYPE>& src,
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
   }
}

}
}
#endif
