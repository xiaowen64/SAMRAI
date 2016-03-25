//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/multiblock/MBDataUtilities.C $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated operations for copying patch data.
//

#ifndef included_pdat_MBDataUtilities_C
#define included_pdat_MBDataUtilities_C

#include "MBDataUtilities.h"

#include "MBUtilities.h"


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
MBDataUtilities<DIM,TYPE>::MBDataUtilities()
{
}

template <int DIM, class TYPE>
MBDataUtilities<DIM,TYPE>::~MBDataUtilities()
{
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for cell data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MBDataUtilities<DIM,TYPE>::translateAndCopyCellData(
   pdat::CellData<DIM,TYPE>& dst,
   const pdat::CellData<DIM,TYPE>& src,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   if (DIM == 2) {
   translateAndCopyArrayData(dst.getArrayData(),
                             src.getArrayData(),
                             shift,
                             rotate);
   } else if (DIM == 3) {
      if (rotate == 0) {
         translateAndCopyArrayData(dst.getArrayData(),
                                   src.getArrayData(),
                                   shift,
                                   rotate);
      } else {

         const typename
         hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier back_rotate =
            hier::MultiblockPatchHierarchy<DIM>::getReverseRotationIdentifier(rotate);
   
         hier::IntVector<DIM> back_shift;
         hier::MultiblockPatchHierarchy<DIM>::calculateReverseShift(
            back_shift, shift, back_rotate);

         TYPE *const dst_ptr = dst.getArrayData().getPointer();
         const TYPE *const src_ptr = src.getArrayData().getPointer();

         const int depth = ((dst.getDepth() < src.getDepth()) ? 
                             dst.getDepth() : src.getDepth());

         const hier::Box<DIM> dst_box = dst.getBox();
         const hier::Box<DIM> src_box = src.getBox();
         const hier::Box<DIM> dst_ghost_box = dst.getGhostBox();
         const hier::Box<DIM> src_ghost_box = src.getGhostBox();

         int dst_lo[DIM];
         int dst_w[DIM];
         int src_glo[DIM];
         int dst_gw[DIM];
         int src_gw[DIM];

         int bshift[DIM];

         for (int nd = 0; nd < DIM; nd++) {
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

         int src_in[DIM];

         /*
          * two_id and one_id used to avoid array-bounds warnings when
          * compiling with DIM < 3.  Since this is within a (DIM == 3)
          * conditional, two_id will always be 2 and one_id will always
          * be 1. zero_id is used to keep the style consistent.
          */
         const int two_id = (DIM > 2) ? 2 : (DIM-1);
         const int one_id = (DIM > 1) ? 1 : 0;
         const int zero_id = 0;
         for (int d = 0; d < depth; d++) {

            for (int i2 = 0; i2 < dst_w[two_id]; i2++) {
               for (int i1 = 0; i1 < dst_w[one_id]; i1++) {
                  for (int i0 = 0; i0 < dst_w[zero_id]; i0++) {

                     dst_ba = (i0 
                             + i1*dst_gw[zero_id]
                             + i2*dst_gw[one_id]*dst_gw[zero_id] 
                              ) + dst_bd;

                     src_in[zero_id] = dst_lo[zero_id] + i0;
                     src_in[one_id] = dst_lo[one_id] + i1;
                     src_in[two_id] = dst_lo[two_id] + i2;

                     hier::MBUtilities<DIM>::rotateIndex(src_in, back_rotate);
                     for (int sb = 0; sb < DIM; sb++) {
                        src_in[sb] += bshift[sb];
                     }

                     src_ba = ((src_in[zero_id] - src_glo[zero_id])
                              + (src_in[one_id] - src_glo[one_id])*
                                src_gw[zero_id]
                              + (src_in[two_id] - src_glo[two_id])*
                                src_gw[one_id]*src_gw[zero_id]
                              ) + src_bd;
 
                     dst_ptr[dst_ba] = src_ptr[src_ba];
                  }
               }
            }

            dst_bd += dst_offset;
            src_bd += src_offset;

         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for node data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MBDataUtilities<DIM,TYPE>::translateAndCopyNodeData(
   pdat::NodeData<DIM,TYPE>& dst,
   const pdat::NodeData<DIM,TYPE>& src,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{

   bool no_rotate;
   if (rotate != 0) {
      no_rotate = false;
   } else {
      no_rotate = true;
   }
     
   if (no_rotate) {
      dst.getArrayData().copy(src.getArrayData(), 
                              dst.getArrayData().getBox(), shift);
   } else if (DIM == 2) {
      hier::Box<DIM> rotatebox(src.getArrayData().getBox());
      int num_rotations = rotate;

      rotatebox.rotate(num_rotations);

      const hier::Box<DIM> copybox = dst.getArrayData().getBox() *
                                (pdat::NodeGeometry<DIM>::toNodeBox(
                                   hier::Box<DIM>::shift(rotatebox, shift)));

      if (!copybox.empty()) {
         TYPE *const dst_ptr = dst.getArrayData().getPointer();
         const TYPE *const src_ptr = src.getArrayData().getPointer();

         const int depth = (dst.getArrayData().getDepth() < 
                            src.getArrayData().getDepth() ?
                            dst.getArrayData().getDepth() : 
                            src.getArrayData().getDepth());

         const int box_w0 = copybox.numberCells(0);

         const int dst_w0 = dst.getArrayData().getBox().numberCells(0);
         int src_w0 = src.getArrayData().getBox().numberCells(0);
         if (num_rotations == 3) {
            src_w0 = -src_w0;
         }
         const int box_w1 = copybox.numberCells(1);

         const int dst_offset = dst.getArrayData().getOffset();
         const int src_offset = src.getArrayData().getOffset();

         int dst_bd = dst.getArrayData().getBox().offset(copybox.lower());
         hier::Index<DIM> src_index(copybox.lower()-shift);

         // rotate src_index 4-num_rotations;
         for (int r = 0; r < 4-num_rotations; r++) {
            hier::Index<DIM> tmp_index(src_index);
            src_index(0) = tmp_index(1);
            src_index(1) = -tmp_index(0);
         }

         int src_bd_orig = src.getArrayData().getBox().offset(src_index);

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
      for (pdat::NodeIterator<DIM> ni(dst.getBox()); ni; ni++) {
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

         for (int d = 0; d < dst.getDepth(); d++) {
            dst(dst_index, d) = src(src_index, d);
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for face data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MBDataUtilities<DIM,TYPE>::translateAndCopyFaceData(
   pdat::FaceData<DIM,TYPE>& dst,
   const pdat::FaceData<DIM,TYPE>& src,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{

   if (rotate == 0) {
      for (int axis = 0; axis < DIM; axis++) {
         hier::Index<DIM> face_shift;
         for (int d = 0; d < DIM; d++) {
            face_shift(d) = shift((axis+d)%DIM);
         }
         translateAndCopyArrayData(dst.getArrayData(axis),
                                   src.getArrayData(axis),
                                   face_shift,
                                   rotate);
      }
   } else if (DIM == 2) {
      for (int axis = 0; axis < DIM; axis++) {
         for (pdat::FaceIterator<DIM> fi(dst.getBox(),axis); fi; fi++) {
            pdat::FaceIndex<DIM> dst_index(fi());
            hier::Index<DIM> dst_xyz_index(dst_index);
            if (axis == 1) {
               dst_xyz_index(0) = dst_index(1);
               dst_xyz_index(1) = dst_index(0);
            }

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

            pdat::FaceIndex<DIM> src_index;
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

            for (int d = 0; d < dst.getDepth(); d++) {
               dst(dst_index, d) = src(src_index, d);
            }
         }
      }
   } else if (DIM == 3) {
      for (int axis = 0; axis < DIM; axis++) {
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
         for (pdat::FaceIterator<DIM> fi(dst.getBox(),axis); fi; fi++) {
            pdat::FaceIndex<DIM> dst_index(fi());
            hier::Index<DIM> dst_xyz_index(dst_index);
            if (axis == 1) {
               dst_xyz_index(0) = dst_index(2);
               dst_xyz_index(1) = dst_index(0);
               dst_xyz_index(2) = dst_index(1);
            } else if (axis == 2) {
               dst_xyz_index(0) = dst_index(1);
               dst_xyz_index(1) = dst_index(2);
               dst_xyz_index(2) = dst_index(0);
            }

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

            pdat::FaceIndex<DIM> src_index;
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
                     src_index(0)++;
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
                     src_index(0)++;
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
                     src_index(0)++;
                     break;

                  default:
                     TBOX_ERROR(" ");
                     break;
               }
            }

            src_index.setAxis(src_axis);

            for (int d = 0; d < dst.getDepth(); d++) {
               dst(dst_index, d) = src(src_index, d);
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for side data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MBDataUtilities<DIM,TYPE>::translateAndCopySideData(
   pdat::SideData<DIM,TYPE>& dst,
   const pdat::SideData<DIM,TYPE>& src,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dst.getDirectionVector() == src.getDirectionVector());
#endif

   hier::IntVector<DIM> dir_vector(dst.getDirectionVector());
   if (rotate == 0) {
      for (int axis = 0; axis < DIM; axis++) {
         if (dir_vector(axis)) {
            translateAndCopyArrayData(dst.getArrayData(axis),
                                      src.getArrayData(axis),
                                      shift,
                                      rotate);
         }
      }
   } else if (DIM == 2) {
      for (int axis = 0; axis < DIM; axis++) {
         if (dir_vector(axis)) {
            for (pdat::SideIterator<DIM> fi(dst.getBox(),axis); fi; fi++) {
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
   
               for (int d = 0; d < dst.getDepth(); d++) {
                  dst(dst_index, d) = src(src_index, d);
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

            for (pdat::SideIterator<DIM> si(dst.getBox(),axis); si; si++) {
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

               for (int d = 0; d < dst.getDepth(); d++) {
                  dst(dst_index, d) = src(src_index, d);
               }
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Translation and copy for edge data                                    *
*                                                                       *
*************************************************************************
*/

template <int DIM, class TYPE>
void MBDataUtilities<DIM,TYPE>::translateAndCopyEdgeData(
   pdat::EdgeData<DIM,TYPE>& dst,
   const pdat::EdgeData<DIM,TYPE>& src,
   const hier::IntVector<DIM>& shift,
   const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   if (rotate == 0) {
      for (int axis = 0; axis < DIM; axis++) {
         translateAndCopyArrayData(dst.getArrayData(axis),
                                   src.getArrayData(axis),
                                   shift,
                                   rotate);
      }
   } else if (DIM == 2) {
      for (int axis = 0; axis < DIM; axis++) {
         for (pdat::EdgeIterator<DIM> fi(dst.getBox(),axis); fi; fi++) {

            pdat::EdgeIndex<DIM> dst_index(fi());
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

            pdat::EdgeIndex<DIM> src_index;
            if (src_axis == 0) {
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

            for (int d = 0; d < dst.getDepth(); d++) {
               dst(dst_index, d) = src(src_index, d);
            }
         }
      }
   } else if (DIM == 3) {
      for (int axis = 0; axis < DIM; axis++) {
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

         for (pdat::EdgeIterator<DIM> fi(dst.getBox(),axis); fi; fi++) {
            pdat::EdgeIndex<DIM> dst_index(fi());

            typename
            hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier back_rotate =
               hier::MultiblockPatchHierarchy<DIM>::
               getReverseRotationIdentifier(rotate);

            hier::Box<DIM> src_box(dst_index, dst_index);

            src_box.rotate(back_rotate);

            hier::IntVector<DIM> back_shift;
            hier::MultiblockPatchHierarchy<DIM>::calculateReverseShift(
               back_shift, shift, back_rotate);

            src_box.shift(back_shift);

            pdat::EdgeIndex<DIM> src_index;
            for (int i = 0; i < DIM; i++) { 
               src_index(i) = src_box.lower()(i);
            }

            switch (rotate) {

               case hier::MultiblockPatchHierarchy<DIM>::IUP_JUP_KUP:
               case hier::MultiblockPatchHierarchy<DIM>::KUP_IUP_JUP:
               case hier::MultiblockPatchHierarchy<DIM>::JUP_KUP_IUP:
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KUP_JUP:
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KUP_JUP_IDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(0)++;
                  }
                  break;
          
               case hier::MultiblockPatchHierarchy<DIM>::JUP_IDOWN_KUP:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JUP_IUP:
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IUP_KDOWN_JUP:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::JUP_IUP_KDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IDOWN_JUP:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JUP_KDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(2)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::JUP_KDOWN_IDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IUP_KUP:
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IUP_KUP_JDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(1)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KUP_JDOWN_IUP:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KUP_IDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IDOWN_JDOWN_KUP:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KUP_IDOWN_JDOWN:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(1)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::JDOWN_KDOWN_IUP:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(1)++; 
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KDOWN_IUP_JDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(1)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IUP_JDOWN_KDOWN:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::JDOWN_IDOWN_KDOWN:
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::KDOWN_JDOWN_IDOWN:
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(1)++;
                  }
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               case hier::MultiblockPatchHierarchy<DIM>::IDOWN_KDOWN_JDOWN:
                  if ( (axis == 1) || (axis == 2) ) {
                     src_index(0)++;
                  }
                  if ( (axis == 0) || (axis == 1) ) {
                     src_index(1)++;
                  }
                  if ( (axis == 0) || (axis == 2) ) {
                     src_index(2)++;
                  }
                  break;

               default:
                  break;

            }

            src_index.setAxis(src_axis);

            //back rotate src_index into src index space
            //back shift src_index to needed location
            //copy data from src_index to dst_index

            for (int d = 0; d < dst.getDepth(); d++) {
               dst(dst_index, d) = src(src_index, d);
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
void MBDataUtilities<DIM,TYPE>::translateAndCopyArrayData(
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
   }
}

}
}
#endif
