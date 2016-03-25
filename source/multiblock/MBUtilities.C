//
// File:	MBUtilities.C
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	utility functions for multiblock
//

#ifndef included_mblk_MBUtilities_C
#define included_mblk_MBUtilities_C

#include "MBUtilities.h"

#include "MBDataUtilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif


namespace SAMRAI {
    namespace mblk {

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
   hier::PatchData<DIM>& dst,
   const hier::PatchData<DIM>& src,
   const hier::IntVector<DIM>& shift,
   const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   const pdat::CellData<DIM,double> *t_cell_dbl_src =
      dynamic_cast<const pdat::CellData<DIM,double> *>(&src);
   const pdat::FaceData<DIM,double> *t_face_dbl_src =
      dynamic_cast<const pdat::FaceData<DIM,double> *>(&src);
   const pdat::EdgeData<DIM,double> *t_edge_dbl_src =
      dynamic_cast<const pdat::EdgeData<DIM,double> *>(&src);
   const pdat::NodeData<DIM,double> *t_node_dbl_src =
      dynamic_cast<const pdat::NodeData<DIM,double> *>(&src);
   const pdat::SideData<DIM,double> *t_side_dbl_src =
      dynamic_cast<const pdat::SideData<DIM,double> *>(&src);
   const pdat::CellData<DIM,int> *t_cell_int_src =
      dynamic_cast<const pdat::CellData<DIM,int> *>(&src);
   const pdat::FaceData<DIM,int> *t_face_int_src =
      dynamic_cast<const pdat::FaceData<DIM,int> *>(&src);
   const pdat::EdgeData<DIM,int> *t_edge_int_src =
      dynamic_cast<const pdat::EdgeData<DIM,int> *>(&src);
   const pdat::NodeData<DIM,int> *t_node_int_src =
      dynamic_cast<const pdat::NodeData<DIM,int> *>(&src);
   const pdat::SideData<DIM,int> *t_side_int_src =
      dynamic_cast<const pdat::SideData<DIM,int> *>(&src);

#ifdef HAVE_BOOL
   const pdat::CellData<DIM,bool> *t_cell_bool_src =
      dynamic_cast<const pdat::CellData<DIM,bool> *>(&src);
   const pdat::FaceData<DIM,bool> *t_face_bool_src =
      dynamic_cast<const pdat::FaceData<DIM,bool> *>(&src);
   const pdat::EdgeData<DIM,bool> *t_edge_bool_src =
      dynamic_cast<const pdat::EdgeData<DIM,bool> *>(&src);
   const pdat::NodeData<DIM,bool> *t_node_bool_src =
      dynamic_cast<const pdat::NodeData<DIM,bool> *>(&src);
   const pdat::SideData<DIM,bool> *t_side_bool_src =
      dynamic_cast<const pdat::SideData<DIM,bool> *>(&src);
#endif

#ifdef HAVE_DCOMPLEX
   const pdat::CellData<DIM,dcomplex> *t_cell_dcmp_src =
      dynamic_cast<const pdat::CellData<DIM,dcomplex> *>(&src);
   const pdat::FaceData<DIM,dcomplex> *t_face_dcmp_src =
      dynamic_cast<const pdat::FaceData<DIM,dcomplex> *>(&src);
   const pdat::EdgeData<DIM,dcomplex> *t_edge_dcmp_src =
      dynamic_cast<const pdat::EdgeData<DIM,dcomplex> *>(&src);
   const pdat::NodeData<DIM,dcomplex> *t_node_dcmp_src =
      dynamic_cast<const pdat::NodeData<DIM,dcomplex> *>(&src);
   const pdat::SideData<DIM,dcomplex> *t_side_dcmp_src =
      dynamic_cast<const pdat::SideData<DIM,dcomplex> *>(&src);
#endif

#ifdef HAVE_FLOAT
   const pdat::CellData<DIM,float> *t_cell_flt_src =
      dynamic_cast<const pdat::CellData<DIM,float> *>(&src);
   const pdat::FaceData<DIM,float> *t_face_flt_src =
      dynamic_cast<const pdat::FaceData<DIM,float> *>(&src);
   const pdat::EdgeData<DIM,float> *t_edge_flt_src =
      dynamic_cast<const pdat::EdgeData<DIM,float> *>(&src);
   const pdat::NodeData<DIM,float> *t_node_flt_src =
      dynamic_cast<const pdat::NodeData<DIM,float> *>(&src);
   const pdat::SideData<DIM,float> *t_side_flt_src =
      dynamic_cast<const pdat::SideData<DIM,float> *>(&src);
#endif

   if (t_cell_dbl_src != NULL) {
      MBDataUtilities<DIM,double>::translateAndCopyCellData(
         (pdat::CellData<DIM,double>&) dst,
         (pdat::CellData<DIM,double>&) src,
         shift,
         rotate);
   } else if (t_node_dbl_src != NULL) {
      MBDataUtilities<DIM,double>::translateAndCopyNodeData(
         (pdat::NodeData<DIM,double>&) dst,
         (pdat::NodeData<DIM,double>&) src,
         shift,
         rotate);
   } else if (t_face_dbl_src != NULL) {
      MBDataUtilities<DIM,double>::translateAndCopyFaceData(
         (pdat::FaceData<DIM,double>&) dst,
         (pdat::FaceData<DIM,double>&) src,
         shift,
         rotate);
   } else if (t_side_dbl_src != NULL) {
      MBDataUtilities<DIM,double>::translateAndCopySideData(
         (pdat::SideData<DIM,double>&) dst,
         (pdat::SideData<DIM,double>&) src,
         shift,
         rotate);
   } else if (t_edge_dbl_src != NULL) {
      MBDataUtilities<DIM,double>::translateAndCopyEdgeData(
         (pdat::EdgeData<DIM,double>&) dst,
         (pdat::EdgeData<DIM,double>&) src,
         shift,
         rotate);
   } else if (t_cell_int_src != NULL) {
      MBDataUtilities<DIM,int>::translateAndCopyCellData(
         (pdat::CellData<DIM,int>&) dst,
         (pdat::CellData<DIM,int>&) src,
         shift,
         rotate);
   } else if (t_node_int_src != NULL) {
      MBDataUtilities<DIM,int>::translateAndCopyNodeData(
         (pdat::NodeData<DIM,int>&) dst,
         (pdat::NodeData<DIM,int>&) src,
         shift,
         rotate);
   } else if (t_face_int_src != NULL) {
      MBDataUtilities<DIM,int>::translateAndCopyFaceData(
         (pdat::FaceData<DIM,int>&) dst,
         (pdat::FaceData<DIM,int>&) src,
         shift,
         rotate);
   } else if (t_side_int_src != NULL) {
      MBDataUtilities<DIM,int>::translateAndCopySideData(
         (pdat::SideData<DIM,int>&) dst,
         (pdat::SideData<DIM,int>&) src,
         shift,
         rotate);
   } else if (t_edge_int_src != NULL) {
      MBDataUtilities<DIM,int>::translateAndCopyEdgeData(
         (pdat::EdgeData<DIM,int>&) dst,
         (pdat::EdgeData<DIM,int>&) src,
         shift,
         rotate);
   }
#ifdef HAVE_BOOL
   else if (t_node_bool_src != NULL) {
      MBDataUtilities<DIM,bool>::translateAndCopyNodeData(
         (pdat::NodeData<DIM,bool>&) dst,
         (pdat::NodeData<DIM,bool>&) src,
         shift,
         rotate);
   } else if (t_face_bool_src != NULL) {
      MBDataUtilities<DIM,bool>::translateAndCopyFaceData(
         (pdat::FaceData<DIM,bool>&) dst,
         (pdat::FaceData<DIM,bool>&) src,
         shift,
         rotate);
   } else if (t_side_bool_src != NULL) {
      MBDataUtilities<DIM,bool>::translateAndCopySideData(
         (pdat::SideData<DIM,bool>&) dst,
         (pdat::SideData<DIM,bool>&) src,
         shift,
         rotate);
   } else if (t_edge_bool_src != NULL) {
      MBDataUtilities<DIM,bool>::translateAndCopyEdgeData(
         (pdat::EdgeData<DIM,bool>&) dst,
         (pdat::EdgeData<DIM,bool>&) src,
         shift,
         rotate);
   } else if (t_cell_bool_src != NULL) {
      MBDataUtilities<DIM,bool>::translateAndCopyCellData(
         (pdat::CellData<DIM,bool>&) dst,
         (pdat::CellData<DIM,bool>&) src,
         shift,
         rotate);
   }
#endif
#ifdef HAVE_DCOMPLEX
else if (t_cell_dcmp_src != NULL) {
      MBDataUtilities<DIM,dcomplex>::translateAndCopyCellData(
         (pdat::CellData<DIM,dcomplex>&) dst,
         (pdat::CellData<DIM,dcomplex>&) src,
         shift,
         rotate);
   } else if (t_node_dcmp_src != NULL) {
      MBDataUtilities<DIM,dcomplex>::translateAndCopyNodeData(
         (pdat::NodeData<DIM,dcomplex>&) dst,
         (pdat::NodeData<DIM,dcomplex>&) src,
         shift,
         rotate);
   } else if (t_face_dcmp_src != NULL) {
      MBDataUtilities<DIM,dcomplex>::translateAndCopyFaceData(
         (pdat::FaceData<DIM,dcomplex>&) dst,
         (pdat::FaceData<DIM,dcomplex>&) src,
         shift,
         rotate);
   } else if (t_side_dcmp_src != NULL) {
      MBDataUtilities<DIM,dcomplex>::translateAndCopySideData(
         (pdat::SideData<DIM,dcomplex>&) dst,
         (pdat::SideData<DIM,dcomplex>&) src,
         shift,
         rotate);
   } else if (t_edge_dcmp_src != NULL) {
      MBDataUtilities<DIM,dcomplex>::translateAndCopyEdgeData(
         (pdat::EdgeData<DIM,dcomplex>&) dst,
         (pdat::EdgeData<DIM,dcomplex>&) src,
         shift,
         rotate);
   }
#endif
#ifdef HAVE_FLOAT
   else if (t_cell_flt_src != NULL) {
      MBDataUtilities<DIM,float>::translateAndCopyCellData(
         (pdat::CellData<DIM,float>&) dst,
         (pdat::CellData<DIM,float>&) src,
         shift,
         rotate);
   } else if (t_node_flt_src != NULL) {
      MBDataUtilities<DIM,float>::translateAndCopyNodeData(
         (pdat::NodeData<DIM,float>&) dst,
         (pdat::NodeData<DIM,float>&) src,
         shift,
         rotate);
   } else if (t_face_flt_src != NULL) {
      MBDataUtilities<DIM,float>::translateAndCopyFaceData(
         (pdat::FaceData<DIM,float>&) dst,
         (pdat::FaceData<DIM,float>&) src,
         shift,
         rotate);
   } else if (t_side_flt_src != NULL) {
      MBDataUtilities<DIM,float>::translateAndCopySideData(
         (pdat::SideData<DIM,float>&) dst,
         (pdat::SideData<DIM,float>&) src,
         shift,
         rotate);
   } else if (t_edge_flt_src != NULL) {
      MBDataUtilities<DIM,float>::translateAndCopyEdgeData(
         (pdat::EdgeData<DIM,float>&) dst,
         (pdat::EdgeData<DIM,float>&) src,
         shift,
         rotate);
   }
#endif
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
   }

   if (DIM == 3) {
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
      assert(axis < DIM);
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
