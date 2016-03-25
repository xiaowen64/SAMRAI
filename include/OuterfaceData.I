//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outerface/OuterfaceData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated outerface centered patch data type
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int OuterfaceData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
TYPE* OuterfaceData<DIM,TYPE>::getPointer(int face_normal,
                                          int side,
                                          int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   return(d_data[face_normal][side].getPointer(d));
}

template<int DIM, class TYPE> inline
const TYPE* OuterfaceData<DIM,TYPE>::getPointer(int face_normal,
                                                int side,
                                                int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   return(d_data[face_normal][side].getPointer(d));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>&
OuterfaceData<DIM,TYPE>::getArrayData(int face_normal,
                                      int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[face_normal][side]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>&
OuterfaceData<DIM,TYPE>::getArrayData(int face_normal,
                                      int side) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((face_normal >= 0) && (face_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[face_normal][side]);
}

template<int DIM, class TYPE> inline
TYPE& OuterfaceData<DIM,TYPE>::operator()(
   const FaceIndex<DIM>& i, int side, int depth)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis][side](i,depth));
}

template<int DIM, class TYPE> inline
const TYPE& OuterfaceData<DIM,TYPE>::operator()(
   const FaceIndex<DIM>& i, int side, int depth) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis][side](i,depth));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
