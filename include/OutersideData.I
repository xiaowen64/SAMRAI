//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outerside/OutersideData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated outerside centered patch data type
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int OutersideData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
TYPE* OutersideData<DIM,TYPE>::getPointer(int side_normal,
                                          int side,
                                          int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[side_normal][side].getPointer(depth));
}

template<int DIM, class TYPE> inline
const TYPE* OutersideData<DIM,TYPE>::getPointer(int side_normal,
                                                int side,
                                                int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[side_normal][side].getPointer(depth));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>& 
OutersideData<DIM,TYPE>::getArrayData(int side_normal,
                                      int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[side_normal][side]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>& 
OutersideData<DIM,TYPE>::getArrayData(int side_normal,
                                      int side) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT((side == 0) || (side == 1));
#endif
   return(d_data[side_normal][side]);
}

template<int DIM, class TYPE> inline
TYPE& OutersideData<DIM,TYPE>::operator()(
   const SideIndex<DIM>& i, int side, int depth)
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
const TYPE& OutersideData<DIM,TYPE>::operator()(
   const SideIndex<DIM>& i, int side, int depth) const
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
