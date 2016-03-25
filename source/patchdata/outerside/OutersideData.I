//
// File:	OutersideData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated outerside centered patch data type
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

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
TYPE *OutersideData<DIM,TYPE>::getPointer(const int axis,
                                            const int side,
                                            const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][side].getPointer(d));
}

template<int DIM, class TYPE> inline
const TYPE *OutersideData<DIM,TYPE>::getPointer(const int axis,
                                                  const int side,
                                                  const int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][side].getPointer(d));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>
&OutersideData<DIM,TYPE>::getArrayData(const int axis,
					 const int side)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
#endif
   return(d_data[axis][side]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>
&OutersideData<DIM,TYPE>::getArrayData(const int axis,
					 const int side) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
#endif
   return(d_data[axis][side]);
}

template<int DIM, class TYPE> inline
TYPE& OutersideData<DIM,TYPE>::operator()(
   const SideIndex<DIM>& i, const int side, const int d)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][side](i,d));
}

template<int DIM, class TYPE> inline
const TYPE& OutersideData<DIM,TYPE>::operator()(
   const SideIndex<DIM>& i, const int side, const int d) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((side == 0) || (side == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][side](i,d));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
