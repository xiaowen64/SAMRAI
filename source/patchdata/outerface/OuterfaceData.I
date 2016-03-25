//
// File:	OuterfaceData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated outerface centered patch data type
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
int OuterfaceData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
TYPE *OuterfaceData<DIM,TYPE>::getPointer(const int axis,
                                            const int face,
                                            const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][face].getPointer(d));
}

template<int DIM, class TYPE> inline
const TYPE *OuterfaceData<DIM,TYPE>::getPointer(const int axis,
                                                  const int face,
                                                  const int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][face].getPointer(d));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>
&OuterfaceData<DIM,TYPE>::getArrayData(const int axis,
					 const int face)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
#endif
   return(d_data[axis][face]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>
&OuterfaceData<DIM,TYPE>::getArrayData(const int axis,
					 const int face) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
#endif
   return(d_data[axis][face]);
}

template<int DIM, class TYPE> inline
TYPE& OuterfaceData<DIM,TYPE>::operator()(
   const FaceIndex<DIM>& i, const int face, const int d)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][face](i,d));
}

template<int DIM, class TYPE> inline
const TYPE& OuterfaceData<DIM,TYPE>::operator()(
   const FaceIndex<DIM>& i, const int face, const int d) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert((face == 0) || (face == 1));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis][face](i,d));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
