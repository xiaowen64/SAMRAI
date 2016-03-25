//
// File:	SideData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated side centered patch data type
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
const hier::IntVector<DIM>& SideData<DIM,TYPE>::getDirectionVector() const
{
   return(d_directions);
}

template<int DIM, class TYPE> inline
int SideData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
TYPE *SideData<DIM,TYPE>::getPointer(const int axis, const int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert(d_directions(axis));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis].getPointer(d));
}

template<int DIM, class TYPE> inline
const TYPE *SideData<DIM,TYPE>::getPointer(const int axis, const int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert(d_directions(axis));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis].getPointer(d));
}

template<int DIM, class TYPE> inline
TYPE& SideData<DIM,TYPE>::operator()(const SideIndex<DIM>& i, const int d)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert(d_directions(axis));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis](i,d));
}

template<int DIM, class TYPE> inline
const TYPE&
SideData<DIM,TYPE>::operator()(const SideIndex<DIM>& i, const int d) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert(d_directions(axis));
   assert((d >= 0) && (d < d_depth));
#endif
   return(d_data[axis](i,d));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>&
SideData<DIM,TYPE>::getArrayData(const int axis)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert(d_directions(axis));
#endif
   return(d_data[axis]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>&
SideData<DIM,TYPE>::getArrayData(const int axis) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((axis >= 0) && (axis < DIM));
   assert(d_directions(axis));
#endif
   return(d_data[axis]);
}

template<int DIM, class TYPE> inline 
void SideData<DIM,TYPE>::copyOnBox(
   const SideData<DIM,TYPE>& src,
   const hier::Box<DIM>& box)
{
   for (int axis = 0; axis < DIM; axis++) {
      const hier::Box<DIM> side_box = SideGeometry<DIM>::toSideBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), side_box);
   }

}
#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
