//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/side/SideData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated side centered patch data type
//


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
TYPE* SideData<DIM,TYPE>::getPointer(int side_normal, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT(d_directions(side_normal));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[side_normal].getPointer(depth));
}

template<int DIM, class TYPE> inline
const TYPE* SideData<DIM,TYPE>::getPointer(int side_normal, int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT(d_directions(side_normal));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[side_normal].getPointer(depth));
}

template<int DIM, class TYPE> inline
TYPE& SideData<DIM,TYPE>::operator()(const SideIndex<DIM>& i, int depth)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT(d_directions(axis));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis](i,depth));
}

template<int DIM, class TYPE> inline
const TYPE&
SideData<DIM,TYPE>::operator()(const SideIndex<DIM>& i, int depth) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT(d_directions(axis));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis](i,depth));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>&
SideData<DIM,TYPE>::getArrayData(int side_normal)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT(d_directions(side_normal));
#endif
   return(d_data[side_normal]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>&
SideData<DIM,TYPE>::getArrayData(int side_normal) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((side_normal >= 0) && (side_normal < DIM));
   TBOX_ASSERT(d_directions(side_normal));
#endif
   return(d_data[side_normal]);
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
