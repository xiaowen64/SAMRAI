//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/edge/EdgeData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated edge centered patch data type
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int EdgeData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
TYPE* EdgeData<DIM,TYPE>::getPointer(int axis, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis].getPointer(depth));
}

template<int DIM, class TYPE> inline
const TYPE* EdgeData<DIM,TYPE>::getPointer(int axis, int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis].getPointer(depth));
}

template<int DIM, class TYPE> inline
TYPE& EdgeData<DIM,TYPE>::operator()(const EdgeIndex<DIM>& i, int depth)
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis](i,depth));
}

template<int DIM, class TYPE> inline
const TYPE&
EdgeData<DIM,TYPE>::operator()(const EdgeIndex<DIM>& i, int depth) const
{
   const int axis = i.getAxis();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data[axis](i,depth));
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>&
EdgeData<DIM,TYPE>::getArrayData(int axis)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
#endif
   return(d_data[axis]);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>&
EdgeData<DIM,TYPE>::getArrayData(int axis) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((axis >= 0) && (axis < DIM));
#endif
   return(d_data[axis]);
}

template<int DIM, class TYPE> inline 
void EdgeData<DIM,TYPE>::copyOnBox(
   const EdgeData<DIM,TYPE>& src,
   const hier::Box<DIM>& box)
{
   for (int axis = 0; axis < DIM; axis++) {
      const hier::Box<DIM> edge_box = EdgeGeometry<DIM>::toEdgeBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), edge_box);
   }

}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
