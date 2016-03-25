//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/node/NodeData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Templated node centered patch data type
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
int NodeData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
ArrayData<DIM,TYPE>& NodeData<DIM,TYPE>::getArrayData()
{
   return(d_data);
}

template<int DIM, class TYPE> inline
const ArrayData<DIM,TYPE>& NodeData<DIM,TYPE>::getArrayData() const
{
   return(d_data);
}

template<int DIM, class TYPE> inline
TYPE* NodeData<DIM,TYPE>::getPointer(int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data.getPointer(depth));
}

template<int DIM, class TYPE> inline
const TYPE* NodeData<DIM,TYPE>::getPointer(int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data.getPointer(depth));
}

template<int DIM, class TYPE> inline
TYPE& NodeData<DIM,TYPE>::operator()(const NodeIndex<DIM>& i, int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data(i,depth));
}

template<int DIM, class TYPE> inline
const TYPE&
NodeData<DIM,TYPE>::operator()(const NodeIndex<DIM>& i, int depth) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
#endif
   return(d_data(i,depth));
}

template<int DIM, class TYPE> inline
void NodeData<DIM,TYPE>::fill(const TYPE& t, int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   d_data.fill(t, d);
}

template<int DIM, class TYPE> inline
void NodeData<DIM,TYPE>::fill(const TYPE& t,
                              const hier::Box<DIM>& box,
                              int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   d_data.fill(t, NodeGeometry<DIM>::toNodeBox(box), d);
}

template<int DIM, class TYPE> inline
void NodeData<DIM,TYPE>::fillAll(const TYPE& t)
{
   d_data.fillAll(t);
}

template<int DIM, class TYPE> inline
void NodeData<DIM,TYPE>::fillAll(const TYPE& t, const hier::Box<DIM>& box)
{
   d_data.fillAll(t, NodeGeometry<DIM>::toNodeBox(box));
}

template<int DIM, class TYPE> inline 
void NodeData<DIM,TYPE>::copyOnBox(
   const NodeData<DIM,TYPE>& src,
   const hier::Box<DIM>& box)
{
   const hier::Box<DIM> node_box = NodeGeometry<DIM>::toNodeBox(box);
   d_data.copy(src.getArrayData(), node_box);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
