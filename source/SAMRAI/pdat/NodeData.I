/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated node centered patch data type
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int NodeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>& NodeData<TYPE>::getArrayData()
{
   return d_data;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>& NodeData<TYPE>::getArrayData() const
{
   return d_data;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * NodeData<TYPE>::getPointer(
   int depth)
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data.getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * NodeData<TYPE>::getPointer(
   int depth) const
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data.getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& NodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data(i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE&
NodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data(i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void NodeData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data.fill(t, d);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void NodeData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data.fill(t, NodeGeometry::toNodeBox(box), d);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void NodeData<TYPE>::fillAll(
   const TYPE& t)
{
   d_data.fillAll(t);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void NodeData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   d_data.fillAll(t, NodeGeometry::toNodeBox(box));
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void NodeData<TYPE>::copyOnBox(
   const NodeData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);
   const hier::Box node_box = NodeGeometry::toNodeBox(box);
   d_data.copy(src.getArrayData(), node_box);
}

}
}
