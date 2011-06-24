/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Templated edge centered patch data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int EdgeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * EdgeData<TYPE>::getPointer(
   int axis,
   int depth)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * EdgeData<TYPE>::getPointer(
   int axis,
   int depth) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & EdgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &
EdgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>&
EdgeData<TYPE>::getArrayData(
   int axis)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   return d_data[axis];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>&
EdgeData<TYPE>::getArrayData(
   int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   return d_data[axis];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void EdgeData<TYPE>::copyOnBox(
   const EdgeData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const hier::Box edge_box = EdgeGeometry::toEdgeBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), edge_box);
   }

}

}
}
