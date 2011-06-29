/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated outerside centered patch data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int OutersideData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * OutersideData<TYPE>::getPointer(
   int side_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal][side].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * OutersideData<TYPE>::getPointer(
   int side_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal][side].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>&
OutersideData<TYPE>::getArrayData(
   int side_normal,
   int side)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[side_normal][side];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>&
OutersideData<TYPE>::getArrayData(
   int side_normal,
   int side) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[side_normal][side];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & OutersideData<TYPE>::operator () (
   const SideIndex& i,
   int side,
   int depth)
{
   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][side](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &OutersideData<TYPE>::operator () (
   const SideIndex& i,
   int side,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][side](i, depth);
}

}
}
