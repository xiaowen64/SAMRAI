/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Templated side centered patch data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const hier::IntVector& SideData<TYPE>::getDirectionVector() const
{
   return d_directions;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int SideData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * SideData<TYPE>::getPointer(
   int side_normal,
   int depth)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * SideData<TYPE>::getPointer(
   int side_normal,
   int depth) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & SideData<TYPE>::operator () (
   const SideIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT(d_directions(axis));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &
SideData<TYPE>::operator () (
   const SideIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT(d_directions(axis));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>&
SideData<TYPE>::getArrayData(
   int side_normal)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));

   return d_data[side_normal];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>&
SideData<TYPE>::getArrayData(
   int side_normal) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));

   return d_data[side_normal];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void SideData<TYPE>::copyOnBox(
   const SideData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const hier::Box side_box = SideGeometry::toSideBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), side_box);
   }

}

}
}
