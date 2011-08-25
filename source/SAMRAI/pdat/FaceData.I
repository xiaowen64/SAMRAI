/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated face centered patch data type
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int FaceData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * FaceData<TYPE>::getPointer(
   int face_normal,
   int depth)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * FaceData<TYPE>::getPointer(
   int face_normal,
   int depth) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& FaceData<TYPE>::operator () (
   const FaceIndex& i,
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
const TYPE&
FaceData<TYPE>::operator () (
   const FaceIndex& i,
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
FaceData<TYPE>::getArrayData(
   int face_normal)
{

   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return d_data[face_normal];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>&
FaceData<TYPE>::getArrayData(
   int face_normal) const
{

   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return d_data[face_normal];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void FaceData<TYPE>::copyOnBox(
   const FaceData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const hier::Box face_box = FaceGeometry::toFaceBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), face_box);
   }

}

}
}
