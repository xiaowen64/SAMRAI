/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Templated outeredge centered patch data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int OuteredgeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool OuteredgeData<TYPE>::dataExists(
   int axis,
   int face_normal) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return d_data[axis][face_normal][0].isInitialized();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * OuteredgeData<TYPE>::getPointer(
   int axis,
   int face_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][face_normal][side].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * OuteredgeData<TYPE>::getPointer(
   int axis,
   int face_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][face_normal][side].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>&
OuteredgeData<TYPE>::getArrayData(
   int axis,
   int face_normal,
   int side)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[axis][face_normal][side];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>&
OuteredgeData<TYPE>::getArrayData(
   int axis,
   int face_normal,
   int side) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[axis][face_normal][side];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & OuteredgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

      if (face_normal != axis) {

         for (int side = 0; side < 2; ++side) {

            if (d_data[axis][face_normal][side].getBox().contains(i)) {
               return d_data[axis][face_normal][side](i, depth);
            }

         }  // iterate over lower/upper sides

      }  // data is undefined when axis == face_normal

   }  // iterate over face normal directions

   TBOX_ERROR("Attempt to access OuteredgeData value with bad index"
      " edge index " << i << " with axis = " << axis << std::endl);
   return d_data[0][0][0](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &OuteredgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

      if (face_normal != axis) {

         for (int side = 0; side < 2; ++side) {

            if (d_data[axis][face_normal][side].getBox().contains(i)) {
               return d_data[axis][face_normal][side](i, depth);
            }

         }  // iterate over lower/upper sides

      }  // data is undefined when axis == face_normal

   }  // iterate over face normal directions

   TBOX_ERROR("Attempt to access OuteredgeData value with bad index"
      " edge index " << i << " with axis = " << axis << std::endl);
   return d_data[0][0][0](i, depth);
}

}
}
