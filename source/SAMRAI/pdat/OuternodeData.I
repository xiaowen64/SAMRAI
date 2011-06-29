/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated outernode centered patch data type 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int OuternodeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool OuternodeData<TYPE>::dataExists(
   int face_normal) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return d_data[face_normal][0].isInitialized();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * OuternodeData<TYPE>::getPointer(
   int face_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal][side].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * OuternodeData<TYPE>::getPointer(
   int face_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal][side].getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>
& OuternodeData<TYPE>::getArrayData(
   int face_normal,
   int side)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[face_normal][side];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>
& OuternodeData<TYPE>::getArrayData(
   int face_normal,
   int side) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[face_normal][side];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & OuternodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth)
{
   for (int d = getDim().getValue() - 1; d >= 0; d--) {
      if (i[d] == d_data[d][0].getBox().lower()[d]) {
         return d_data[d][0](i, depth);
      }
      if (i[d] == d_data[d][1].getBox().upper()[d]) {
         return d_data[d][1](i, depth);
      }
   }

   /*
    * The following lines should only be executed if there's a bug
    * in the Outernode datatype.
    */
   TBOX_ERROR("Bad index used to access outernode data\n"
      << "Given index is not an outernode of this instance.\n");
   return d_data[0][0](i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &OuternodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   for (int d = getDim() - 1; d >= 0; d--) {
      if (i[d] == d_data[d][0].getBox().lower()[d]) {
         return d_data[d][0](i, depth);
      }
      if (i[d] == d_data[d][1].getBox().upper()[d]) {
         return d_data[d][1](i, depth);
      }
   }
   /*
    * The following lines should only be executed if there's a bug
    * in the Outernode datatype.
    */
   TBOX_ERROR("Bad index used to access outernode data\n"
      << "Given index is not an outernode of this instance.\n");
   return d_data[0][0](i, depth);
}

}
}
