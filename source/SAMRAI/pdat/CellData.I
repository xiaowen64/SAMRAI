/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated cell centered patch data type
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int CellData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ArrayData<TYPE>& CellData<TYPE>::getArrayData()
{
   return d_data;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const ArrayData<TYPE>& CellData<TYPE>::getArrayData() const
{
   return d_data;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * CellData<TYPE>::getPointer(
   int depth)
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data.getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * CellData<TYPE>::getPointer(
   int depth) const
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data.getPointer(depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& CellData<TYPE>::operator () (
   const CellIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data(i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE& CellData<TYPE>::operator () (
   const CellIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data(i, depth);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void CellData<TYPE>::fill(
   const TYPE& t,
   int d)
{

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data.fill(t, d);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void CellData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data.fill(t, box, d);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void CellData<TYPE>::fillAll(
   const TYPE& t)
{
   d_data.fillAll(t);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD void CellData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   d_data.fillAll(t, box);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD void CellData<TYPE>::copyOnBox(
   const CellData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);
   const hier::Box cell_box = CellGeometry::toCellBox(box);
   d_data.copy(src.getArrayData(), cell_box);
}

}
}
