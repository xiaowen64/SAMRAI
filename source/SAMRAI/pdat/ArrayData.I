/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Templated array data structure supporting patch data types
 *
 ************************************************************************/

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
ArrayData<TYPE>::getDim() const
{
   return d_dim;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const hier::Box&
ArrayData<TYPE>::getBox() const
{
   return d_box;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool
ArrayData<TYPE>::isInitialized() const
{
   return d_depth > 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int
ArrayData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int
ArrayData<TYPE>::getOffset() const
{
   return d_offset;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE*
ArrayData<TYPE>::getPointer(
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_array.getPointer(d * d_offset);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE*
ArrayData<TYPE>::getPointer(
   int d) const
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_array.getPointer(d * d_offset);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE&
ArrayData<TYPE>::operator () (
   const hier::Index& i,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   const int index = d_box.offset(i) + d * d_offset;

   TBOX_ASSERT((index >= 0) && (index < d_depth * d_offset));

   return d_array[index];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE&
ArrayData<TYPE>::operator () (
   const hier::Index& i,
   int d) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   const int index = d_box.offset(i) + d * d_offset;

   TBOX_ASSERT((index >= 0) && (index < d_depth * d_offset));

   return d_array[index];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
size_t
ArrayData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   return tbox::MemoryUtilities::align(box.size() * depth * sizeof(TYPE));
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int
ArrayData<TYPE>::getDataStreamSize(
   const hier::BoxContainer& boxes,
   const hier::IntVector& source_shift) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(source_shift);
#endif

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, source_shift);

   const int nelements = boxes.getTotalSizeOfBoxes();

   if (typeid(TYPE) == typeid(bool)) {
      return tbox::MessageStream::getSizeof<bool>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(char)) {
      return tbox::MessageStream::getSizeof<char>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(dcomplex)) {
      return tbox::MessageStream::getSizeof<dcomplex>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(double)) {
      return tbox::MessageStream::getSizeof<double>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(float)) {
      return tbox::MessageStream::getSizeof<float>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(int)) {
      return tbox::MessageStream::getSizeof<int>(d_depth * nelements);
   }

   TBOX_ERROR("ArrayData::getDataStreamSize() -- Invalid type" << std::endl);
   return 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool
ArrayData<TYPE>::canEstimateStreamSizeFromBox()
{
   if ((typeid(TYPE) == typeid(bool))
       | (typeid(TYPE) == typeid(char))
       | (typeid(TYPE) == typeid(double))
       | (typeid(TYPE) == typeid(float))
       | (typeid(TYPE) == typeid(int))
       | (typeid(TYPE) == typeid(dcomplex))) {
      return true;
   } else {
      return false;
   }
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool
ArrayData<TYPE>::isValid()
{
   return !d_box.isEmpty();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
