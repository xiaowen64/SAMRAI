//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/array/ArrayData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated array data structure supporting patch data types
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline 
const hier::Box<DIM>& ArrayData<DIM,TYPE>::getBox() const
{
   return(d_box);
}

template<int DIM, class TYPE> inline
bool ArrayData<DIM,TYPE>::isInitialized() const
{
   return(d_depth > 0);
}

template<int DIM, class TYPE> inline
int ArrayData<DIM,TYPE>::getDepth() const
{
   return(d_depth);
}

template<int DIM, class TYPE> inline
int ArrayData<DIM,TYPE>::getOffset() const
{
   return(d_offset);
}

template<int DIM, class TYPE> inline
TYPE *ArrayData<DIM,TYPE>::getPointer(int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   return(d_array.getPointer(d * d_offset));
}

template<int DIM, class TYPE> inline
const TYPE *ArrayData<DIM,TYPE>::getPointer(int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   return(d_array.getPointer(d * d_offset));
}

template<int DIM, class TYPE> inline
TYPE& ArrayData<DIM,TYPE>::operator()(const hier::Index<DIM>& i, 
                                      int d)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   const int index = d_box.offset(i) + d*d_offset;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((index >= 0) && (index < d_depth*d_offset));
#endif
   return(d_array[index]);
}

template<int DIM, class TYPE> inline
const TYPE&
ArrayData<DIM,TYPE>::operator()(const hier::Index<DIM>& i, 
                                int d) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((d >= 0) && (d < d_depth));
#endif
   const int index = d_box.offset(i) + d*d_offset;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((index >= 0) && (index < d_depth*d_offset));
#endif
   return(d_array[index]);
}


template<int DIM, class TYPE> inline
bool ArrayData<DIM,TYPE>::isStandardType()
{
   if ( ( typeid(TYPE) == typeid(bool) ) | 
	( typeid(TYPE) == typeid(char) ) | 
	( typeid(TYPE) == typeid(double) ) | 
	( typeid(TYPE) == typeid(float) ) | 
	( typeid(TYPE) == typeid(int) ) ) {
      return true;
   } else {
      return false;
   }
}

template<int DIM, class TYPE>  inline
size_t ArrayData<DIM,TYPE>::getSizeOfData(
   const hier::Box<DIM>& box, 
   int depth)
{
   return(tbox::Arena::align(box.size() * depth * sizeof(TYPE)));
}

template<int DIM, class TYPE> inline
int ArrayData<DIM,TYPE>::getDataStreamSize(
   const hier::BoxList<DIM>& boxes,
   const hier::IntVector<DIM>& source_shift) const
{
   (void) source_shift;

   const int nelements = boxes.getTotalSizeOfBoxes();

   if (typeid(TYPE) == typeid(bool)) {
      return(tbox::AbstractStream::sizeofBool(d_depth * nelements));
   } else if (typeid(TYPE) == typeid(char)) {
      return(tbox::AbstractStream::sizeofChar(d_depth * nelements));
   } else if (typeid(TYPE) == typeid(dcomplex)) {
      return(tbox::AbstractStream::sizeofDoubleComplex(d_depth * nelements));
   } else if (typeid(TYPE) == typeid(double)) {
      return(tbox::AbstractStream::sizeofDouble(d_depth * nelements));
   } else if (typeid(TYPE) == typeid(float)) {
      return(tbox::AbstractStream::sizeofFloat(d_depth * nelements));
   } else if (typeid(TYPE) == typeid(int)) {
      return(tbox::AbstractStream::sizeofInt(d_depth * nelements));
   } 

   TBOX_ERROR("ArrayData::getDataStreamSize() -- Invalid type" << std::endl);
   return 0;
}

template<int DIM, class TYPE> inline
bool ArrayData<DIM,TYPE>::canEstimateStreamSizeFromBox()
{
  if ( ( typeid(TYPE) == typeid(bool) ) |
       ( typeid(TYPE) == typeid(char) ) |
       ( typeid(TYPE) == typeid(double) ) |
       ( typeid(TYPE) == typeid(float) ) |
       ( typeid(TYPE) == typeid(int) ) |
       ( typeid(TYPE) == typeid(dcomplex) ) ) {
      return true;
   } else {
      return false;
   }
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
