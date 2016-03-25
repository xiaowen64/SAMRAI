//
// File:	Array.I
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	A simple array template class
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


#ifndef NULL
#define NULL (0)
#endif

template <class TYPE>
inline
Array<TYPE>::Array()
{
   d_standard_type = false;
   d_objects  = (TYPE *) NULL;
   d_counter  = (ReferenceCounter *) NULL;
   d_elements = 0;
}

template <class TYPE>
inline
Array<TYPE>::Array(const Array<TYPE>& rhs)
{
   d_standard_type = rhs.d_standard_type;
   d_objects  = rhs.d_objects;
   d_counter  = rhs.d_counter;
   d_elements = rhs.d_elements;
   if (d_counter) d_counter->addReference();
}

template <class TYPE>
inline
Array<TYPE>::~Array()
{
   if (d_counter && d_counter->deleteReference()) deleteObjects();
}

template <class TYPE>
inline
TYPE& Array<TYPE>::operator[](const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_elements));
#endif
   return(d_objects[i]);
}

template <class TYPE>
inline
const TYPE& Array<TYPE>::operator[](const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_elements));
#endif
   return(d_objects[i]);
}

template <class TYPE>
inline
void Array<TYPE>::setNull()
{
   if (d_counter && d_counter->deleteReference()) deleteObjects();
   d_objects  = (TYPE *) NULL;
   d_counter  = (ReferenceCounter *) NULL;
   d_elements = 0;
}

template <class TYPE>
inline
bool Array<TYPE>::isNull() const
{
   return(!d_objects);
}

template <class TYPE>
inline
TYPE *Array<TYPE>::getPointer(const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_elements));
#endif
   return(&d_objects[i]);
}

template <class TYPE>
inline
const TYPE *Array<TYPE>::getPointer(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((i >= 0) && (i < d_elements));
#endif
   return(&d_objects[i]);
}

template <class TYPE>
inline
int Array<TYPE>::getSize() const
{
   return(d_elements);
}

template <class TYPE>
inline
int Array<TYPE>::size() const
{
   return(d_elements);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

