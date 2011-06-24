/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   A simple array template class 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Array<TYPE>::Array()
{
   d_objects = (TYPE *)NULL;
   d_counter = (ReferenceCounter *)NULL;
   d_elements = 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Array<TYPE>::Array(
   const Array<TYPE>& rhs)
{
   d_objects = rhs.d_objects;
   d_counter = rhs.d_counter;
   d_elements = rhs.d_elements;
   if (d_counter) d_counter->addReference();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Array<TYPE>::~Array()
{
   if (d_counter && d_counter->deleteReference()) deleteObjects();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & Array<TYPE>::operator [] (
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < d_elements));

   return d_objects[i];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &Array<TYPE>::operator [] (
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < d_elements));

   return d_objects[i];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void Array<TYPE>::setNull()
{
   if (d_counter && d_counter->deleteReference()) deleteObjects();
   d_objects = (TYPE *)NULL;
   d_counter = (ReferenceCounter *)NULL;
   d_elements = 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void Array<TYPE>::clear()
{
   if (d_counter && d_counter->deleteReference()) deleteObjects();
   d_objects = (TYPE *)NULL;
   d_counter = (ReferenceCounter *)NULL;
   d_elements = 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool Array<TYPE>::isNull() const
{
   return !d_objects;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool Array<TYPE>::empty() const
{
   return !d_objects;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * Array<TYPE>::getPointer(
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < d_elements));

   return &d_objects[i];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * Array<TYPE>::getPointer(
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < d_elements));

   return &d_objects[i];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int Array<TYPE>::getSize() const
{
   return d_elements;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int Array<TYPE>::size() const
{
   return d_elements;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
size_t Array<TYPE>::align(
   const size_t bytes)
{
   size_t aligned = bytes + ALLOCATION_ALIGNMENT - 1;
   aligned -= aligned % ALLOCATION_ALIGNMENT;
   return aligned;
}

}
}
