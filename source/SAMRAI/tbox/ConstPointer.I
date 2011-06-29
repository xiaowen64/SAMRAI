/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart const pointer template class with RTTI 
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ConstPointer<TYPE>::ConstPointer()
{
   d_object = (const TYPE *)NULL;
   d_counter = (ReferenceCounter *)NULL;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ConstPointer<TYPE>::ConstPointer(
   const ConstPointer<TYPE>& ptr):
   ConstPointerBase()
{
   d_object = ptr.d_object;
   d_counter = ptr.d_counter;
   if (d_counter) d_counter->addReference();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ConstPointer<TYPE>::~ConstPointer()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ConstPointer<TYPE>&
ConstPointer<TYPE>::operator = (
   const ConstPointer<TYPE>& ptr)
{
   if (this != &ptr) {
      if (d_counter && d_counter->deleteReference()) deleteObject();
      d_object = ptr.d_object;
      d_counter = ptr.d_counter;
      if (d_counter) d_counter->addReference();
   }
   return *this;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool
ConstPointer<TYPE>::operator == (
   const ConstPointer<TYPE>& rhs) const
{
   return d_object == rhs.d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool
ConstPointer<TYPE>::operator != (
   const ConstPointer<TYPE>& rhs) const
{
   return d_object != rhs.d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void ConstPointer<TYPE>::setNull()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
   d_object = (TYPE *)NULL;
   d_counter = (ReferenceCounter *)NULL;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * ConstPointer<TYPE>::operator -> () const
{
   return d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &ConstPointer<TYPE>::operator * () const
{
   return *d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ConstPointer<TYPE>::operator const TYPE* () const
{
   return d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool ConstPointer<TYPE>::isNull() const
{
   return !d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ConstPointer<TYPE>::operator bool () const
{
   return d_object ? true : false;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool ConstPointer<TYPE>::operator ! () const
{
   return !d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE * ConstPointer<TYPE>::getPointer() const
{
   return d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ReferenceCounter * ConstPointer<TYPE>::getReferenceCounter() const
{
   return d_counter;
}

}
}
