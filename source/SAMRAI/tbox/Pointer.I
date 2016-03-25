/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer template class with RTTI
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::Pointer():
   d_object(0),
   d_counter(0)
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::Pointer(
   const Pointer<TYPE>& ptr):
   PointerBase(),
   d_object(ptr.d_object),
   d_counter(ptr.d_counter)
{
   if (d_counter) d_counter->addReference();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void Pointer<TYPE>::deleteObject()
{
   delete d_object;
   delete d_counter;

#ifdef DEBUG_INITIALIZE_UNDEFINED
   // For a little performance gain only do this when
   // debugging.
   d_object = 0;
   d_counter = 0;
#endif
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::~Pointer()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>& Pointer<TYPE>::operator = (
   const Pointer<TYPE>& ptr)
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
bool Pointer<TYPE>::operator == (
   const Pointer<TYPE>& rhs) const
{
   return d_object == rhs.d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool Pointer<TYPE>::operator != (
   const Pointer<TYPE>& rhs) const
{
   return d_object != rhs.d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void Pointer<TYPE>::setNull()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
   d_object = 0;
   d_counter = 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * Pointer<TYPE>::operator -> () const
{
   return d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& Pointer<TYPE>::operator * () const
{
   return *d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::operator TYPE * () const
{
   return d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool Pointer<TYPE>::isNull() const
{
   return !d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::operator bool () const
{
   return d_object ? true : false;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool Pointer<TYPE>::operator ! () const
{
   return !d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * Pointer<TYPE>::getPointer() const
{
   return d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ReferenceCounter * Pointer<TYPE>::getReferenceCounter() const
{
   return d_counter;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const DescribedClass * Pointer<TYPE>::getSubclassPointer() const
{
   // the explicit cast is needed by the brain-damaged SGI C++ compiler
   return (DescribedClass *)d_object;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ReferenceCounter * Pointer<TYPE>::getSubclassReferenceCounter() const
{
   return (ReferenceCounter *)d_counter;
}

}
}
