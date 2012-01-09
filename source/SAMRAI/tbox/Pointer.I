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
   d_object(ptr.d_object),
   d_counter(ptr.d_counter)
{
   if (d_counter) d_counter->addReference();
}

template<typename TYPE>
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

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::~Pointer()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
void Pointer<TYPE>::reset()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
   d_object = 0;
   d_counter = 0;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * Pointer<TYPE>::operator -> () const
{
   return d_object;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& Pointer<TYPE>::operator * () const
{
   return *d_object;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
Pointer<TYPE>::operator bool () const
{
   return d_object ? true : false;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
TYPE * Pointer<TYPE>::get() const
{
   return d_object;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
ReferenceCounter * Pointer<TYPE>::getReferenceCounter() const
{
   return d_counter;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
const DescribedClass * Pointer<TYPE>::getSubclassPointer() const
{
   // the explicit cast is needed by the brain-damaged SGI C++ compiler
   return (DescribedClass *)d_object;
}

template<typename TYPE>
SAMRAI_INLINE_KEYWORD
ReferenceCounter * Pointer<TYPE>::getSubclassReferenceCounter() const
{
   return (ReferenceCounter *)d_counter;
}

}
}
