//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/Pointer.I $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: A smart pointer template class with RTTI
//

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
Pointer<TYPE>::Pointer()
{
   d_object  = (TYPE *) NULL;
   d_counter = (ReferenceCounter *) NULL;
}

template <class TYPE>
inline
Pointer<TYPE>::Pointer(const Pointer<TYPE>& ptr)
{
   d_object  = ptr.d_object;
   d_counter = ptr.d_counter;
   if (d_counter) d_counter->addReference();
}

template <class TYPE>
inline
Pointer<TYPE>::~Pointer()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
}

template <class TYPE>
inline
Pointer<TYPE>& Pointer<TYPE>::operator=(const Pointer<TYPE>& ptr)
{
   if (this != &ptr) {
      if (d_counter && d_counter->deleteReference()) deleteObject();
      d_object  = ptr.d_object;
      d_counter = ptr.d_counter;
      if (d_counter) d_counter->addReference();
   }
   return(*this);
}

template <class TYPE>
inline
bool Pointer<TYPE>::operator==(const Pointer<TYPE>& rhs) const
{
   return(d_object == rhs.d_object);
}

template <class TYPE>
inline
bool Pointer<TYPE>::operator!=(const Pointer<TYPE>& rhs) const
{
   return(d_object != rhs.d_object);
}

template <class TYPE>
inline
void Pointer<TYPE>::setNull()
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
   d_object  = (TYPE *) NULL;
   d_counter = (ReferenceCounter *) NULL;
}

template <class TYPE>
inline
TYPE *Pointer<TYPE>::operator->() const
{
   return(d_object);
}

template <class TYPE>
inline
TYPE& Pointer<TYPE>::operator*() const
{
   return(*d_object);
}

template <class TYPE>
inline
Pointer<TYPE>::operator TYPE *() const
{
   return(d_object);
}

template <class TYPE>
inline
bool Pointer<TYPE>::isNull() const
{
   return(!d_object);
}

template <class TYPE>
inline
Pointer<TYPE>::operator bool() const
{
   return(d_object ? true : false);
}

template <class TYPE>
inline
bool Pointer<TYPE>::operator!() const
{
   return(!d_object);
}

template <class TYPE>
inline
TYPE *Pointer<TYPE>::getPointer() const
{
   return(d_object);
}

template <class TYPE>
inline
ReferenceCounter *Pointer<TYPE>::getReferenceCounter() const
{
   return(d_counter);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
