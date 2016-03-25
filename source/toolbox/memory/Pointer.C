//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/Pointer.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: A smart pointer template class with RTTI
//

#include <typeinfo>

#include "tbox/Pointer.h"
#include "tbox/Arena.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Pointer.I"
#endif

namespace SAMRAI {
   namespace tbox {


template <class TYPE>
Pointer<TYPE>::Pointer(TYPE *ptr, const bool managed)
{
   d_object = ptr;
   if (d_object && managed) {
      d_counter = new ReferenceCounter;
   } else {
      d_counter = (ReferenceCounter *) NULL;
   }
}

template <class TYPE>
Pointer<TYPE>::Pointer(TYPE *ptr,
                                 const Pointer<Arena>& pool)
{
   d_object = ptr;
   if (d_object) {
      d_counter = new ReferenceCounter(pool.getPointer(),
                                            pool.getReferenceCounter());
   } else {
      d_counter = (ReferenceCounter *) NULL;
   }
}

template <class TYPE>
Pointer<TYPE>::Pointer(const PointerBase& ptr)
{
   const DescribedClass *sub_ptr = ptr.getSubclassPointer();
   if(sub_ptr) {
      d_object = (TYPE *)dynamic_cast<const TYPE *>(sub_ptr);
   } else {
      d_object = NULL;
   }

   if (d_object) {
      d_counter = ptr.getSubclassReferenceCounter();
      if (d_counter) d_counter->addReference();
   } else {
      d_counter = (ReferenceCounter *) NULL;
   }
}

template <class TYPE>
Pointer<TYPE>& Pointer<TYPE>::operator=(TYPE *ptr)
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
   d_object = ptr;
   if (d_object) {
      d_counter = new ReferenceCounter;
   } else {
      d_counter = (ReferenceCounter *) NULL;
   }
   return(*this);
}

template <class TYPE>
Pointer<TYPE>& Pointer<TYPE>::operator=(const PointerBase& ptr)
{
   if (this != &ptr) {
      if (d_counter && d_counter->deleteReference()) deleteObject();


      const DescribedClass *sub_ptr = ptr.getSubclassPointer();
      if(sub_ptr) {
	 d_object = (TYPE *)dynamic_cast<const TYPE *>(sub_ptr);
      } else {
	 d_object = NULL;
      }

      if (d_object) {
	 d_counter = ptr.getSubclassReferenceCounter();
	 if (d_counter) d_counter->addReference();
      } else {
	 d_counter = (ReferenceCounter *) NULL;
      }
   }
   return(*this);
}

template <class TYPE>
void Pointer<TYPE>::deleteObject()
{
   if (d_counter->getArena()) {
      d_object->~TYPE();
      d_counter->getArena()->free(d_object);
   } else {
      delete d_object;
   }
   delete d_counter;

   d_object  = (TYPE *) NULL;
   d_counter = (ReferenceCounter *) NULL;
}

template <class TYPE>
const DescribedClass *Pointer<TYPE>::getSubclassPointer() const
{
   // the explicit cast is needed by the brain-damaged SGI C++ compiler
   return((DescribedClass *) d_object);
}

template <class TYPE>
ReferenceCounter *Pointer<TYPE>::getSubclassReferenceCounter() const
{
   return((ReferenceCounter *) d_counter);
}

}
}
