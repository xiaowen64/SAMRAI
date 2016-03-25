//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/ConstPointer.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: A smart const pointer template class with RTTI
//


#include <typeinfo>

#include "tbox/ConstPointer.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/ConstPointer.I"
#endif

namespace SAMRAI {
   namespace tbox {


template <class TYPE>
ConstPointer<TYPE>::ConstPointer(const TYPE *ptr, const bool managed)
{
   d_object = ptr;
   if (d_object && managed) {
      d_counter = new ReferenceCounter;
   } else {
      d_counter = (ReferenceCounter *) NULL;
   }
}

template <class TYPE>
ConstPointer<TYPE>::ConstPointer(
   const TYPE *ptr,
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
ConstPointer<TYPE>::ConstPointer(const ConstPointerBase& ptr)
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
ConstPointer<TYPE>& ConstPointer<TYPE>::operator=(const TYPE *ptr)
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
ConstPointer<TYPE>&
ConstPointer<TYPE>::operator=(const ConstPointerBase& ptr)
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
void ConstPointer<TYPE>::deleteObject()
{
   if (d_counter->getArena()) {
      d_object->~TYPE();
      d_counter->getArena()->free((void *) d_object);
   } else {
      delete d_object;
   }
   delete d_counter;

   d_object  = (TYPE *) NULL;
   d_counter = (ReferenceCounter *) NULL;
}

template <class TYPE>
const DescribedClass *ConstPointer<TYPE>::getSubclassPointer() const
{
    return((const DescribedClass *)d_object);
}


template <class TYPE>
ReferenceCounter *
ConstPointer<TYPE>::getSubclassReferenceCounter() const
{
   return((ReferenceCounter *) d_counter);
}

}
}
