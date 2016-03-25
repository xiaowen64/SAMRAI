//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/Array.C $
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2037 $
// Modified:	$LastChangedDate: 2008-03-05 15:54:45 -0800 (Wed, 05 Mar 2008) $
// Description:	A simple array template class
//

#include "tbox/Array.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Array.I"
#endif

namespace SAMRAI {
   namespace tbox {


/*
 * Default assume Array is not a standard type
 */
template <class TYPE> const bool Array<TYPE>::s_standard_type = false;

template <class TYPE>
Array<TYPE>::Array(const int n)
{
   if (n > 0) {
      d_objects  = new TYPE[n];
      d_counter  = new ReferenceCounter;
      d_elements = n;
   } else {
      d_objects  = (TYPE *) NULL;
      d_counter  = (ReferenceCounter *) NULL;
      d_elements = 0;
   }
}

template <class TYPE>
Array<TYPE>::Array(const int n, const Pointer<Arena>& pool)
{
   if (n > 0) {
      Arena *arena = pool.getPointer();
      if (!arena) {
         d_objects  = new TYPE[n];
         d_counter  = new ReferenceCounter;
         d_elements = n;
      } else {
         d_objects  = allocateObjects(n, arena);
         d_counter  =
            new ReferenceCounter(arena, pool.getReferenceCounter());
         d_elements = n;
      }
   } else {
      d_objects  = (TYPE *) NULL;
      d_counter  = (ReferenceCounter *) NULL;
      d_elements = 0;
   }
}

template <class TYPE>
Array<TYPE>& Array<TYPE>::operator=(const Array<TYPE>& rhs)
{
   if (this != &rhs) {
      if (d_counter && d_counter->deleteReference()) deleteObjects();
      d_objects  = rhs.d_objects;
      d_counter  = rhs.d_counter;
      d_elements = rhs.d_elements;
      if (d_counter) d_counter->addReference();
   }
   return(*this);
}

template <class TYPE>
TYPE *Array<TYPE>::allocateObjects(const int n, Arena *arena)
{
   TYPE *ptr = (TYPE *) ::operator new(n*sizeof(TYPE), arena);
   if (!s_standard_type) {
      for (int i = 0; i < n; i++) {
         (void) new (&ptr[i]) TYPE;
      }
   }
   return(ptr);
}

template <class TYPE>
void Array<TYPE>::deleteObjects()
{
   if (d_counter->getArena()) {
      if (!s_standard_type) {
         for (int i = 0; i < d_elements; i++) {
            d_objects[i].~TYPE();
         }
      }
      d_counter->getArena()->free(d_objects);
   } else {
      delete [] d_objects;
   }
   delete d_counter;

   d_objects  = (TYPE *) NULL;
   d_counter  = (ReferenceCounter *) NULL;
   d_elements = 0;
}

template <class TYPE>
void Array<TYPE>::resizeArray(const int n)
{
   if (n != d_elements) {
      Array<TYPE> array(n);
      const int s = (d_elements < n ? d_elements : n);
      for (int i = 0; i < s; i++) {
         array.d_objects[i] = d_objects[i];
      }
      this->operator=(array);
   }
}

template <class TYPE>
void Array<TYPE>::resizeArray(
   const int n, const Pointer<Arena>& pool)
{
   Array<TYPE> array(n, pool);
   const int s = (d_elements < n ? d_elements : n);
   for (int i = 0; i < s; i++) {
      array.d_objects[i] = d_objects[i];
   }
   this->operator=(array);
}

}
}
