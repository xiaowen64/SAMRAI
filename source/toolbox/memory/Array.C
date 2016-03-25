//
// File:	Array.C
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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


template <class TYPE>
Array<TYPE>::Array(const int n, const bool standard_type)
{
   d_standard_type = standard_type;
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
Array<TYPE>::Array(const int n, const Pointer<Arena>& pool,
                             const bool standard_type)
{
   d_standard_type = standard_type;
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
      d_standard_type = rhs.d_standard_type;
      if (d_counter) d_counter->addReference();
   }
   return(*this);
}

template <class TYPE>
TYPE *Array<TYPE>::allocateObjects(const int n, Arena *arena)
{
   TYPE *ptr = (TYPE *) ::operator new(n*sizeof(TYPE), arena);
   if (!d_standard_type) {
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
      if (!d_standard_type) {
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
      Array<TYPE> array(n, d_standard_type);
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
   Array<TYPE> array(n, pool, d_standard_type);
   const int s = (d_elements < n ? d_elements : n);
   for (int i = 0; i < s; i++) {
      array.d_objects[i] = d_objects[i];
   }
   this->operator=(array);
}

}
}
