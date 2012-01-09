/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer template class with RTTI
 *
 ************************************************************************/

#ifndef included_tbox_Pointer_C
#define included_tbox_Pointer_C

#include <typeinfo>

#include "SAMRAI/tbox/Pointer.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/tbox/Pointer.I"
#endif

namespace SAMRAI {
namespace tbox {

template<typename TYPE>
Pointer<TYPE>::Pointer(
   TYPE* ptr,
   const bool managed):
   d_object(ptr)
{
   if (d_object && managed) {
      d_counter = new ReferenceCounter;
   } else {
      d_counter = (ReferenceCounter *)NULL;
   }
}

template<typename TYPE>
Pointer<TYPE>& Pointer<TYPE>::operator = (
   TYPE* ptr)
{
   if (d_object != ptr) {
      if (d_counter && d_counter->deleteReference()) deleteObject();
      d_object = ptr;
      if (d_object) {
         d_counter = new ReferenceCounter;
      } else {
         d_counter = (ReferenceCounter *)NULL;
      }
   }
   return *this;
}

template<typename TYPE>
void Pointer<TYPE>::reset(
   TYPE* ptr)
{
   if (d_counter && d_counter->deleteReference()) deleteObject();
   d_object = ptr;
   if (d_object) {
      d_counter = new ReferenceCounter;
   } else {
      d_counter = (ReferenceCounter *)NULL;
   }
}

}
}

#endif
