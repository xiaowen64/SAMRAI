/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer template class with RTTI
 *
 ************************************************************************/

#ifndef included_tbox_Pointer
#define included_tbox_Pointer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/ReferenceCounter.h"
#include "SAMRAI/tbox/PointerBase.h"

namespace SAMRAI {
namespace tbox {

/**
 * Class Pointer<TYPE> defines a smart pointer to TYPE.  It frees the
 * user from explicitly deleting and tracking aliases for object pointers.
 * It manages all reference counting and deallocation of the pointer (even
 * if the data was originally allocated from a memory arena).  When the
 * reference count on a Pointer<TYPE> object goes to zero, the object
 * is automatically deallocated.  A block with a references count and arena
 * pointer is allocated for all non-NULL pointers.  These reference counted
 * blocks are freed at the end of the lifetime of the pointer.
 *
 * Non-const pointers can be created only from other non-const pointers.
 * The non-const and const pointer classes have been designed so that an
 * attempted conversion from a const pointer into a non-const pointer causes
 * a compile-time error.  The const convention for this class is that any
 * member function that does not change the pointer object itself may be
 * const.  Note that this means that a const pointer does not mean that the
 * pointed-to object is const; these semantics require the const pointer class.
 *
 * Class Pointer<TYPE> performs type-checking when assigning pointers
 * of different TYPEs.  If a bad type conversion is performed, then the
 * destination pointer is set to NULL.
 *
 * @see tbox::DescribedClass
 * @see tbox::Array
 * @see tbox::ConstPointer
 * @see tbox::PointerBase
 * @see tbox::ReferenceCounter
 */

template<class TYPE>
class Pointer:public PointerBase
{
public:
   /**
    * The default constructor creates a null pointer.
    */
   Pointer();

   /**
    * Create a smart pointer with value ptr.  If managed is true, then
    * deallocation of the object pointed to by ptr will be taken care of
    * by the smart pointer.  This form assumes the pointer was allocated
    * using the standard new operator.
    */
   explicit Pointer(
      TYPE * ptr,
      const bool managed = true);

   /**
    * The pointer const constructor creates a smart pointer reference
    * aliased to the argument.
    */
   Pointer(
      const Pointer<TYPE>& ptr);

   /**
    * Create a pointer by attempting to type-cast the argument to TYPE.
    * If the type-cast fails, then the destination pointer will be set
    * to NULL.
    */
   Pointer(
      const PointerBase& ptr);

   /**
    * The pointer destructor frees the pointer data if the reference
    * count drops to zero.  The object is deallocated from the memory
    * pool (if it was specified in the constructor call).
    */
   ~Pointer();

   /**
    * Smart pointer assignment.  The left hand side points to the
    * right hand side and the reference count is incremented by one.
    */
   Pointer<TYPE>&
   operator = (
      const Pointer<TYPE>& ptr);

   /**
    * Create a managed smart pointer with value ptr.  The object pointed
    * to by ptr will be deallocated via delete when the reference count
    * goes to zero.
    */
   Pointer<TYPE>&
   operator = (
      TYPE * ptr);

   /**
    * Attempt to convert the argument pointer to a Pointer<TYPE>.
    * If the type conversion fails, then the destination pointer will be
    * set to NULL.
    */
   Pointer<TYPE>&
   operator = (
      const PointerBase& ptr);

   /**
    * Check whether two smart pointers point to the same object.
    */
   bool
   operator == (
      const Pointer<TYPE>& rhs) const;

   /**
    * Check whether two smart pointers point to different objects.
    */
   bool
   operator != (
      const Pointer<TYPE>& rhs) const;

   /**
    * Delegate member operations to the pointed-to object.  C++ defines
    * the ``->'' operator in a funny way to support delegation.  The
    * statement ptr->foo() acts as if ptr where actually a pointer
    * to an object with member function foo() instead of a class that
    * holds that pointer.  This member function is const since it cannot
    * change the pointer, although the pointed-to object may change.
    */
   TYPE *
   operator -> () const;

   /**
    * Dereference the smart pointer.  This member function is const since
    * it cannot change the pointer, although the pointed-to object may
    * change.
    */
   TYPE&
   operator * () const;

   /**
    * Implicit conversion of the smart pointer to the pointed-to object.
    * This conversion operator is const since it cannot change the pointer,
    * although the pointed-to object may change.
    */
   operator TYPE
   * () const;

   /**
    * Explicitly convert the smart pointer to the pointed-to object.
    * This member function is const since it cannot change the pointer,
    * although the pointed-to object may change.
    */
   TYPE *
   getPointer() const;

   /**
    * Set the pointed-to object to ptr.
    *
    * NOTE: reset name comes from boost library.
    */
   void
   reset(
      TYPE * ptr);

   /**
    * Check whether the smart pointer points to NULL.
    */
   bool
   isNull() const;

   /**
    * Return true if the pointer is non-NULL.
    */
   operator bool () const;

   /**
    * Return true if the pointer is NULL and false otherwise.  This operator
    * mimics the semantics of !p applied to a (regular) pointer p.
    */
   bool
   operator ! () const;

   /**
    * Set the smart pointer to NULL.
    */
   void
   setNull();

   /**
    * Return a pointer to the internal reference counter.  This routine
    * should not be called by the casual user.
    */
   ReferenceCounter *
   getReferenceCounter() const;

private:
   /*!
    * @brief Unused conversion, meant to prevent unintended conversion
    * of Pointer to int, which can cause serious hidden bugs.
    */
   operator int () const;

   void
   deleteObject();
   ReferenceCounter *
   getSubclassReferenceCounter() const;

   const DescribedClass *
   getSubclassPointer() const;

   TYPE* d_object;
   ReferenceCounter* d_counter;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/Pointer.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SAMRAI/tbox/Pointer.C"
#endif

#endif
