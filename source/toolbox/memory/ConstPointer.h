//
// File:	ConstPointer.h
// Package:	SAMRAI toolbox for memory management
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A smart const pointer template class with RTTI
//

#ifndef included_tbox_ConstPointer
#define included_tbox_ConstPointer

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_ConstPointerBase
#include "tbox/ConstPointerBase.h"
#endif
#ifndef included_tbox_ReferenceCounter
#include "tbox/ReferenceCounter.h"
#endif


namespace SAMRAI {
   namespace tbox {


class Arena;
template <class TYPE> class Pointer;

/**
 * Class ConstPointer<TYPE> defines a smart const pointer to TYPE.
 * It frees the user from explicitly deleting and tracking aliases for
 * object pointers.  It manages all reference counting and deallocation
 * of the pointer (even if the data was originally allocated from a memory
 * arena).  When the reference count on a ConstPointer<TYPE> object
 * goes to zero, the object is automatically deallocated.  A block with a
 * references count and arena pointer is allocated for all non-NULL pointers.
 * These reference counted blocks are freed at the end of the lifetime of
 * the pointer.
 *
 * Const pointers can be created from either const or non-const pointers.
 * The non-const and const pointer classes have been designed so that an
 * attempted conversion from a const pointer into a non-const pointer causes
 * a compile-time error.
 *
 * Class ConstPointer<TYPE> performs type-checking when assigning
 * pointers of different TYPEs.  If a bad type conversion is performed,
 * then the destination pointer is set to NULL.  
 *
 * @see tbox::Array
 * @see tbox::ConstPointerBase
 * @see tbox::ReferenceCounter
 * @see tbox::PointerBase
 * @see tbox::Pointer
 */

template <class TYPE>
class ConstPointer : public ConstPointerBase
{
public:
   /**
    * The default constructor creates a null pointer.
    */
   ConstPointer();

   /**
    * Create a smart pointer with value ptr.  If managed is true, then
    * deallocation of the object pointed to by ptr will be taken care of
    * by the smart pointer.  This form assumes the pointer was allocated
    * using the standard new operator.
    */
   ConstPointer(const TYPE *ptr, const bool managed = true);

   /**
    * Create a smart pointer for which the data was allocated from pool.
    * When the pointer is destroyed, the object is deallocated from the
    * specified memory pool.
    */
   ConstPointer(const TYPE *ptr, const Pointer<Arena>& pool);

   /**
    * The pointer const constructor creates a smart pointer reference
    * aliased to the argument.
    */
   ConstPointer(const ConstPointer<TYPE>& ptr);

   /**
    * Create a pointer by attempting to type-cast the argument to TYPE.
    * If the type-cast fails, then the destination pointer will be set
    * to NULL.
    */
   ConstPointer(const ConstPointerBase& ptr);

   /**
    * The pointer destructor frees the pointer data if the reference
    * count drops to zero.  The object is deallocated from the memory
    * pool (if it was specified in the constructor call).
    */
   ~ConstPointer();

   /**
    * Smart pointer assignment.  The left hand side points to the
    * right hand side and the reference count is incremented by one.
    */
   ConstPointer<TYPE>& operator=(const ConstPointer<TYPE>& ptr);

   /**
    * Create a managed smart pointer with value ptr.  The object pointed
    * to by ptr will be deallocated via delete when the reference count
    * goes to zero.
    */
   ConstPointer<TYPE>& operator=(const TYPE *ptr);

   /**
    * Attempt to convert the argument pointer to a ConstPointer<TYPE>.
    * If the type conversion fails, then the destination pointer will be set
    * to NULL.
    */
   ConstPointer<TYPE>& operator=(const ConstPointerBase& ptr);

   /**
    * Check whether two smart pointers point to the same object.
    */
   bool operator==(const ConstPointer<TYPE>& rhs) const;

   /**
    * Check whether two smart pointers point to different objects.
    */
   bool operator!=(const ConstPointer<TYPE>& rhs) const;

   /**
    * Delegate member operations to the pointed-to object.  C++ defines
    * the ``->'' operator in a funny way to support delegation.  The
    * statement ptr->foo() acts as if ptr where actually a pointer
    * to an object with member function foo() instead of a class that
    * holds that pointer.
    */
   const TYPE *operator->() const;

   /**
    * Dereference the smart pointer.  The pointer returned is a const
    * pointer to the object.
    */
   const TYPE& operator*() const;

   /**
    * Implicit conversion of the smart pointer to the pointed-to object.
    * The pointer returned is a const pointer to the object.
    */
   operator const TYPE *() const;

   /**
    * Explicitly convert the smart pointer to the pointed-to object.
    * The pointer returned is a const pointer to the object.
    */
   const TYPE *getPointer() const;

   /**
    * Check whether the smart pointer points to NULL.
    */
   bool isNull() const;

   /**
    * Return true if the pointer is non-NULL.
    */
   operator bool() const;

   /**
    * Return true if the pointer is NULL and false otherwise.  This operator
    * mimics the semantics of !p applied to a (regular) pointer p.
    */
   bool operator!() const;

   /**
    * Set the smart pointer to NULL.
    */
   void setNull();

   /**
    * Return a pointer to the internal reference counter.  This routine
    * should not be called by the casual user.
    */
   ReferenceCounter *getReferenceCounter() const;

private:
   void deleteObject();
   ReferenceCounter *getSubclassReferenceCounter() const;
   const DescribedClass *getSubclassPointer() const;

   const TYPE *d_object;
   ReferenceCounter *d_counter;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/ConstPointer.I"
#endif
#endif
