/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A simple doubly-linked list template class
 *
 ************************************************************************/

#ifndef included_tbox_List
#define included_tbox_List

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <cstddef>

namespace SAMRAI {
namespace tbox {

template<class TYPE>
class List;
template<class TYPE>
class ListNode;
template<class TYPE>
class ListIterator;

/**
 * Class List is a template-based doubly-linked list container.
 * It is simpler that the STL List class, which is not yet delivered
 * with all compilers.  In the future this class may be replaced with
 * the STL version.
 *
 * The List item in the list must define the assignment operator "=" and
 * the default const copy constructor.  This list implementation uses copy
 * to add the list item to the linked list.  If this is too expensive for
 * the object to be placed in the list, then a smart pointer to the
 * object should be used instead.
 */

template<class TYPE>
class List
{
public:
   /**
    * The iterator for class List.  Because of problems with some
    * compilers and nested templates, the iterator for List is a
    * separate class called ListIterator<TYPE> that is typedef-ed
    * to List<TYPE>::Iterator.  Use the List<TYPE>::Iterator
    * form, since the other form may disappear as compilers get better
    * and can accept nested template classes.
    */
   typedef ListIterator<TYPE> Iterator;

   /**
    * Create a list with no elements.
    */
   List();

   /**
    * The const constructor copies list items from the source list to
    * the newly constructed list.  The source list is not modified.
    */
   List(
      const List<TYPE>& list);

   /**
    * The assignment operator copies the values from the right hand side
    * list into the list on the left side of the equals sign.  The source
    * list is not modified.  All previous values in the destination list
    * are lost.
    */
   List<TYPE>&
   operator = (
      const List<TYPE>& list);

   /**
    * Deallocate a list.  If items remain in the list, they are deallocated.
    */
   virtual ~List();

   /**
    * Check whether a list is empty (has no elements).
    */
   bool
   isEmpty() const;

   /**
    * Return the number of items in the list.
    */
   int
   getNumberOfItems() const;

   /**
    * Return the number of items in the list.  Identical to getNumberOfItems(),
    * but this method is common to several container classes.
    */
   int
   size() const;

   /**
    * Add a new item to the head of the list.
    */
   void
   addItem(
      const TYPE& item);

   /**
    * Add a new item to the list before the item pointed to by the
    * specified iterator.
    */
   void
   addItemBefore(
      ListIterator<TYPE>& iter,
      const TYPE& item);

   /**
    * Add a new item to the list after the item pointed to by the
    * specified iterator.
    */
   void
   addItemAfter(
      ListIterator<TYPE>& iter,
      const TYPE& item);

   /**
    * Add a new item to the tail of the list.
    */
   void
   appendItem(
      const TYPE& item);

   /**
    * Copy list items from the argument list.  The new items are
    * appended to the end of the current list.  The argument list
    * is not modified.
    */
   void
   copyItems(
      const List<TYPE>& list);

   /**
    * Catenate list items from the argument list to the end of the
    * current list.  Unlike copyItems(), the argument list is set
    * to null.
    */
   void
   catenateItems(
      List<TYPE>& list);

   /**
    * Catenate list items from the argument list to the front of the
    * current list.  Unlike copyItems(), the argument list is set
    * to null.
    */
   void
   catenateItemsAtFront(
      List<TYPE>& list);

   /**
    * Empty the list.  All list items are deallocated.
    */
   void
   clearItems();

   /**
    * Return a reference to the first item in the list.  This operation is
    * not defined if the list is empty.  This member function is const since
    * it cannot change the list, although the item in the list may change.
    */
   TYPE&
   getFirstItem() const;

   /**
    * Return a reference to the last item in the list.  This operation is
    * not defined if the list is empty.  This member function is const since
    * it cannot change the list, although the item in the list may change.
    */
   TYPE&
   getLastItem() const;

   /**
    * Remove (deallocate) the first item in the list.
    */
   void
   removeFirstItem();

   /**
    * Remove (deallocate) the last item in the list.
    */
   void
   removeLastItem();

   /**
    * Remove (deallocate) the item pointed to by the iterator.  After
    * performing the remove, the iterator becomes invalid.  Iterators
    * that point to other elements in the list remain valid after the
    * deletion of an item.
    */
   void
   removeItem(
      ListIterator<TYPE>& iter);

   /**
    * Reverse all elements in the list.  This operation is performed in-place
    * by swapping next and previous list element pointers.
    */
   void
   reverse();

   /**
    * Swap the contents with that of the argument.
    *
    * This method is O(1) complexity.
    */
   void
   swap(
      List& r);

   /**
    * Return a List<TYPE>::Iterator that points to the start of the list.
    */
   typename List<TYPE>::Iterator
   listStart() const;

   /**
    * Return a List<TYPE>::Iterator that points to the end of the list.
    */
   typename List<TYPE>::Iterator
   listEnd() const;

private:
   friend class ListIterator<TYPE>;

   /**
    * @brief Cleanup static caches.
    *
    * Free cached list nodes.  A cached list of list nodes is used to speed
    * allocation and deallocation of list items.  This call will return this
    * cached memory to the free memory pool.
    */
   static void
   finalizeCallback();

   int d_number_items;
   ListNode<TYPE>* d_list_head;
   ListNode<TYPE>* d_list_tail;

   /**
    * @brief Has shutdown handler been initialized.
    *
    * This should be checked and set in every ctor.
    */
   static bool s_initialized;

   /**
    * @brief Initialize shutdown handler.
    */
   static bool
   initializeHandler();

};

/**
 * Class ListNode holds items for the linked list.  This class should
 * be defined inside the List class, but nested template classes can
 * cause some compilers to barf.  List nodes should never be seen by the
 * user of a list.
 *
 * There are no public functions for ListNode.
 *
 * @see tbox::List
 */

template<class TYPE>
class ListNode
{
private:
   static void
   freeCachedListItems();

   friend class ListIterator<TYPE>;
   friend class List<TYPE>;

   ListNode(
      const TYPE& t,
      ListNode<TYPE> * n,
      ListNode<TYPE> * p);
   ~ListNode();

   void *
   operator new (
      size_t bytes);
   void
   operator delete (
      void* what);

   static ListNode<TYPE>* s_free_list;
   static int s_num_free;
   static int s_max_free;
   static bool s_registered_callback;

   TYPE d_item;
   ListNode<TYPE>* d_next;
   ListNode<TYPE>* d_prev;
};

/**
 * Class ListIterator provides methods for stepping through lists.
 * This class definition should be nested within List, but nested template
 * classes can cause some compilers to barf.  The user should access this
 * class through the name List<TYPE>::Iterator rather than through
 * ListIterator<TYPE>, since the implementation may change as compilers
 * mature.
 *
 * The list iterator should be used as follows:
 * \verbatim
 * List<TYPE> list;
 * ...
 * for (List<TYPE>::Iterator l(list); l; l++) {
 *    ... = l();
 * }
 * \endverbatim
 *
 * @see tbox::List
 */

template<class TYPE>
class ListIterator
{
public:
   /**
    * Default constructor for ListIterator.  This must be initialized
    * before it can be used to iterate over a list.
    */
   ListIterator();

   /**
    * Create a list iterator pointing to the beginning of the specified list.
    */
   explicit ListIterator(
      const List<TYPE>& list);

   /**
    * Copy constructor for ListIterator.  Make this iterator point to the
    * same list and item as the argument.
    */
   ListIterator(
      const ListIterator<TYPE>& iterator);

   /**
    * Make the list iterator point to the same list and position as
    * the argument.
    */
   ListIterator<TYPE>&
   operator = (
      const ListIterator<TYPE>& iterator);

   /**
    * Destructor for class ListIterator.
    */
   ~ListIterator();

   /**
    * Return the current item from the list.  The iterator must point
    * to a valid item in the list; otherwise, this function may barf.
    */
   TYPE&
   operator * ();

   /**
    * Return a const references to the current item from the list.
    * The iterator must point to a valid item in the list; otherwise,
    * this function may barf.
    */
   const TYPE&
   operator * () const;

   /**
    * Return the current item from the list.  The iterator must point
    * to a valid item in the list; otherwise, this function may barf.
    */
   TYPE&
   operator () ();

   /**
    * Return a const reference to the current item from the list.
    * The iterator must point to a valid item in the list; otherwise,
    * this function may barf.
    */
   const TYPE&
   operator () () const;

   /**
    * Return true if the iterator points to a valid list item.
    */
   operator bool () const;

#ifndef LACKS_BOOL_VOID_RESOLUTION
   /**
    * Return a non-NULL value if the iterator points to a valid list item.
    */
   operator const void
   * () const;
#endif

   /**
    * Return whether the iterator points to a valid list item.  This
    * operator mimics the !p operation applied to a pointer p.
    */
   bool
   operator ! () const;

   /**
    * Increment the iterator to point to the next element in the list.
    */
   void
   operator ++ (
      int);

   /**
    * Decrement the iterator to point to the previous element in the list.
    */
   void
   operator -- (
      int);

   /**
    * Rewind the iterator to point to the first element in the list.
    */
   void
   rewindIterator();

   /**
    * Fast forward the iterator to point to the last element in the list.
    */
   void
   fastforwardIterator();

   /**
    * Test two iterators for equality (pointing to the same list item).
    * This does not check whether two list items have the same value, but
    * whether they are the same item.  Thus, if a and b are iterators,
    * "a == b" (object sameness) is more restrictive than "a() == b()"
    * (object equality).
    */
   bool
   operator == (
      const ListIterator<TYPE>& iterator) const;

   /**
    * Test two iterators for inequality (pointing to different list items).
    * See the comments for operator== for more information about equality
    * for list items.
    */
   bool
   operator != (
      const ListIterator<TYPE>& iterator) const;

private:
   friend class List<TYPE>;

   ListIterator(
      List<TYPE>* list,
      ListNode<TYPE>* node);

   List<TYPE>* d_list;
   ListNode<TYPE>* d_node;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/List.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SAMRAI/tbox/List.C"
#endif

#endif
