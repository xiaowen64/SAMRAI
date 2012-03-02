/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A simple doubly-linked list template class
 *
 ************************************************************************/

#ifndef included_tbox_List_C
#define included_tbox_List_C

#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

template<class TYPE>
ListNode<TYPE>* ListNode<TYPE>::s_free_list = NULL;
template<class TYPE>
int ListNode<TYPE>::s_num_free = 0;
template<class TYPE>
int ListNode<TYPE>::s_max_free = 1000;
template<class TYPE>
bool ListNode<TYPE>::s_registered_callback = false;

template<class TYPE>
bool List<TYPE>::s_initialized = false;

template<class TYPE>
ListNode<TYPE>::ListNode(
   const TYPE& t,
   ListNode<TYPE>* n,
   ListNode<TYPE>* p):
   d_item(t),
   d_next(n),
   d_prev(p)
{
}

template<class TYPE>
ListNode<TYPE>::~ListNode()
{
}

template<class TYPE>
void*
ListNode<TYPE>::operator new (
   size_t bytes)
{
   if (s_free_list) {
      ListNode<TYPE>* node = s_free_list;
      s_free_list = s_free_list->d_next;
      --s_num_free;
      return node;
   } else {
      return ::operator new (
                bytes);
   }
}

template<class TYPE>
void
ListNode<TYPE>::operator delete (
   void* what)
{
   ListNode<TYPE>* node = (ListNode<TYPE> *)what;
   node->d_next = s_free_list;
   s_free_list = node;
   if (++s_num_free > s_max_free) {
      int cut_to = s_max_free / 2;
      while (s_num_free > cut_to) {
         void * byebye = s_free_list;
         s_free_list = s_free_list->d_next;
         --s_num_free
         ;
         ::operator delete (
            byebye);
      }
   }
}

template<class TYPE>
void
ListNode<TYPE>::freeCachedListItems()
{
   while (s_free_list) {
      void * byebye = s_free_list;
      s_free_list = s_free_list->d_next
      ;
      ::operator delete (
         byebye);
   }
   s_free_list = NULL;
   s_num_free = 0;
}

template<class TYPE>
ListIterator<TYPE>::ListIterator()
{
   d_list = ((List<TYPE> *)NULL);
   d_node = ((ListNode<TYPE> *)NULL);
}

template<class TYPE>
ListIterator<TYPE>::ListIterator(
   const List<TYPE>& list)
{
   d_list = (List<TYPE> *) & list;
   d_node = list.d_list_head;
}

template<class TYPE>
ListIterator<TYPE>::ListIterator(
   List<TYPE>* list,
   ListNode<TYPE>* node)
{
   d_list = list;
   d_node = node;
}

template<class TYPE>
ListIterator<TYPE>::ListIterator(
   const ListIterator<TYPE>& iter)
{
   d_list = iter.d_list;
   d_node = iter.d_node;
}

template<class TYPE>
ListIterator<TYPE>::~ListIterator()
{
}

template<class TYPE>
ListIterator<TYPE>&
ListIterator<TYPE>::operator = (
   const ListIterator<TYPE>& iter)
{
   d_list = iter.d_list;
   d_node = iter.d_node;
   return *this;
}

template<class TYPE>
TYPE&
ListIterator<TYPE>::operator * ()
{
   return d_node->d_item;
}

template<class TYPE>
const TYPE&
ListIterator<TYPE>::operator * () const
{
   return d_node->d_item;
}

template<class TYPE>
TYPE&
ListIterator<TYPE>::operator () ()
{
   return d_node->d_item;
}

template<class TYPE>
const TYPE&
ListIterator<TYPE>::operator () () const
{
   return d_node->d_item;
}

template<class TYPE>
ListIterator<TYPE>::operator bool () const
{
   return d_node ? true : false;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<class TYPE>
ListIterator<TYPE>::operator const void * () const
{
   return d_node ? this : NULL;
}
#endif

template<class TYPE>
bool
ListIterator<TYPE>::operator ! () const
{
   return d_node ? false : true;
}

template<class TYPE>
void
ListIterator<TYPE>::operator ++ (
   int)
{
   if (d_node) d_node = d_node->d_next;
}

template<class TYPE>
void
ListIterator<TYPE>::operator -- (
   int)
{
   if (d_node) d_node = d_node->d_prev;
}

template<class TYPE>
void
ListIterator<TYPE>::rewindIterator()
{
   if (d_list) d_node = d_list->d_list_head;
}

template<class TYPE>
void
ListIterator<TYPE>::fastforwardIterator()
{
   if (d_list) d_node = d_list->d_list_tail;
}

template<class TYPE>
bool
ListIterator<TYPE>::operator == (
   const ListIterator<TYPE>& i) const
{
   return d_node == i.d_node;
}

template<class TYPE>
bool
ListIterator<TYPE>::operator != (
   const ListIterator<TYPE>& i) const
{
   return d_node != i.d_node;
}

template<class TYPE>
List<TYPE>::List(
   const List<TYPE>& list):
   d_number_items(0),
   d_list_head((ListNode<TYPE> *)NULL),
   d_list_tail((ListNode<TYPE> *)NULL)
{
   if (!s_initialized) {
      s_initialized = List<TYPE>::initializeHandler();
   }

   copyItems(list);
}

template<class TYPE>
List<TYPE>::List():
   d_number_items(0),
   d_list_head((ListNode<TYPE> *)NULL),
   d_list_tail((ListNode<TYPE> *)NULL)
{
   if (!s_initialized) {
      s_initialized = List<TYPE>::initializeHandler();
   }
}

template<class TYPE>
List<TYPE>::~List()
{
   clearItems();
}

template<class TYPE>
List<TYPE>&
List<TYPE>::operator = (
   const List<TYPE>& list)
{
   if (this != &list) {
      clearItems();
      copyItems(list);
   }
   return *this;
}

template<class TYPE>
void
List<TYPE>::addItem(
   const TYPE& item)
{

   ListNode<TYPE>* new_item =
      new ListNode<TYPE>(item, d_list_head, NULL);
   if (d_list_head) d_list_head->d_prev = new_item;
   d_list_head = new_item;
   if (!d_list_tail) d_list_tail = new_item;
   d_number_items++;
}

template<class TYPE>
void
List<TYPE>::addItemBefore(
   ListIterator<TYPE>& iter,
   const TYPE& item)
{
   if ((iter.d_list == this) && iter.d_node) {
      ListNode<TYPE>* new_item =
         new ListNode<TYPE>(item, iter.d_node, iter.d_node->d_prev);
      if (iter.d_node->d_prev == NULL) {
         d_list_head = new_item;
      } else {
         iter.d_node->d_prev->d_next = new_item;
      }
      iter.d_node->d_prev = new_item;
      d_number_items++;
   } else {
      addItem(item);
   }
}

template<class TYPE>
void
List<TYPE>::addItemAfter(
   ListIterator<TYPE>& iter,
   const TYPE& item)
{
   if ((iter.d_list == this) && iter.d_node) {
      ListNode<TYPE>* new_item =
         new ListNode<TYPE>(item, iter.d_node->d_next, iter.d_node);
      if (iter.d_node->d_next == NULL) {
         d_list_tail = new_item;
      } else {
         iter.d_node->d_next->d_prev = new_item;
      }
      iter.d_node->d_next = new_item;
      d_number_items++;
   } else {
      appendItem(item);
   }
}

template<class TYPE>
void
List<TYPE>::appendItem(
   const TYPE& item)
{
   ListNode<TYPE>* new_item =
      new ListNode<TYPE>(item, NULL, d_list_tail);
   if (d_list_tail) d_list_tail->d_next = new_item;
   d_list_tail = new_item;
   if (!d_list_head) d_list_head = new_item;
   d_number_items++;
}

template<class TYPE>
void
List<TYPE>::copyItems(
   const List<TYPE>& list)
{
   for (Iterator l(list); l; l++) {
      appendItem(l());
   }
}

template<class TYPE>
void
List<TYPE>::catenateItems(
   List<TYPE>& list)
{
   if (!list.isEmpty()) {
      if (isEmpty()) {
         d_list_head = list.d_list_head;
         d_list_tail = list.d_list_tail;
      } else {
         d_list_tail->d_next = list.d_list_head;
         list.d_list_head->d_prev = d_list_tail;
         d_list_tail = list.d_list_tail;
      }
      d_number_items += list.d_number_items;
      list.d_list_head = list.d_list_tail = NULL;
      list.d_number_items = 0;
   }
}

template<class TYPE>
void
List<TYPE>::catenateItemsAtFront(
   List<TYPE>& list)
{
   if (!list.isEmpty()) {
      if (isEmpty()) {
         d_list_head = list.d_list_head;
         d_list_tail = list.d_list_tail;
      } else {
         d_list_head->d_prev = list.d_list_tail;
         list.d_list_tail->d_next = d_list_head;
         d_list_head = list.d_list_head;
      }
      d_number_items += list.d_number_items;
      list.d_list_head = list.d_list_tail = NULL;
      list.d_number_items = 0;
   }
}

template<class TYPE>
void
List<TYPE>::removeFirstItem()
{
   if (!isEmpty()) {
      ListNode<TYPE>* node = d_list_head;
      d_list_head = d_list_head->d_next;
      if ((--d_number_items) > 0) {
         d_list_head->d_prev = NULL;
      } else {
         d_list_tail = NULL;
      }
      delete node;
   }
}

template<class TYPE>
void
List<TYPE>::removeLastItem()
{
   if (!isEmpty()) {
      ListNode<TYPE>* node = d_list_tail;
      d_list_tail = d_list_tail->d_prev;
      if ((--d_number_items) > 0) {
         d_list_tail->d_next = NULL;
      } else {
         d_list_head = NULL;
      }
      delete node;
   }
}

template<class TYPE>
void
List<TYPE>::clearItems()
{
   while (d_list_head) {
      ListNode<TYPE>* byebye = d_list_head;
      d_list_head = d_list_head->d_next;
      delete byebye;
   }
   d_list_head = d_list_tail = NULL;
   d_number_items = 0;
}

template<class TYPE>
void
List<TYPE>::removeItem(
   ListIterator<TYPE>& iter)
{
   if ((iter.d_list == this) && iter.d_node) {

      if ((d_list_head == iter.d_node) && (d_list_tail == iter.d_node)) {
         d_list_head = d_list_tail = NULL;

      } else if (d_list_head == iter.d_node) {
         d_list_head = iter.d_node->d_next;
         iter.d_node->d_next->d_prev = NULL;

      } else if (d_list_tail == iter.d_node) {
         d_list_tail = iter.d_node->d_prev;
         iter.d_node->d_prev->d_next = NULL;

      } else {
         iter.d_node->d_next->d_prev = iter.d_node->d_prev;
         iter.d_node->d_prev->d_next = iter.d_node->d_next;
      }

      d_number_items--;
      delete iter.d_node;
      iter.d_node = NULL;
   }
}

template<class TYPE>
void
List<TYPE>::reverse()
{
   ListNode<TYPE>* ptr = d_list_head;
   while (ptr) {
      ListNode<TYPE>* next = ptr->d_next;
      ptr->d_next = ptr->d_prev;
      ptr->d_prev = next;
      ptr = next;
   }
   ptr = d_list_head;
   d_list_head = d_list_tail;
   d_list_tail = ptr;
}

template<class TYPE>
void
List<TYPE>::swap(
   List& r)
{
   int tmpi = d_number_items;
   d_number_items = r.d_number_items;
   r.d_number_items = tmpi;

   ListNode<TYPE>* tmph = d_list_head;
   d_list_head = r.d_list_head;
   r.d_list_head = tmph;

   ListNode<TYPE>* tmpt = d_list_tail;
   d_list_tail = r.d_list_tail;
   r.d_list_tail = tmpt;
}

/*
 * Initialize shutdown handler.
 */
template<class TYPE>
bool
List<TYPE>::initializeHandler(
   void) {
   static StartupShutdownManager::Handler handler(
      0,
      0,
      0,
      List::finalizeCallback,
      StartupShutdownManager::priorityList);

   return true;
}

template<class TYPE>
bool
List<TYPE>::isEmpty() const
{
   return d_number_items == 0;
}

template<class TYPE>
int
List<TYPE>::getNumberOfItems() const
{
   return d_number_items;
}

template<class TYPE>
int
List<TYPE>::size() const
{
   return d_number_items;
}

template<class TYPE>
TYPE&
List<TYPE>::getFirstItem() const
{
   return d_list_head->d_item;
}

template<class TYPE>
TYPE&
List<TYPE>::getLastItem() const
{
   return d_list_tail->d_item;
}

template<class TYPE>
void
List<TYPE>::finalizeCallback()
{
   ListNode<TYPE>::freeCachedListItems();
}

template<class TYPE>
ListIterator<TYPE>
List<TYPE>::listStart() const
{
   return ListIterator<TYPE>((List<TYPE> *) this, d_list_head);
}

template<class TYPE>
ListIterator<TYPE>
List<TYPE>::listEnd() const
{
   return ListIterator<TYPE>((List<TYPE> *) this, d_list_tail);
}

}
}

#endif
