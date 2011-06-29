/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A simple doubly-linked list template class 
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
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
SAMRAI_INLINE_KEYWORD
ListNode<TYPE>::~ListNode()
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::ListIterator()
{
   d_list = ((List<TYPE> *)NULL);
   d_node = ((ListNode<TYPE> *)NULL);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::ListIterator(
   const List<TYPE>& list)
{
   d_list = (List<TYPE> *) & list;
   d_node = list.d_list_head;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::ListIterator(
   List<TYPE>* list,
   ListNode<TYPE>* node)
{
   d_list = list;
   d_node = node;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::ListIterator(
   const ListIterator<TYPE>& iter)
{
   d_list = iter.d_list;
   d_node = iter.d_node;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>&
ListIterator<TYPE>::operator = (
   const ListIterator<TYPE>& iter)
{
   d_list = iter.d_list;
   d_node = iter.d_node;
   return *this;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::~ListIterator()
{
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & ListIterator<TYPE>::operator * ()
{
   return d_node->d_item;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &ListIterator<TYPE>::operator * () const
{
   return d_node->d_item;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE & ListIterator<TYPE>::operator () ()
{
   return d_node->d_item;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const TYPE &ListIterator<TYPE>::operator () () const
{
   return d_node->d_item;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::operator bool () const
{
   return d_node ? true : false;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE>::operator const void* () const
{
   return d_node ? this : NULL;
}
#endif

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool ListIterator<TYPE>::operator ! () const
{
   return d_node ? false : true;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void ListIterator<TYPE>::operator ++ (
   int)
{
   if (d_node) d_node = d_node->d_next;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void ListIterator<TYPE>::operator -- (
   int)
{
   if (d_node) d_node = d_node->d_prev;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void ListIterator<TYPE>::rewindIterator()
{
   if (d_list) d_node = d_list->d_list_head;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void ListIterator<TYPE>::fastforwardIterator()
{
   if (d_list) d_node = d_list->d_list_tail;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool ListIterator<TYPE>::operator == (
   const ListIterator<TYPE>& i) const
{
   return d_node == i.d_node;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool ListIterator<TYPE>::operator != (
   const ListIterator<TYPE>& i) const
{
   return d_node != i.d_node;
}

/**
 * @brief Initialize shutdown handler.
 */
template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool List<TYPE>::initializeHandler(
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
SAMRAI_INLINE_KEYWORD
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
SAMRAI_INLINE_KEYWORD
List<TYPE>::~List()
{
   clearItems();
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
bool List<TYPE>::isEmpty() const
{
   return d_number_items == 0;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int List<TYPE>::getNumberOfItems() const
{
   return d_number_items;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int List<TYPE>::size() const
{
   return d_number_items;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& List<TYPE>::getFirstItem() const
{
   return d_list_head->d_item;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
TYPE& List<TYPE>::getLastItem() const
{
   return d_list_tail->d_item;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE> List<TYPE>::listStart() const
{
   return ListIterator<TYPE>((List<TYPE> *) this, d_list_head);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
ListIterator<TYPE> List<TYPE>::listEnd() const
{
   return ListIterator<TYPE>((List<TYPE> *) this, d_list_tail);
}

}
}
