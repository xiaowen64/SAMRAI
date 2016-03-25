//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/base/List.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	A simple doubly-linked list template class
//

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {

template <class TYPE>
inline
ListNode<TYPE>::ListNode(const TYPE& t,
                                   ListNode<TYPE> *n,
                                   ListNode<TYPE> *p)
:  d_item(t), d_next(n), d_prev(p)
{
}

template <class TYPE>
inline
ListNode<TYPE>::~ListNode()
{
}

template <class TYPE>
inline
ListIterator<TYPE>::ListIterator()
{
   d_list = ((List<TYPE> *) NULL);
   d_node = ((ListNode<TYPE> *) NULL);
}

template <class TYPE>
inline
ListIterator<TYPE>::ListIterator(const List<TYPE>& list)
{
   d_list = (List<TYPE> *) &list;
   d_node = list.d_list_head;
}

template <class TYPE>
inline
ListIterator<TYPE>::ListIterator(List<TYPE>* list,
                                           ListNode<TYPE> *node)
{
   d_list = list;
   d_node = node;
}

template <class TYPE>
inline
ListIterator<TYPE>::ListIterator(const ListIterator<TYPE>& iter)
{
   d_list = iter.d_list;
   d_node = iter.d_node;
}

template <class TYPE>
inline
ListIterator<TYPE>&
ListIterator<TYPE>::operator=(const ListIterator<TYPE>& iter)
{
   d_list = iter.d_list;
   d_node = iter.d_node;
   return(*this);
}

template <class TYPE>
inline
ListIterator<TYPE>::~ListIterator()
{
}

template <class TYPE>
inline
TYPE& ListIterator<TYPE>::operator*()
{
   return(d_node->d_item);
}

template <class TYPE>
inline
const TYPE& ListIterator<TYPE>::operator*() const
{
   return(d_node->d_item);
}

template <class TYPE>
inline
TYPE& ListIterator<TYPE>::operator()()
{
   return(d_node->d_item);
}

template <class TYPE>
inline
const TYPE& ListIterator<TYPE>::operator()() const
{
   return(d_node->d_item);
}

template <class TYPE>
inline
ListIterator<TYPE>::operator bool() const
{
   return(d_node ? true : false);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template <class TYPE>
inline
ListIterator<TYPE>::operator const void*() const
{
   return(d_node ? this : NULL);
}
#endif

template <class TYPE>
inline
bool ListIterator<TYPE>::operator!() const
{
   return(d_node ? false : true);
}

template <class TYPE>
inline
void ListIterator<TYPE>::operator++(int)
{
   if (d_node) d_node = d_node->d_next;
}

template <class TYPE>
inline
void ListIterator<TYPE>::operator--(int)
{
   if (d_node) d_node = d_node->d_prev;
}

template <class TYPE>
inline
void ListIterator<TYPE>::rewindIterator()
{
   if (d_list) d_node = d_list->d_list_head;
}

template <class TYPE>
inline
void ListIterator<TYPE>::fastforwardIterator()
{
   if (d_list) d_node = d_list->d_list_tail;
}

template <class TYPE>
inline
bool ListIterator<TYPE>::operator==(const ListIterator<TYPE>& i) const
{
   return(d_node == i.d_node);
}

template <class TYPE>
inline
bool ListIterator<TYPE>::operator!=(const ListIterator<TYPE>& i) const
{
   return(d_node != i.d_node);
}

template <class TYPE>
inline
void List<TYPE>::freeCachedListItems()
{
   ListNode<TYPE>::freeCachedListItems();
}

template <class TYPE>
inline
List<TYPE>::List()
:  d_number_items(0),
   d_list_head((ListNode<TYPE> *) NULL),
   d_list_tail((ListNode<TYPE> *) NULL)
{
}

template <class TYPE>
inline
List<TYPE>::~List()
{
   clearItems();
}

template <class TYPE>
inline
bool List<TYPE>::isEmpty() const
{
   return(d_number_items == 0);
}

template <class TYPE>
inline
int List<TYPE>::getNumberItems() const
{
   return(d_number_items);
}

template <class TYPE>
inline
int List<TYPE>::size() const
{
   return(d_number_items);
}

template <class TYPE>
inline
TYPE& List<TYPE>::getFirstItem() const
{
   return(d_list_head->d_item);
}

template <class TYPE>
inline
TYPE& List<TYPE>::getLastItem() const
{
   return(d_list_tail->d_item);
}

template <class TYPE>
inline
ListIterator<TYPE> List<TYPE>::listStart() const
{
   return(ListIterator<TYPE>((List<TYPE> *) this, d_list_head));
}

template <class TYPE>
inline
ListIterator<TYPE> List<TYPE>::listEnd() const
{
   return(ListIterator<TYPE>((List<TYPE> *) this, d_list_tail));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

