/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexDataNode<TYPE, BOX_GEOMETRY>::IndexDataNode():
   d_index(tbox::Dimension::getInvalidDimension())
{
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexDataNode<TYPE, BOX_GEOMETRY>::IndexDataNode(
   const hier::Index& index,
   const int offset,
   TYPE& t,
   IndexDataNode<TYPE, BOX_GEOMETRY>* n,
   IndexDataNode<TYPE, BOX_GEOMETRY>* p):
   d_index(index),
   d_offset(offset),
   d_item(&t),
   d_next(n),
   d_prev(p)
{
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexDataNode<TYPE, BOX_GEOMETRY>::~IndexDataNode()
{
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator()
{
   d_index_data = ((IndexData<TYPE, BOX_GEOMETRY> *)NULL);
   d_node = ((IndexDataNode<TYPE, BOX_GEOMETRY> *)NULL);
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   const IndexData<TYPE, BOX_GEOMETRY>& index_data)
{
   d_index_data = (IndexData<TYPE, BOX_GEOMETRY> *) & index_data;
   d_node = d_index_data->d_list_head;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   IndexData<TYPE, BOX_GEOMETRY>* index_data,
   IndexDataNode<TYPE, BOX_GEOMETRY>* node)
{
   d_index_data = index_data;
   d_node = node;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   const IndexIterator<TYPE, BOX_GEOMETRY>& iter)
{
   d_index_data = iter.d_index_data;
   d_node = iter.d_node;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator = (
   const IndexIterator<TYPE, BOX_GEOMETRY>& iter)
{
   d_index_data = iter.d_index_data;
   d_node = iter.d_node;
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::~IndexIterator()
{
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
TYPE& IndexIterator<TYPE, BOX_GEOMETRY>::operator * ()
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const TYPE& IndexIterator<TYPE, BOX_GEOMETRY>::operator * () const
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
TYPE& IndexIterator<TYPE, BOX_GEOMETRY>::getItem()
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const hier::Index& IndexIterator<TYPE, BOX_GEOMETRY>::getIndex() const
{
   return d_node->d_index;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
TYPE& IndexIterator<TYPE, BOX_GEOMETRY>::operator () ()
{
   return getItem();
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const TYPE& IndexIterator<TYPE, BOX_GEOMETRY>::operator () () const
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::operator bool () const
{
   return d_node ? true : false;
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexIterator<TYPE, BOX_GEOMETRY>::operator const void * () const
{
   return d_node ? this : NULL;
}
#endif

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool IndexIterator<TYPE, BOX_GEOMETRY>::operator ! () const
{
   return d_node ? false : true;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void IndexIterator<TYPE, BOX_GEOMETRY>::operator ++ (
   int)
{
   if (d_node) d_node = d_node->d_next;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void IndexIterator<TYPE, BOX_GEOMETRY>::operator -- (
   int)
{
   if (d_node) d_node = d_node->d_prev;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void IndexIterator<TYPE, BOX_GEOMETRY>::rewindIterator()
{
   d_node = d_index_data->d_list_head;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void IndexIterator<TYPE, BOX_GEOMETRY>::fastforwardIterator()
{
   d_node = d_index_data->d_list_tail;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool IndexIterator<TYPE, BOX_GEOMETRY>::operator == (
   const IndexIterator<TYPE, BOX_GEOMETRY>& i) const
{
   return d_node == i.d_node;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool IndexIterator<TYPE, BOX_GEOMETRY>::operator != (
   const IndexIterator<TYPE, BOX_GEOMETRY>& i) const
{
   return d_node != i.d_node;
}

template<class TYPE, class BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
IndexDataNode<TYPE, BOX_GEOMETRY>& IndexIterator<TYPE, BOX_GEOMETRY>::getNode()
{
   return *d_node;
}

}
}
