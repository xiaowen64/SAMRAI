//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/index/IndexData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 2249 $
// Modified:	$LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description: hier::Patch data structure for irregular grid data
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexDataNode<DIM,TYPE,BOX_GEOMETRY>::IndexDataNode()
{
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexDataNode<DIM,TYPE,BOX_GEOMETRY>::IndexDataNode(
   const hier::Index<DIM>& index,
   const int offset,
   TYPE& t,
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *n,
   IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *p)
   :   d_index(index), d_offset(offset), d_item(&t), d_next(n), d_prev(p)
{
}


template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexDataNode<DIM,TYPE,BOX_GEOMETRY>::~IndexDataNode()
{
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::IndexIterator()
{
   d_index_data = ((IndexData<DIM,TYPE,BOX_GEOMETRY> *) NULL);
   d_node = ((IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *) NULL);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::IndexIterator(const IndexData<DIM,TYPE,BOX_GEOMETRY>& index_data)
{
   d_index_data = (IndexData<DIM,TYPE,BOX_GEOMETRY> *) &index_data;
   d_node = d_index_data -> d_list_head;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::IndexIterator(IndexData<DIM,TYPE,BOX_GEOMETRY>* index_data,
                                           IndexDataNode<DIM,TYPE,BOX_GEOMETRY> *node)
{
   d_index_data = index_data;
   d_node = node;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::IndexIterator(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iter)
{
   d_index_data = iter.d_index_data;
   d_node = iter.d_node;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>&
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator=(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& iter)
{
   d_index_data = iter.d_index_data;
   d_node = iter.d_node;
   return(*this);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::~IndexIterator()
{
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
TYPE& IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator*()
{
   return(*d_node->d_item);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
const TYPE& IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator*() const
{
   return(*d_node->d_item);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
TYPE& IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator()()
{
   return(*d_node->d_item);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
const TYPE& IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator()() const
{
   return(*d_node->d_item);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator bool() const
{
   return(d_node ? true : false);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator const void*() const
{
   return(d_node ? this : NULL);
}
#endif

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
bool IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator!() const
{
   return(d_node ? false : true);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
void IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator++(int)
{
   if (d_node) d_node = d_node->d_next;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
void IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator--(int)
{
   if (d_node) d_node = d_node->d_prev;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
void IndexIterator<DIM,TYPE,BOX_GEOMETRY>::rewindIterator()
{
   d_node = d_index_data->d_list_head;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
void IndexIterator<DIM,TYPE,BOX_GEOMETRY>::fastforwardIterator()
{
   d_node = d_index_data->d_list_tail;
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
bool IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator==(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& i) const
{
   return(d_node == i.d_node);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
bool IndexIterator<DIM,TYPE,BOX_GEOMETRY>::operator!=(const IndexIterator<DIM,TYPE,BOX_GEOMETRY>& i) const
{
   return(d_node != i.d_node);
}

template<int DIM, class TYPE, class BOX_GEOMETRY>
inline
IndexDataNode<DIM,TYPE,BOX_GEOMETRY> &IndexIterator<DIM,TYPE,BOX_GEOMETRY>::getNode()
{
   return(*d_node);
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
