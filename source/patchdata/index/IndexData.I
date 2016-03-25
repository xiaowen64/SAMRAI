//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/index/IndexData.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: hier::Patch data structure for irregular grid data
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM, class TYPE> inline
IndexIterator<DIM,TYPE>&
IndexIterator<DIM,TYPE>::operator=(
   const IndexIterator<DIM,TYPE>& iter)
{
   d_iterator = iter.d_iterator;
   return(*this);
}

template<int DIM, class TYPE> inline
IndexIterator<DIM,TYPE>::~IndexIterator()
{
}

template<int DIM, class TYPE> inline
TYPE& IndexIterator<DIM,TYPE>::operator*()
{
   return(*(d_index_data->getItem(d_iterator())));
}

template<int DIM, class TYPE> inline
TYPE& IndexIterator<DIM,TYPE>::operator()()
{
   return(*(d_index_data->getItem(d_iterator())));
}

template<int DIM, class TYPE> inline
TYPE& IndexIterator<DIM,TYPE>::getItem()
{
   return(*(d_index_data->getItem(d_iterator())));
}

template<int DIM, class TYPE> inline
const hier::Index<DIM>& IndexIterator<DIM,TYPE>::getIndex() const
{
   return(d_iterator());
}

template<int DIM, class TYPE> inline
IndexIterator<DIM,TYPE>::operator bool() const
{
   return((bool) d_iterator);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM, class TYPE> inline
IndexIterator<DIM,TYPE>::operator const void*() const
{
   return((const void *) d_iterator);
}
#endif

template<int DIM, class TYPE> inline
bool IndexIterator<DIM,TYPE>::operator!() const
{
   return(!IndexIterator<DIM,TYPE>::operator bool());
}

template<int DIM, class TYPE> inline
void IndexIterator<DIM,TYPE>::operator++(int)
{
   d_iterator++;
}

template<int DIM, class TYPE> inline
bool IndexIterator<DIM,TYPE>::operator==(
   const IndexIterator<DIM,TYPE>& iter) const
{
   return(d_iterator == iter.d_iterator);
}

template<int DIM, class TYPE> inline
bool IndexIterator<DIM,TYPE>::operator!=(
   const IndexIterator<DIM,TYPE>& iter) const
{
   return(d_iterator != iter.d_iterator);
}
#ifdef ENABLE_CONST_ITERATOR
// ConstIndexIterator

template<int DIM, class TYPE> inline
ConstIndexIterator<DIM,TYPE>&
ConstIndexIterator<DIM,TYPE>::operator=(
   const ConstIndexIterator<DIM,TYPE>& iter)
{
   d_iterator = iter.d_iterator;
   return(*this);
}

template<int DIM, class TYPE> inline
ConstIndexIterator<DIM,TYPE>&
ConstIndexIterator<DIM,TYPE>::operator=(
   const IndexIterator<DIM,TYPE>& iter)
{
   d_iterator = iter.d_iterator;
   return(*this);
}

template<int DIM, class TYPE> inline
ConstIndexIterator<DIM,TYPE>::~ConstIndexIterator()
{
}

template<int DIM, class TYPE> inline
const TYPE& ConstIndexIterator<DIM,TYPE>::operator*()
{
   return(*(d_index_data->getItem(d_iterator())));
}

template<int DIM, class TYPE> inline
const TYPE& ConstIndexIterator<DIM,TYPE>::operator()()
{
   return(*(d_index_data->getItem(d_iterator())));
}

template<int DIM, class TYPE> inline
const TYPE& ConstIndexIterator<DIM,TYPE>::getItem()
{
   return(*(d_index_data->getItem(d_iterator())));
}

template<int DIM, class TYPE> inline
const hier::Index<DIM>& ConstIndexIterator<DIM,TYPE>::getIndex() const
{
   return(d_iterator());
}

template<int DIM, class TYPE> inline
ConstIndexIterator<DIM,TYPE>::operator bool() const
{
   return((bool) d_iterator);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM, class TYPE> inline
ConstIndexIterator<DIM,TYPE>::operator const void*() const
{
   return((const void *) d_iterator);
}
#endif

template<int DIM, class TYPE> inline
bool ConstIndexIterator<DIM,TYPE>::operator!() const
{
   return(!ConstIndexIterator<DIM,TYPE>::operator bool());
}

template<int DIM, class TYPE> inline
void ConstIndexIterator<DIM,TYPE>::operator++(int)
{
   d_iterator++;
}

template<int DIM, class TYPE> inline
bool ConstIndexIterator<DIM,TYPE>::operator==(
   const ConstIndexIterator<DIM,TYPE>& iter) const
{
   return(d_iterator == iter.d_iterator);
}

template<int DIM, class TYPE> inline
bool ConstIndexIterator<DIM,TYPE>::operator!=(
   const ConstIndexIterator<DIM,TYPE>& iter) const
{
   return(d_iterator != iter.d_iterator);
}
#endif
#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
