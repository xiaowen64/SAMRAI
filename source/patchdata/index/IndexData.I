//
// File:	IndexData.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	0.1
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
