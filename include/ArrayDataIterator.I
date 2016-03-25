//
// File:	ArrayDataIterator.I
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Iterator for array patch data types
//

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
ArrayDataIterator<DIM>::ArrayDataIterator()
{
}

template<int DIM> inline
ArrayDataIterator<DIM>::ArrayDataIterator(const hier::Box<DIM>& box)
:  d_index(box.lower()), d_box(box)
{
}

template<int DIM> inline
ArrayDataIterator<DIM>::ArrayDataIterator(const ArrayDataIterator<DIM>& iter)
:  d_index(iter.d_index), d_box(iter.d_box)
{
}

template<int DIM> inline
ArrayDataIterator<DIM>&
ArrayDataIterator<DIM>::operator=(const ArrayDataIterator<DIM>& iter)
{
   d_index = iter.d_index;
   d_box   = iter.d_box;
   return(*this);
}

template<int DIM> inline
ArrayDataIterator<DIM>::~ArrayDataIterator()
{
}

template<int DIM> inline
const hier::Index<DIM>& ArrayDataIterator<DIM>::operator*() const
{
   return(d_index);
}

template<int DIM> inline
const hier::Index<DIM>& ArrayDataIterator<DIM>::operator()() const
{
   return(d_index);
}

template<int DIM> inline
ArrayDataIterator<DIM>::operator bool() const
{
   bool retval = true;
   for (int i = 0; i < DIM; i++) {
      if (d_index(i) > d_box.upper(i)) {
         retval = false;
         break;
      }
   }

   return(retval);
}

#ifndef LACKS_BOOL_VOID_RESOLUTION
template<int DIM> inline
ArrayDataIterator<DIM>::operator const void*() const
{
   return(ArrayDataIterator<DIM>::operator bool() ? this : NULL);
}
#endif

template<int DIM> inline
bool ArrayDataIterator<DIM>::operator!() const
{
   return(!ArrayDataIterator<DIM>::operator bool());
}
   
template<int DIM> inline
void ArrayDataIterator<DIM>::operator++(int)
{
   d_index(0)++;
   for (int i = 0; i < DIM-1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i+1)++;
      } else {
         break;
      }
   }
}

template<int DIM> inline
bool ArrayDataIterator<DIM>::operator==(const ArrayDataIterator<DIM>& iter) const
{
   return(d_index == iter.d_index);
}

template<int DIM> inline
bool ArrayDataIterator<DIM>::operator!=(const ArrayDataIterator<DIM>& iter) const
{
   return(d_index != iter.d_index);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
